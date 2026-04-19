import matplotlib.pyplot as plt
import pandas as pd
import re
import requests
from xml.etree import ElementTree as ET
import numpy as np

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
SEARCH_URL = BASE_URL + "esearch.fcgi"
FETCH_URL = BASE_URL + "efetch.fcgi"
EMAIL = "a.ekergard@proton.com"

def search_pubmed(query: str, max_results: int = 1000) -> list[str]:
    params = {"db": "pubmed", "term": query, "retmax": max_results, "retmode": "xml", "email": EMAIL}
    response = requests.get(SEARCH_URL, params=params)
    response.raise_for_status()
    root = ET.fromstring(response.text)
    return [id_elem.text for id_elem in root.findall(".//IdList/Id")]

def fetch_article_abstracts(pmids: list[str]) -> list[str]:
    if not pmids:
        return []
    params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml", "rettype": "abstract", "email": EMAIL}
    response = requests.get(FETCH_URL, params=params)
    response.raise_for_status()
    root = ET.fromstring(response.text)
    abstracts = []
    for article in root.findall(".//PubmedArticle"):
        abstract_text = " ".join(text_elem.text.strip() for text_elem in article.findall(".//AbstractText") if text_elem.text)
        if abstract_text:
            abstracts.append(abstract_text)
    return abstracts

# --- Regex-generering ---
def generate_regex_patterns(text: str) -> list[str]:
    words = text.split()
    patterns = []

    if len(words) == 1:
        return [
            rf'\b{re.escape(words[0])}\b',
            rf'\b{re.escape(words[0][0])}\.\s*{re.escape(words[0][1:])}\b',
            rf'\b{re.escape(words[0][0])}\s*{re.escape(words[0][1:])}\b'
        ]

    first_word = words[0]
    rest = " ".join(words[1:])
    patterns = [
        rf'\b{re.escape(text)}\b',
        rf'\b{re.escape(first_word[0])}\.\s*{re.escape(rest)}\b',
        rf'\b{re.escape(first_word[0])}\s*{re.escape(rest)}\b',
        rf'\b{re.escape(first_word[0])}\.{re.escape(rest)}\b',
    ]

    for word in words:
        patterns.append(rf'\b{re.escape(word)}\b')

    return patterns

def find_disease(abstracts: list[str], microbe_patterns: list[str], fishs: list[str], microbe_name: str) -> dict[str, dict[str, int]]:
    “”"Analyze abstracts to count the number of mentions of the microorganism in connection with each fish, summed across all pattern matches.
    Args:
        abstracts (list[str]): List of abstracts to analyze
        microbe_patterns (list[str]): List of regex patterns for microorganisms
        fishs (list[str]): List of fish species
        microbe_name (str): Name of the microorganism/disease used as a key in the results dictionary

    Returns:
        dict[str, dict[str, int]]: A dictionary containing a single microorganism name and counts per fish, as well as the total
    “”"
    


counts = {fish: 0 for fish in fishs}
    counts["total"] = 0

    if not abstracts:
        return {}

    combined_pattern = "|".join(sorted(set(microbe_patterns), key=len, reverse=True))

    for abstract in abstracts:
        for match in re.finditer(combined_pattern, abstract, re.IGNORECASE):
            full_match = match.group(0).strip()
            counts["total"] += 1
            for fish in fishs:
                if re.search(
                    rf"{re.escape(full_match)}.*{re.escape(fish)}|{re.escape(fish)}.*{re.escape(full_match)}",
                    abstract,
                    re.IGNORECASE,
                ):
                    counts[fish] += 1

    return {microbe_name: counts} if counts["total"] > 0 else {}

# --- Visualisering ---
def plot_microbe_mentions(data: dict[str, dict[str, int]], fishs: list[str]) -> None:
    rows = []
    for bacteria, counts in data.items():
        total_mentions = counts.get("total", 0)
        for fish in fishs:
            rows.append({
                "bacteria": bacteria,
                "fish": fish,
                "count": counts.get(fish, 0),
                "total": total_mentions
            })
    df = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_title("Mentions of Microbes in Abstracts", fontsize=16)
    ax.set_xlabel("Microbe", fontsize=14)
    ax.set_ylabel("Count", fontsize=14)

    colors = plt.cm.tab10(np.linspace(0, 1, len(fishs) + 1))


    bacteria_order = sorted(df["bacteria"].unique())
    total_df = df[df["fish"] == fishs[0]].drop_duplicates(subset="bacteria").set_index("bacteria").reindex(bacteria_order).reset_index()


    ax.bar(
        total_df["bacteria"],
        total_df["total"],
        color=colors[0],
        label="Total Mentions",
        alpha=0.8,
        width=0.6
    )


    x = np.arange(len(bacteria_order))
    width = 0.6 / (len(fishs) + 1)  

    for i, fish in enumerate(fishs):
        fish_df = df[df["fish"] == fish].set_index("bacteria").reindex(bacteria_order).reset_index()
        ax.bar(
            x + (i + 1) * width,
            fish_df["count"],
            color=colors[i + 1],
            label=f"Mentions with {fish}",
            alpha=0.8,
            width=width
        )

    ax.set_xticks(x + width * (len(fishs) / 2))
    ax.set_xticklabels(bacteria_order, rotation=45, ha="right", fontsize=10)
    ax.legend(fontsize=12, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.yticks(fontsize=12)
    plt.tight_layout()
    fig.savefig("microbe_mentions.png", dpi=300, bbox_inches="tight")
    plt.show()

def main():
    user_input = input("Please provide disease of interest: ").strip()
    if not user_input:
        raise ValueError("Disease input cannot be empty.")

    fishs = ['tilapia', 'salmon', 'shrimp']
    query = f"{user_input} aquaculture"

    pmids = search_pubmed(query, max_results=100)
    abstracts = fetch_article_abstracts(pmids)
    microbe_patterns = generate_regex_patterns(user_input)
    data = find_disease(abstracts, microbe_patterns, fishs, user_input)

    if not data:
        print("No mentions found for the provided disease.")
    else:
        plot_microbe_mentions(data, fishs)

if __name__ == "__main__":
    main()
