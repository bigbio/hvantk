"""
Run ClinGen geneset extraction with MONDO ontology-based categorization.

This script demonstrates how to use the MONDO disease ontology to properly
categorize diseases based on their ontological relationships (is_a hierarchy),
rather than simple keyword matching.

Usage:
    python run_ontology_categorization.py

Requirements:
    - Hail installed and configured
    - ClinGen raw CSV file (downloaded)
    - MONDO OBO file (downloaded)
"""

import logging
import json
from pathlib import Path
from typing import Dict, Set

import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Paths
SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "data"
OUTPUT_DIR = SCRIPT_DIR / "results"
INPUT_CSV = DATA_DIR / "clingen_gene_disease.csv"
MONDO_OBO = DATA_DIR / "mondo.obo"


def ensure_data_files():
    """Ensure ClinGen and MONDO data files are available."""
    import subprocess

    DATA_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_CSV.exists():
        print("Downloading ClinGen data...")
        subprocess.run(
            [
                "curl",
                "-L",
                "-k",
                "https://search.clinicalgenome.org/kb/gene-validity/download",
                "-o",
                str(INPUT_CSV),
            ],
            check=True,
        )

    if not MONDO_OBO.exists():
        print("Downloading MONDO ontology...")
        subprocess.run(
            [
                "curl",
                "-L",
                "-k",
                "https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.obo",
                "-o",
                str(MONDO_OBO),
            ],
            check=True,
        )


def main():
    """Main workflow with ontology-based categorization."""
    import hail as hl
    from hvantk.skills.clingen.streamer import ClinGenStreamer
    from hvantk.core.utils.mondo_parser import MONDO_DISEASE_CATEGORIES

    # Ensure data files exist
    ensure_data_files()

    # Initialize Hail
    logger.info("Initializing Hail...")
    hl.init(quiet=True)

    # Ensure output directory exists
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Paths
    output_ht = str(OUTPUT_DIR / "clingen_full.ht")
    mondo_obo = str(MONDO_OBO)

    # Step 1: Ensure the ClinGen Hail Table exists.
    # Builds go through the plugin system; this example focuses on the
    # ontology-based categorization that consumes a built artifact.
    print("\n" + "=" * 70)
    print("Step 1: ClinGen Hail Table")
    print("=" * 70)

    if not Path(output_ht).exists():
        print(f"ERROR: ClinGen Hail Table not found at {output_ht}")
        print("Build it first with:")
        print(
            f"  hvantk reprocess clingen:gene_disease "
            f"--raw-dir {DATA_DIR} --output {output_ht} --skip-download"
        )
        return
    logger.info(f"Using existing ClinGen table: {output_ht}")

    # Step 2: Initialize ClinGenStreamer
    print("\n" + "=" * 70)
    print("Step 2: Loading ClinGen data and MONDO ontology")
    print("=" * 70)

    streamer = ClinGenStreamer(output_ht, init_hail=False)
    streamer.setup()

    # Step 3: Ontology-based categorization
    print("\n" + "=" * 70)
    print("Step 3: Ontology-based disease categorization")
    print("=" * 70)

    print("\nUsing MONDO ontology hierarchy to categorize diseases...")
    print("This maps each disease to its ontological ancestors in MONDO.")
    print("\nAvailable categories:")
    for mondo_id, name in sorted(MONDO_DISEASE_CATEGORIES.items(), key=lambda x: x[1]):
        print(f"  {mondo_id}: {name}")

    # Get categorization results
    results = streamer.categorize_by_ontology(
        ontology=mondo_obo,
        min_classification="Limited",  # Include Limited and above
    )

    # Get summary DataFrame
    summary_df = streamer.categorize_by_ontology_summary(
        ontology=mondo_obo,
        min_classification="Limited",
    )

    print("\n" + "=" * 70)
    print("Ontology-Based Categorization Results")
    print("=" * 70)
    print("\n" + summary_df.to_string(index=False))

    # Save summary
    summary_path = OUTPUT_DIR / "ontology_category_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"\nSummary saved to: {summary_path}")

    # Step 4: Compare with keyword-based approach
    print("\n" + "=" * 70)
    print("Step 4: Comparison with keyword-based categorization")
    print("=" * 70)

    # Keyword-based categories (from previous approach)
    keyword_categories = {
        "Cancer/Tumor": ["cancer", "carcinoma", "tumor", "neoplasm", "melanoma"],
        "Cardiovascular": ["cardiac", "heart", "arrhythmia", "cardiomyopathy"],
        "Neurological": ["neurological", "brain", "epilepsy", "seizure", "ataxia"],
        "Metabolic": ["metabolic", "diabetes", "hypercholesterolemia"],
    }

    keyword_results = streamer.aggregate_by_disease_category(keyword_categories)

    print("\nKeyword-based vs Ontology-based gene counts:")
    print("-" * 60)
    print(f"{'Category':<25} {'Keyword':<12} {'Ontology':<12} {'Difference':<12}")
    print("-" * 60)

    # Map ontology categories to comparable keyword categories
    category_mapping = {
        "Cancer/Tumor": ["cancer", "neoplasm"],
        "Cardiovascular": ["cardiovascular disease", "cardiogenetic disease"],
        "Neurological": [
            "nervous system disease",
            "hereditary neurological disease",
            "neurodevelopmental disorder",
        ],
        "Metabolic": ["metabolic disease"],
    }

    for kw_cat, onto_cats in category_mapping.items():
        kw_count = len(keyword_results.get(kw_cat, set()))

        # Combine ontology categories
        onto_genes = set()
        for onto_cat in onto_cats:
            if onto_cat in results:
                onto_genes.update(results[onto_cat]["genes"])
        onto_count = len(onto_genes)

        diff = onto_count - kw_count
        diff_str = f"+{diff}" if diff > 0 else str(diff)
        print(f"{kw_cat:<25} {kw_count:<12} {onto_count:<12} {diff_str:<12}")

    # Step 5: Detailed results for specific categories
    print("\n" + "=" * 70)
    print("Step 5: Detailed results for selected categories")
    print("=" * 70)

    selected_categories = [
        "cardiovascular disease",
        "neoplasm",
        "nervous system disease",
        "eye disease",
        "hereditary disease",
    ]

    for category in selected_categories:
        if category in results:
            data = results[category]
            print(f"\n{category.upper()}")
            print("-" * 50)
            print(f"  Genes: {len(data['genes'])}")
            print(f"  Diseases: {len(data['diseases'])}")
            print(f"  Sample genes: {', '.join(sorted(data['genes'])[:15])}")
            if len(data["genes"]) > 15:
                print(f"    ... and {len(data['genes']) - 15} more")

    # Step 6: Save detailed results
    print("\n" + "=" * 70)
    print("Step 6: Saving detailed results")
    print("=" * 70)

    # Save genes per category
    genes_by_category = {cat: sorted(data["genes"]) for cat, data in results.items()}
    genes_path = OUTPUT_DIR / "genes_by_ontology_category.json"
    with open(genes_path, "w") as f:
        json.dump(genes_by_category, f, indent=2)
    print(f"Genes by category saved to: {genes_path}")

    # Save diseases per category
    diseases_by_category = {
        cat: sorted(data["diseases"]) for cat, data in results.items()
    }
    diseases_path = OUTPUT_DIR / "diseases_by_ontology_category.json"
    with open(diseases_path, "w") as f:
        json.dump(diseases_by_category, f, indent=2)
    print(f"Diseases by category saved to: {diseases_path}")

    # Final summary
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)

    total_categorized_genes = set()
    for data in results.values():
        if data != results.get("uncategorized", {}):
            total_categorized_genes.update(data["genes"])

    uncategorized_genes = results.get("uncategorized", {}).get("genes", set())

    print(f"Total genes categorized: {len(total_categorized_genes)}")
    print(f"Genes in 'uncategorized': {len(uncategorized_genes)}")
    print(f"Categories with genes: {sum(1 for d in results.values() if d['genes'])}")

    print(f"\nOutput files:")
    print(f"  - {summary_path}")
    print(f"  - {genes_path}")
    print(f"  - {diseases_path}")

    # Stop Hail
    hl.stop()


if __name__ == "__main__":
    main()
