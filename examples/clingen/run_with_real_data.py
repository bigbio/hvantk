"""
Run ClinGen geneset extraction with real (downloaded) data.

This script uses the downloaded ClinGen Gene-Disease Validity CSV file
instead of the test data.
"""

import logging
import json
from pathlib import Path
from typing import Dict, Set

import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Paths
SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "data"
OUTPUT_DIR = SCRIPT_DIR / "results"
INPUT_CSV = DATA_DIR / "clingen_gene_disease.csv"


def main():
    """Main workflow with real ClinGen data."""
    import hail as hl
    from hvantk.tables.table_builders import create_clingen_gene_disease_tb
    from hvantk.data.clingen_streamer import ClinGenStreamer

    # Initialize Hail
    logger.info("Initializing Hail...")
    hl.init(quiet=True)

    # Ensure output directory exists
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Paths
    input_csv = str(INPUT_CSV)
    output_ht = str(OUTPUT_DIR / "clingen_full.ht")

    if not INPUT_CSV.exists():
        print(f"ERROR: ClinGen data file not found at {INPUT_CSV}")
        print("Please download it first with:")
        print('  curl -L -k "https://search.clinicalgenome.org/kb/gene-validity/download" -o examples/clingen/data/clingen_gene_disease.csv')
        return

    # Step 1: Build the ClinGen Hail Table
    print("\n" + "=" * 70)
    print("Step 1: Building ClinGen Hail Table from real data")
    print("=" * 70)
    logger.info(f"Building ClinGen table from: {input_csv}")
    create_clingen_gene_disease_tb(
        input_path=input_csv,
        output_path=output_ht,
        key_by="gene_disease",
        overwrite=True,
    )
    logger.info(f"ClinGen table saved to: {output_ht}")

    # Step 2: Extract genesets per disease using ClinGenStreamer
    print("\n" + "=" * 70)
    print("Step 2: Extracting genesets per disease")
    print("=" * 70)
    streamer = ClinGenStreamer(output_ht, init_hail=False)
    streamer.setup()

    genesets_per_disease = streamer.get_geneset_per_disease()
    logger.info(f"Found {len(genesets_per_disease)} unique diseases")

    # Show sample of diseases
    print(f"\nFound {len(genesets_per_disease)} unique diseases")
    print("\nSample of genesets (first 15 diseases):")
    print("-" * 50)
    for i, (disease, genes) in enumerate(sorted(genesets_per_disease.items())[:15]):
        genes_str = ", ".join(sorted(genes)[:5])
        if len(genes) > 5:
            genes_str += f", ... (+{len(genes)-5} more)"
        print(f"  {disease}:")
        print(f"    Genes ({len(genes)}): {genes_str}")

    # Save genesets to JSON
    genesets_json = {k: sorted(v) for k, v in genesets_per_disease.items()}
    genesets_path = OUTPUT_DIR / "genesets_per_disease_full.json"
    with open(genesets_path, "w") as f:
        json.dump(genesets_json, f, indent=2)
    print(f"\nAll genesets saved to: {genesets_path}")

    # Step 3: Summarize genes per disease category
    print("\n" + "=" * 70)
    print("Step 3: Summary of genes per disease category")
    print("=" * 70)

    # Define disease categories
    categories = {
        "Cancer/Tumor": ["cancer", "carcinoma", "tumor", "neoplasm", "melanoma", "leukemia", "lymphoma"],
        "Cardiovascular": ["cardiac", "heart", "arrhythmia", "cardiomyopathy", "aortic", "vascular"],
        "Neurological": ["neurological", "brain", "epilepsy", "seizure", "ataxia", "neuropathy", "encephalopathy"],
        "Metabolic": ["metabolic", "diabetes", "hypercholesterolemia", "lipid", "obesity"],
        "Developmental": ["developmental", "intellectual disability", "autism", "congenital"],
        "Immunological": ["immune", "immunodeficiency", "autoimmune", "inflammatory"],
        "Hereditary Syndromes": ["syndrome", "hereditary", "familial"],
        "Eye/Vision": ["retinal", "retinopathy", "macular", "blindness", "optic", "eye"],
        "Hearing": ["hearing", "deafness", "auditory"],
        "Kidney/Renal": ["kidney", "renal", "nephropathy", "nephrotic"],
    }

    category_genesets = streamer.aggregate_by_disease_category(categories)

    # Build summary DataFrame
    summary_data = []
    for category, genes in sorted(category_genesets.items(), key=lambda x: -len(x[1])):
        summary_data.append({
            "category": category,
            "n_genes": len(genes),
            "sample_genes": ", ".join(sorted(genes)[:10]) + ("..." if len(genes) > 10 else ""),
        })

    df = pd.DataFrame(summary_data)
    print("\n" + df.to_string(index=False))

    # Save summary
    summary_path = OUTPUT_DIR / "category_summary_full.tsv"
    df.to_csv(summary_path, sep="\t", index=False)
    print(f"\nSummary saved to: {summary_path}")

    # Step 4: Classification distribution
    print("\n" + "=" * 70)
    print("Step 4: Classification distribution")
    print("=" * 70)
    classification_df = streamer.classification_summary()
    print("\n" + classification_df.to_string(index=False))

    # Save classification distribution
    class_path = OUTPUT_DIR / "classification_distribution_full.tsv"
    classification_df.to_csv(class_path, sep="\t", index=False)
    print(f"\nClassification distribution saved to: {class_path}")

    # Step 5: Compute overall statistics
    print("\n" + "=" * 70)
    print("Step 5: Dataset Statistics")
    print("=" * 70)
    stats = streamer.compute_stats()

    print(f"\nTotal gene-disease associations: {stats['total_associations']:,}")
    print(f"Unique genes: {stats['unique_genes']:,}")
    print(f"Unique diseases: {stats['unique_diseases']:,}")
    print(f"Last updated: {stats['last_update']}")

    print("\nTop 10 genes by number of associated diseases:")
    for gene, n_diseases in stats['top_genes_by_diseases']:
        print(f"  {gene}: {n_diseases} diseases")

    # Final summary
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"Total unique diseases: {len(genesets_per_disease)}")
    total_genes = set()
    for genes in genesets_per_disease.values():
        total_genes.update(genes)
    print(f"Total unique genes: {len(total_genes)}")
    print(f"\nOutput files:")
    print(f"  - {genesets_path}")
    print(f"  - {summary_path}")
    print(f"  - {class_path}")

    # Stop Hail
    hl.stop()


if __name__ == "__main__":
    main()
