#!/usr/bin/env python
"""
Prepare Gene Set Collections for PSROC Multi-Group Analysis

This script demonstrates how to extract named gene set collections from
ClinGen disease data for use with PSROC's --gene-sets option.

The script uses hvantk's streamer layer (Layer 2) to transform a ClinGen
Hail Table into a gene set collection JSON file that PSROC can consume.

Usage:
    # From a pre-built ClinGen Hail Table
    python examples/psroc/prepare_gene_sets.py \
        --clingen-ht /data/tables/clingen.ht \
        --output /data/gene_sets/disease_categories.json

    # From a ClinGen CSV (builds the table first)
    python examples/psroc/prepare_gene_sets.py \
        --clingen-csv /data/clingen/gene_curation_list.csv \
        --output /data/gene_sets/disease_categories.json

    # Then run PSROC with the gene sets:
    hvantk psroc \
        --gene-sets /data/gene_sets/disease_categories.json \
        --clinvar-ht /data/tables/clinvar_grch38.ht \
        --dbnsfp-ht /data/tables/dbnsfp_grch38.ht \
        --scores "CADD_phred,REVEL_score,MetaLR_score" \
        --output-dir /results/psroc_by_disease
"""

import argparse
import json
import sys
from pathlib import Path


# Default disease category definitions (keyword-based matching)
DEFAULT_CATEGORIES = {
    "cardiac": [
        "cardiomyopathy",
        "arrhythmia",
        "long_qt",
        "brugada",
        "channelopathy",
    ],
    "neurological": [
        "epilepsy",
        "neuropathy",
        "ataxia",
        "encephalopathy",
        "neurodegeneration",
    ],
    "cancer": [
        "cancer",
        "tumor",
        "neoplasm",
        "lynch",
        "polyposis",
    ],
    "metabolic": [
        "metabolic",
        "lysosomal",
        "glycogen",
        "mitochondrial",
    ],
    "connective_tissue": [
        "marfan",
        "ehlers",
        "osteogenesis",
        "connective",
    ],
}


def prepare_from_hail_table(clingen_ht_path: str, categories: dict) -> dict:
    """Extract gene sets from a ClinGen Hail Table using keyword matching.

    Parameters
    ----------
    clingen_ht_path : str
        Path to pre-built ClinGen Hail Table.
    categories : dict
        Mapping of category name to list of search terms.

    Returns
    -------
    dict
        Mapping of category name to set of gene symbols.
    """
    from hvantk.data.clingen_streamer import ClinGenStreamer

    streamer = ClinGenStreamer(table_path=clingen_ht_path)
    streamer.setup()

    gene_sets = streamer.aggregate_by_disease_category(
        categories=categories,
        min_classification="Definitive",
    )

    return gene_sets


def prepare_from_csv(clingen_csv_path: str, categories: dict, tmp_dir: str) -> dict:
    """Build a ClinGen Hail Table from CSV, then extract gene sets.

    Parameters
    ----------
    clingen_csv_path : str
        Path to ClinGen gene curation list CSV.
    categories : dict
        Mapping of category name to list of search terms.
    tmp_dir : str
        Directory for the intermediate Hail Table.

    Returns
    -------
    dict
        Mapping of category name to set of gene symbols.
    """
    from hvantk.core.hail_context import init_hail
    from hvantk.tables.table_builders import create_clingen_gene_disease_tb

    init_hail(quiet=True)

    ht_path = str(Path(tmp_dir) / "clingen.ht")
    create_clingen_gene_disease_tb(
        input_path=clingen_csv_path,
        output_path=ht_path,
        overwrite=True,
    )
    print(f"  Built ClinGen HT: {ht_path}")

    return prepare_from_hail_table(ht_path, categories)


def save_gene_sets(gene_sets: dict, output_path: str) -> None:
    """Save gene sets as a GeneSetCollection JSON file.

    Parameters
    ----------
    gene_sets : dict
        Mapping of category name to set of gene symbols.
    output_path : str
        Path to write the JSON file.
    """
    from hvantk.utils.gene_sets import load_gene_sets_from_dict

    # Convert sets to lists for the loader
    gene_sets_lists = {name: sorted(genes) for name, genes in gene_sets.items()}
    collection = load_gene_sets_from_dict(gene_sets_lists, source="clingen_disease")
    collection.save(output_path)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare gene set collections for PSROC multi-group analysis"
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--clingen-ht",
        type=str,
        help="Path to pre-built ClinGen Hail Table",
    )
    source.add_argument(
        "--clingen-csv",
        type=str,
        help="Path to ClinGen gene curation list CSV (will build HT first)",
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output path for gene set collection JSON",
    )
    parser.add_argument(
        "--categories-json",
        type=str,
        default=None,
        help="Custom categories JSON file (default: built-in disease categories)",
    )
    parser.add_argument(
        "--tmp-dir",
        type=str,
        default="/tmp/psroc_gene_sets",
        help="Temporary directory for intermediate tables",
    )

    args = parser.parse_args()

    # Load categories
    if args.categories_json:
        with open(args.categories_json) as f:
            categories = json.load(f)
        print(f"Loaded {len(categories)} custom categories from {args.categories_json}")
    else:
        categories = DEFAULT_CATEGORIES
        print(f"Using {len(categories)} default disease categories")

    # Extract gene sets
    print("\nExtracting gene sets from ClinGen data...")
    if args.clingen_ht:
        gene_sets = prepare_from_hail_table(args.clingen_ht, categories)
    else:
        gene_sets = prepare_from_csv(args.clingen_csv, categories, args.tmp_dir)

    # Report
    print("\nGene sets extracted:")
    for name, genes in sorted(gene_sets.items()):
        print(f"  {name}: {len(genes)} genes")

    # Filter out empty groups
    non_empty = {name: genes for name, genes in gene_sets.items() if genes}
    if len(non_empty) < len(gene_sets):
        empty = set(gene_sets) - set(non_empty)
        print(f"\nWarning: {len(empty)} empty group(s) removed: {', '.join(empty)}")

    if not non_empty:
        print("\nError: No gene sets with genes found. Check categories and data.")
        return 1

    # Save
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    save_gene_sets(non_empty, str(output_path))

    print(f"\nGene set collection saved to: {output_path}")
    print(f"\nNext step — run PSROC:")
    print(f"  hvantk psroc \\")
    print(f"    --gene-sets {output_path} \\")
    print(f"    --clinvar-ht /data/tables/clinvar_grch38.ht \\")
    print(f"    --dbnsfp-ht /data/tables/dbnsfp_grch38.ht \\")
    print(f"    --scores \"CADD_phred,REVEL_score,MetaLR_score\" \\")
    print(f"    --output-dir /results/psroc_by_disease")

    return 0


if __name__ == "__main__":
    sys.exit(main())
