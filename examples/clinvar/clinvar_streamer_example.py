# Example usage of Clinvar Data Streamer
# Demonstrates how to use the streaming architecture for different scenarios

import logging
from hvantk.core.streamers.clinvar import (
    ClinvarDataStreamer,
    create_clinvar_training_set_streamer,
)
import hail as hl
from pathlib import Path

# General ClinVar path variable (adjust as needed for your environment)
CLINVAR_TEST_VCF = Path(
    "../hvantk/tests/testdata/raw/clinvar/clinvar_20220403_chr20.vcf.bgz"
)
# Primary path alias used across examples (can be reassigned if needed)
clinvar_path = CLINVAR_TEST_VCF

logger = logging.getLogger(__name__)


def example_basic_usage():
    """Example of basic Clinvar training set generation"""
    print("=== Basic Clinvar Training Set Generation ===")

    # Sample CHD genes
    gene_set = {"GATA4", "NKX2-5", "TBX5", "NOTCH1", "CHD7"}

    # Create processor
    processor = create_clinvar_training_set_streamer(
        clinvar_path=str(clinvar_path),
        output_dir="./data/training_set",
        gene_set=gene_set,
    )

    # Process data
    result = processor.process()

    if result:
        print(f"Generated training set with {result.count()} variants")
    else:
        print("No training set generated")


def example_custom_streamer():
    """Example of using individual streamers with custom processing"""
    print("=== Custom Streamer Usage ===")

    # Create custom CHD gene set
    custom_gene_set = {
        "GATA4",
        "NKX2-5",
        "TBX5",
        "NOTCH1",
        "CHD7",
        "TBX1",
        "MYH6",
        "ACTC1",
        "MYH7",
        "TNNT2",
        "SCN5A",
        "KCNQ1",
        "KCNH2",
        "RYR2",
        "PKP2",
    }

    # Create individual streamer
    streamer = ClinvarDataStreamer(
        clinvar_path=str(clinvar_path),
        gene_set=custom_gene_set,
        chunk_size=5000,  # Smaller chunks for demo
    )

    # Manual processing
    streamer.setup()

    try:
        chunk_count = 0
        total_variants = 0

        for chunk in streamer.stream():
            chunk_count += 1
            chunk_size = chunk.count()
            total_variants += chunk_size
            print(f"Processed chunk {chunk_count}: {chunk_size} variants")

            # Could do additional processing here per chunk
            # e.g., apply filters, compute statistics, etc.

        print(f"Total processed: {total_variants} variants in {chunk_count} chunks")

    finally:
        streamer.teardown()


def example_with_validation():
    """Example with data validation and quality checks"""
    print("=== Training Set Generation with Validation ===")

    processor = create_clinvar_training_set_streamer(
        clinvar_path=str(clinvar_path), output_dir="./data/training_set"
    )

    # Process with validation
    result = processor.process()

    if result:
        # Validate results
        total_count = result.count()
        tp_count = result.filter(result.rf_label == "TP").count()
        tn_count = result.filter(result.rf_label == "TN").count()

        print(f"Training Set Validation:")
        print(f"  Total variants: {total_count}")
        print(f"  True Positives: {tp_count} ({tp_count/total_count*100:.1f}%)")
        print(f"  True Negatives: {tn_count} ({tn_count/total_count*100:.1f}%)")

        # Check for balanced dataset
        if abs(tp_count - tn_count) / total_count > 0.3:
            print("  WARNING: Dataset is highly imbalanced!")
        else:
            print("  Dataset balance looks reasonable")

        # Show gene distribution
        gene_counts = result.group_by(result.gene).aggregate(
            n_variants=hl.agg.count(),
            n_tp=hl.agg.count_where(result.rf_label == "TP"),
            n_tn=hl.agg.count_where(result.rf_label == "TN"),
            tp_fraction=hl.agg.fraction(result.rf_label == "TP"),
            tn_fraction=hl.agg.fraction(result.rf_label == "TN"),
        )
        # Order by the per-gene variant count (field name avoids method collision)
        top_genes = gene_counts.order_by(hl.desc(gene_counts.n_variants)).take(10)

        print("  Top 10 genes by variant count (with TP/TN stats):")
        for gene_row in top_genes:
            tp_pct = (
                f"{gene_row.tp_fraction:.2%}"
                if gene_row.tp_fraction is not None
                else "NA"
            )
            tn_pct = (
                f"{gene_row.tn_fraction:.2%}"
                if gene_row.tn_fraction is not None
                else "NA"
            )
            print(
                f"    {gene_row.gene}: total={gene_row.n_variants} TP={gene_row.n_tp} TN={gene_row.n_tn} "
                f"TP%={tp_pct} TN%={tn_pct}"
            )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    print("Clinvar Data Streamer Examples")
    print("=" * 50)

    # Note: These examples assume you have Clinvar data available
    # Uncomment the examples you want to run:

    # example_basic_usage()
    # example_custom_streamer()
    example_with_validation()

    print("\nTo run these examples, ensure you have:")
    print("1. Clinvar VCF file available")
    print("2. Proper Hail environment setup")
    print("3. Sufficient memory for processing")
    print("\nUncomment the example functions above to run them.")
