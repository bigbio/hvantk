#!/usr/bin/env python
"""
PSROC Example: Running Prediction Score ROC Analysis

This example demonstrates how to use the PSROC pipeline to evaluate
variant pathogenicity prediction scores against ClinVar truth labels
using synthetic test data.

The example covers:
1. Building Hail Tables from synthetic TSV data
2. Running the PSROC pipeline via Python API
3. Interpreting results and generated plots

Requirements:
    - hvantk installed with Hail: poetry install
    - Run from the repository root

Usage:
    python examples/psroc/run_psroc_example.py [--output-dir DIR]

Output:
    - ROC curve plots for each prediction score
    - AUC comparison bar chart
    - Missingness summary
    - Summary dashboard
    - Metrics JSON file
"""

import argparse
import sys
from pathlib import Path


def main(output_dir: str = "/tmp/psroc_example") -> int:
    """Run the PSROC example with synthetic data."""

    print("=" * 70)
    print("PSROC Example: Prediction Score ROC Analysis")
    print("=" * 70)
    print()

    # Locate test data
    repo_root = Path(__file__).parent.parent.parent
    testdata_dir = repo_root / "hvantk" / "tests" / "testdata" / "psroc"

    clinvar_tsv = testdata_dir / "synthetic_clinvar.tsv"
    dbnsfp_tsv = testdata_dir / "synthetic_dbnsfp.tsv"

    if not clinvar_tsv.exists() or not dbnsfp_tsv.exists():
        print("ERROR: Synthetic test data not found.")
        print("Please run: python hvantk/tests/testdata/psroc/generate_synthetic_data.py")
        return 1

    print(f"Test data directory: {testdata_dir}")
    print(f"Output directory: {output_dir}")
    print()

    # Import hvantk modules (requires Hail)
    print("Loading hvantk modules...")
    try:
        from hvantk.core.hail_context import init_hail
        from hvantk.psroc import PSROCConfig, PSROCPipeline
    except ImportError as e:
        print(f"ERROR: Failed to import hvantk modules: {e}")
        print("Make sure hvantk is installed: poetry install")
        return 1

    # Initialize Hail
    print("Initializing Hail...")
    init_hail(quiet=True)

    # Import the data generation module
    sys.path.insert(0, str(testdata_dir))
    from generate_synthetic_data import build_hail_tables

    # Build Hail Tables from TSV files
    print()
    print("-" * 70)
    print("Step 1: Building Hail Tables from synthetic TSV data")
    print("-" * 70)

    ht_dir = Path(output_dir) / "tables"
    clinvar_ht, dbnsfp_ht = build_hail_tables(
        str(clinvar_tsv),
        str(dbnsfp_tsv),
        str(ht_dir),
    )

    # Configure PSROC pipeline
    print()
    print("-" * 70)
    print("Step 2: Configuring PSROC pipeline")
    print("-" * 70)

    config = PSROCConfig(
        genes=["BRCA1", "BRCA2", "TP53"],
        clinvar_ht=clinvar_ht,
        dbnsfp_ht=dbnsfp_ht,
        scores=["CADD_phred", "REVEL_score", "MetaLR_score", "VEST4_score"],
        output_dir=str(Path(output_dir) / "results"),
        output_prefix="psroc_example",
        max_missingness=0.3,
        threshold_method="youden",
        min_stars=0,  # Include all variants for this example
        generate_plots=True,
        export_tsv=True,
        overwrite=True,
    )

    print(f"Genes: {config.genes}")
    print(f"Scores: {config.scores}")
    print(f"Max missingness: {config.max_missingness}")
    print(f"Threshold method: {config.threshold_method}")

    # Run pipeline
    print()
    print("-" * 70)
    print("Step 3: Running PSROC pipeline")
    print("-" * 70)

    pipeline = PSROCPipeline(config)
    result = pipeline.run()

    # Display results
    print()
    print("-" * 70)
    print("Step 4: Results Summary")
    print("-" * 70)
    print()

    print(f"Total variants analyzed: {result.n_total}")
    print(f"  - Pathogenic: {result.n_pathogenic}")
    print(f"  - Benign: {result.n_benign}")
    print(f"  - Excluded (VUS): {result.n_excluded}")
    print()

    print("Scores included in analysis:")
    for name in result.scores_included:
        roc = result.metrics[name]
        miss = result.missingness[name]
        print(f"  - {name}:")
        print(f"      AUC: {roc.auc:.3f}")
        print(f"      Optimal threshold: {roc.optimal_threshold:.3f}")
        print(f"      Sensitivity: {roc.sensitivity_at_optimal:.3f}")
        print(f"      Specificity: {roc.specificity_at_optimal:.3f}")
        print(f"      Missingness: {miss.missingness_rate:.1%}")
        print()

    if result.scores_excluded:
        print("Scores EXCLUDED (high missingness):")
        for name in result.scores_excluded:
            miss = result.missingness[name]
            print(f"  - {name}: {miss.missingness_rate:.1%} missing")
        print()

    # Show output files
    print("-" * 70)
    print("Generated Output Files")
    print("-" * 70)

    results_dir = Path(output_dir) / "results"
    print(f"\nOutput directory: {results_dir}")
    print()

    plots_dir = results_dir / "plots"
    if plots_dir.exists():
        print("Plots:")
        for plot_file in sorted(plots_dir.glob("*.png")):
            print(f"  - {plot_file.name}")
        print()

    print("Data files:")
    for data_file in sorted(results_dir.glob("*.json")) + sorted(results_dir.glob("*.tsv")):
        print(f"  - {data_file.name}")

    print()
    print("=" * 70)
    print("PSROC Example Complete!")
    print("=" * 70)
    print()
    print(f"View the dashboard plot: {results_dir / 'plots' / 'psroc_example_dashboard.png'}")

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run PSROC example with synthetic test data"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="/tmp/psroc_example",
        help="Output directory for results (default: /tmp/psroc_example)",
    )
    args = parser.parse_args()

    sys.exit(main(output_dir=args.output_dir))
