#!/usr/bin/env python
"""Basic ancestry inference example.

This script demonstrates the basic usage of the ancestry inference pipeline:
1. Loading query and reference MatrixTables
2. Running the inference pipeline with default parameters
3. Examining the results
4. Generating visualizations and reports

Usage:
    python basic_inference.py --query-mt /path/to/query.mt \
                              --reference-mt /path/to/reference.mt \
                              --output-dir /path/to/output

For testing without real data, the script can generate synthetic data.
"""

import argparse
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def create_synthetic_data(output_dir: Path):
    """Create synthetic MatrixTables for testing.

    Uses Hail's Balding-Nichols population genetics model to generate
    synthetic genotype data with known population structure.

    Parameters
    ----------
    output_dir : Path
        Directory to save synthetic MatrixTables.

    Returns
    -------
    tuple
        (query_mt_path, reference_mt_path)
    """
    import hail as hl

    logger.info("Generating synthetic data for testing...")

    query_path = str(output_dir / "synthetic_query.mt")
    reference_path = str(output_dir / "synthetic_reference.mt")

    # Generate query cohort (single population, unknown ancestry)
    logger.info("Creating synthetic query MatrixTable (100 samples)")
    query_mt = hl.balding_nichols_model(
        n_populations=1,
        n_samples=100,
        n_variants=5000,
        n_partitions=4,
    )
    query_mt = query_mt.annotate_cols(
        s=hl.str("query_") + hl.str(query_mt.sample_idx)
    )
    query_mt = query_mt.key_cols_by("s")
    query_mt = query_mt.key_rows_by("locus", "alleles")
    query_mt.write(query_path, overwrite=True)

    # Generate reference panel (3 populations with known labels)
    logger.info("Creating synthetic reference MatrixTable (150 samples, 3 populations)")
    ref_mt = hl.balding_nichols_model(
        n_populations=3,
        n_samples=150,  # 50 per population
        n_variants=5000,
        pop_dist=[1 / 3, 1 / 3, 1 / 3],
        fst=[0.1, 0.1, 0.1],  # Population differentiation
        n_partitions=4,
    )

    # Add ancestry labels based on population assignment
    pop_labels = hl.literal(["EUR", "AFR", "EAS"])
    ref_mt = ref_mt.annotate_cols(
        s=hl.str("ref_") + hl.str(ref_mt.sample_idx),
        ancestry=pop_labels[ref_mt.pop],
    )
    ref_mt = ref_mt.key_cols_by("s")
    ref_mt = ref_mt.key_rows_by("locus", "alleles")
    ref_mt.write(reference_path, overwrite=True)

    logger.info(f"Synthetic query MT: {query_path}")
    logger.info(f"Synthetic reference MT: {reference_path}")

    return query_path, reference_path


def run_ancestry_inference_example(
    query_mt_path: str,
    reference_mt_path: str,
    output_dir: Path,
    ancestry_col: str = "ancestry",
):
    """Run the ancestry inference pipeline.

    Parameters
    ----------
    query_mt_path : str
        Path to query MatrixTable.
    reference_mt_path : str
        Path to reference MatrixTable with ancestry labels.
    output_dir : Path
        Output directory for results.
    ancestry_col : str
        Column name containing ancestry labels in reference MT.
    """
    import hail as hl
    from hvantk.ancestry import run_ancestry_inference, PipelineConfig

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load MatrixTables
    logger.info(f"Loading query MatrixTable: {query_mt_path}")
    query_mt = hl.read_matrix_table(query_mt_path)
    logger.info(f"  Samples: {query_mt.count_cols()}")
    logger.info(f"  Variants: {query_mt.count_rows()}")

    logger.info(f"Loading reference MatrixTable: {reference_mt_path}")
    reference_mt = hl.read_matrix_table(reference_mt_path)
    logger.info(f"  Samples: {reference_mt.count_cols()}")
    logger.info(f"  Variants: {reference_mt.count_rows()}")

    # Configure pipeline (using defaults with some adjustments for synthetic data)
    config = PipelineConfig(
        min_af=0.05,  # Slightly higher for synthetic data
        min_call_rate=0.95,
        skip_ld_pruning=True,  # Skip for small synthetic dataset
        n_pcs=10,
        n_pcs_classify=5,
        min_prob=0.75,
        validate_model=True,
        min_shared_variants=100,  # Lower threshold for synthetic data
    )

    # Run ancestry inference
    logger.info("Running ancestry inference pipeline...")
    result = run_ancestry_inference(
        query_mt=query_mt,
        reference_mt=reference_mt,
        ancestry_col=ancestry_col,
        config=config,
    )

    # =====================
    # Examine Results
    # =====================

    # Get predictions DataFrame
    predictions_df = result.get_predictions_df()
    query_predictions = result.get_query_predictions()

    logger.info("\n" + "=" * 60)
    logger.info("RESULTS SUMMARY")
    logger.info("=" * 60)

    # Print prediction summary
    logger.info(f"\nTotal samples: {len(predictions_df)}")
    logger.info(f"Query samples: {len(query_predictions)}")

    # Ancestry distribution
    logger.info("\nAncestry distribution (query samples):")
    for ancestry, count in query_predictions["predicted_ancestry"].value_counts().items():
        pct = 100 * count / len(query_predictions)
        logger.info(f"  {ancestry}: {count} ({pct:.1f}%)")

    # Cross-validation accuracy
    accuracy = result.get_accuracy()
    if accuracy:
        logger.info(f"\nCross-validation accuracy: {accuracy:.2%}")

    # Variance explained by PCs
    var_explained = result.variance_explained()
    logger.info(f"\nVariance explained by PC1: {var_explained[0]:.2%}")
    logger.info(f"Variance explained by PC1-5: {sum(var_explained[:5]):.2%}")

    # =====================
    # Save Results
    # =====================

    logger.info("\n" + "=" * 60)
    logger.info("SAVING OUTPUTS")
    logger.info("=" * 60)

    # Save predictions Table
    predictions_path = str(output_dir / "predictions.ht")
    result.predictions.write(predictions_path, overwrite=True)
    logger.info(f"Saved predictions: {predictions_path}")

    # Export to TSV
    tsv_path = output_dir / "predictions.tsv"
    predictions_df.to_csv(tsv_path, sep="\t", index=False)
    logger.info(f"Exported TSV: {tsv_path}")

    # Generate HTML report
    report_path = output_dir / "ancestry_report.html"
    result.generate_report(report_path)
    logger.info(f"Generated report: {report_path}")

    # =====================
    # Create Visualizations
    # =====================

    logger.info("\n" + "=" * 60)
    logger.info("CREATING VISUALIZATIONS")
    logger.info("=" * 60)

    # PCA scatter plot
    fig = result.plot_pca(pc_x=1, pc_y=2)
    pca_path = output_dir / "pca_plot.png"
    fig.savefig(pca_path, dpi=150, bbox_inches="tight")
    logger.info(f"Saved PCA plot: {pca_path}")

    # Two-panel PCA plot
    fig = result.plot_pca_panel()
    panel_path = output_dir / "pca_panel.png"
    fig.savefig(panel_path, dpi=150, bbox_inches="tight")
    logger.info(f"Saved PCA panel: {panel_path}")

    # Ancestry proportions
    fig = result.plot_ancestry_proportions()
    props_path = output_dir / "ancestry_proportions.png"
    fig.savefig(props_path, dpi=150, bbox_inches="tight")
    logger.info(f"Saved ancestry proportions: {props_path}")

    # Probability distribution
    fig = result.plot_probability_distribution()
    prob_path = output_dir / "probability_distribution.png"
    fig.savefig(prob_path, dpi=150, bbox_inches="tight")
    logger.info(f"Saved probability distribution: {prob_path}")

    # Confusion matrix (if validation was performed)
    fig = result.plot_confusion_matrix()
    if fig:
        cm_path = output_dir / "confusion_matrix.png"
        fig.savefig(cm_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved confusion matrix: {cm_path}")

    logger.info("\n" + "=" * 60)
    logger.info("EXAMPLE COMPLETE")
    logger.info("=" * 60)
    logger.info(f"All outputs saved to: {output_dir}")

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Basic ancestry inference example",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--query-mt",
        type=str,
        default=None,
        help="Path to query MatrixTable. If not provided, generates synthetic data.",
    )
    parser.add_argument(
        "--reference-mt",
        type=str,
        default=None,
        help="Path to reference MatrixTable. If not provided, generates synthetic data.",
    )
    parser.add_argument(
        "--ancestry-col",
        type=str,
        default="ancestry",
        help="Column name containing ancestry labels in reference MT",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./ancestry_example_output",
        help="Output directory for results",
    )
    parser.add_argument(
        "--use-synthetic",
        action="store_true",
        help="Force use of synthetic data even if paths provided",
    )

    args = parser.parse_args()
    output_dir = Path(args.output_dir)

    # Initialize Hail
    import hail as hl
    from hvantk.core.hail_context import init_hail

    logger.info("Initializing Hail...")
    try:
        init_hail()
    except AssertionError:
        logger.debug("Hail already initialized")

    # Determine data source
    if args.use_synthetic or (args.query_mt is None and args.reference_mt is None):
        logger.info("Using synthetic data for demonstration")
        output_dir.mkdir(parents=True, exist_ok=True)
        query_mt_path, reference_mt_path = create_synthetic_data(output_dir)
    else:
        if args.query_mt is None or args.reference_mt is None:
            parser.error("Both --query-mt and --reference-mt required unless using synthetic data")
        query_mt_path = args.query_mt
        reference_mt_path = args.reference_mt

    # Run example
    run_ancestry_inference_example(
        query_mt_path=query_mt_path,
        reference_mt_path=reference_mt_path,
        output_dir=output_dir,
        ancestry_col=args.ancestry_col,
    )


if __name__ == "__main__":
    main()
