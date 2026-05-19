#!/usr/bin/env python
"""Custom reference panel example.

This script demonstrates how to:
1. Prepare a custom reference panel for ancestry inference
2. Handle different ancestry label formats
3. Configure pipeline parameters for specific use cases
4. Combine results with downstream analysis

Usage:
    python custom_reference.py --reference-mt /path/to/custom_reference.mt \
                               --query-mt /path/to/query.mt \
                               --ancestry-col population \
                               --output-dir /path/to/output
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def prepare_reference_panel(
    mt_path: str,
    ancestry_col: str,
    population_mapping: Optional[Dict[str, str]] = None,
    min_samples_per_pop: int = 20,
    output_path: Optional[str] = None,
):
    """Prepare a reference panel for ancestry inference.

    This function helps standardize a reference panel by:
    - Validating the ancestry column exists
    - Optionally remapping population labels
    - Filtering populations with insufficient samples
    - Computing basic statistics

    Parameters
    ----------
    mt_path : str
        Path to input reference MatrixTable.
    ancestry_col : str
        Column name containing ancestry labels.
    population_mapping : dict, optional
        Mapping from original labels to standardized labels.
        E.g., {"CEU": "EUR", "YRI": "AFR", "CHB": "EAS"}
    min_samples_per_pop : int
        Minimum samples required per population.
    output_path : str, optional
        Path to save prepared MatrixTable. If None, returns without saving.

    Returns
    -------
    hl.MatrixTable
        Prepared reference MatrixTable.
    """
    import hail as hl

    logger.info(f"Loading reference MatrixTable: {mt_path}")
    mt = hl.read_matrix_table(mt_path)

    # Validate ancestry column exists
    if ancestry_col not in mt.col:
        available_cols = list(mt.col)
        raise ValueError(
            f"Ancestry column '{ancestry_col}' not found. "
            f"Available columns: {available_cols}"
        )

    # Get population statistics
    pop_counts = mt.aggregate_cols(hl.agg.counter(mt[ancestry_col]))
    logger.info(f"Original populations: {len(pop_counts)}")
    for pop, count in sorted(pop_counts.items(), key=lambda x: -x[1]):
        logger.info(f"  {pop}: {count} samples")

    # Apply population mapping if provided
    if population_mapping:
        logger.info("Applying population label mapping...")
        mapping_expr = hl.literal(population_mapping)
        original_ancestry = mt[ancestry_col]

        # Map labels, keeping original if not in mapping
        mt = mt.annotate_cols(
            **{
                ancestry_col: hl.coalesce(
                    mapping_expr.get(original_ancestry),
                    original_ancestry,
                )
            }
        )

        # Log new population counts
        new_counts = mt.aggregate_cols(hl.agg.counter(mt[ancestry_col]))
        logger.info(f"Populations after mapping: {len(new_counts)}")
        for pop, count in sorted(new_counts.items(), key=lambda x: -x[1]):
            logger.info(f"  {pop}: {count} samples")

    # Filter populations with insufficient samples
    current_counts = mt.aggregate_cols(hl.agg.counter(mt[ancestry_col]))
    pops_to_keep = [
        pop for pop, count in current_counts.items() if count >= min_samples_per_pop
    ]

    if len(pops_to_keep) < len(current_counts):
        removed = set(current_counts.keys()) - set(pops_to_keep)
        logger.warning(
            f"Removing populations with <{min_samples_per_pop} samples: {removed}"
        )
        mt = mt.filter_cols(hl.literal(pops_to_keep).contains(mt[ancestry_col]))

    # Final statistics
    n_samples = mt.count_cols()
    n_variants = mt.count_rows()
    final_counts = mt.aggregate_cols(hl.agg.counter(mt[ancestry_col]))

    logger.info("\nPrepared reference panel:")
    logger.info(f"  Total samples: {n_samples}")
    logger.info(f"  Total variants: {n_variants}")
    logger.info(f"  Populations: {len(final_counts)}")

    # Save if output path provided
    if output_path:
        logger.info(f"Saving prepared reference: {output_path}")
        mt.write(output_path, overwrite=True)

    return mt


def run_with_custom_reference(
    query_mt_path: str,
    reference_mt_path: str,
    ancestry_col: str,
    output_dir: Path,
    population_mapping: Optional[Dict[str, str]] = None,
    custom_colors: Optional[Dict[str, str]] = None,
):
    """Run ancestry inference with a custom reference panel.

    Parameters
    ----------
    query_mt_path : str
        Path to query MatrixTable.
    reference_mt_path : str
        Path to reference MatrixTable.
    ancestry_col : str
        Column name with ancestry labels.
    output_dir : Path
        Output directory.
    population_mapping : dict, optional
        Population label mapping.
    custom_colors : dict, optional
        Custom colors for populations in plots.
    """
    import hail as hl
    from hvantk.algorithms.ancestry import run_ancestry_inference, PipelineConfig

    output_dir.mkdir(parents=True, exist_ok=True)

    # Prepare reference panel
    logger.info("Preparing reference panel...")
    reference_mt = prepare_reference_panel(
        mt_path=reference_mt_path,
        ancestry_col=ancestry_col,
        population_mapping=population_mapping,
        min_samples_per_pop=20,
    )

    # Load query MatrixTable
    logger.info(f"Loading query MatrixTable: {query_mt_path}")
    query_mt = hl.read_matrix_table(query_mt_path)

    # Configure pipeline for custom reference
    config = PipelineConfig(
        # Variant filtering - adjust based on your data
        min_af=0.01,
        max_af=0.99,
        min_call_rate=0.98,
        # LD pruning
        ld_r2=0.2,
        ld_window=500000,
        # PCA - may need more PCs for diverse references
        n_pcs=20,
        n_pcs_classify=10,
        # Classification
        n_estimators=100,
        min_prob=0.75,
        # Validation
        validate_model=True,
        n_cv_folds=5,
    )

    # Run inference
    logger.info("Running ancestry inference...")
    result = run_ancestry_inference(
        query_mt=query_mt,
        reference_mt=reference_mt,
        ancestry_col=ancestry_col,
        config=config,
    )

    # Print results
    predictions = result.get_query_predictions()
    logger.info("\nResults:")
    logger.info(f"  Query samples: {len(predictions)}")
    logger.info(f"  CV Accuracy: {result.get_accuracy():.2%}")
    logger.info("\n  Ancestry distribution:")
    for pop, count in predictions["predicted_ancestry"].value_counts().items():
        pct = 100 * count / len(predictions)
        logger.info(f"    {pop}: {count} ({pct:.1f}%)")

    # Save outputs
    logger.info("\nSaving outputs...")

    # Predictions
    result.predictions.write(str(output_dir / "predictions.ht"), overwrite=True)

    # TSV export
    result.get_predictions_df().to_csv(
        output_dir / "predictions.tsv", sep="\t", index=False
    )

    # Generate report
    result.generate_report(output_dir / "ancestry_report.html")

    # Custom PCA plot with custom colors
    if custom_colors:
        fig = result.plot_pca(colors=custom_colors)
    else:
        fig = result.plot_pca()
    fig.savefig(output_dir / "pca_plot.png", dpi=150, bbox_inches="tight")

    logger.info(f"\nAll outputs saved to: {output_dir}")

    return result


def create_synthetic_custom_reference(output_dir: Path):
    """Create synthetic data with custom population labels for testing."""
    import hail as hl

    logger.info("Creating synthetic data with custom labels...")

    # Query data
    query_path = str(output_dir / "query.mt")
    query_mt = hl.balding_nichols_model(
        n_populations=1,
        n_samples=50,
        n_variants=3000,
        n_partitions=4,
    )
    query_mt = query_mt.annotate_cols(s=hl.str("sample_") + hl.str(query_mt.sample_idx))
    query_mt = query_mt.key_cols_by("s")
    query_mt = query_mt.key_rows_by("locus", "alleles")
    query_mt.write(query_path, overwrite=True)

    # Reference data with custom population labels
    ref_path = str(output_dir / "reference.mt")
    ref_mt = hl.balding_nichols_model(
        n_populations=4,
        n_samples=120,  # 30 per population
        n_variants=3000,
        pop_dist=[0.25, 0.25, 0.25, 0.25],
        fst=[0.1, 0.1, 0.1, 0.1],
        n_partitions=4,
    )

    # Use custom population labels (like HapMap)
    custom_labels = hl.literal(["CEU", "YRI", "CHB", "GIH"])
    ref_mt = ref_mt.annotate_cols(
        s=hl.str("ref_") + hl.str(ref_mt.sample_idx),
        population=custom_labels[ref_mt.pop],  # Note: using 'population' not 'ancestry'
    )
    ref_mt = ref_mt.key_cols_by("s")
    ref_mt = ref_mt.key_rows_by("locus", "alleles")
    ref_mt.write(ref_path, overwrite=True)

    return query_path, ref_path


def main():
    parser = argparse.ArgumentParser(
        description="Ancestry inference with custom reference panel",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--query-mt",
        type=str,
        help="Path to query MatrixTable",
    )
    parser.add_argument(
        "--reference-mt",
        type=str,
        help="Path to reference MatrixTable",
    )
    parser.add_argument(
        "--ancestry-col",
        type=str,
        default="population",
        help="Column name with ancestry labels",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./custom_reference_output",
        help="Output directory",
    )
    parser.add_argument(
        "--use-synthetic",
        action="store_true",
        help="Use synthetic data for testing",
    )

    args = parser.parse_args()
    output_dir = Path(args.output_dir)

    # Initialize Hail
    import hail as hl
    from hvantk.core.utils.hail_context import init_hail

    logger.info("Initializing Hail...")
    try:
        init_hail()
    except AssertionError:
        logger.debug("Hail already initialized")

    # Example population mapping (HapMap to super-populations)
    population_mapping = {
        "CEU": "EUR",  # Utah residents, European ancestry
        "TSI": "EUR",  # Toscani, Italian
        "GBR": "EUR",  # British
        "FIN": "EUR",  # Finnish
        "YRI": "AFR",  # Yoruba, Nigerian
        "LWK": "AFR",  # Luhya, Kenyan
        "ASW": "AFR",  # African Americans
        "CHB": "EAS",  # Han Chinese, Beijing
        "JPT": "EAS",  # Japanese, Tokyo
        "CHS": "EAS",  # Southern Han Chinese
        "GIH": "SAS",  # Gujarati Indians
        "PJL": "SAS",  # Punjabi
        "BEB": "SAS",  # Bengali
        "MXL": "AMR",  # Mexican ancestry
        "PUR": "AMR",  # Puerto Rican
        "CLM": "AMR",  # Colombian
    }

    # Custom colors for visualization
    custom_colors = {
        "EUR": "#1f77b4",  # Blue
        "AFR": "#ff7f0e",  # Orange
        "EAS": "#2ca02c",  # Green
        "SAS": "#d62728",  # Red
        "AMR": "#9467bd",  # Purple
        "unassigned": "#7f7f7f",  # Gray
    }

    # Determine data source
    if args.use_synthetic or (args.query_mt is None and args.reference_mt is None):
        logger.info("Using synthetic data for demonstration")
        output_dir.mkdir(parents=True, exist_ok=True)
        query_mt_path, reference_mt_path = create_synthetic_custom_reference(output_dir)
        # For synthetic data, use simpler mapping
        population_mapping = {
            "CEU": "EUR",
            "YRI": "AFR",
            "CHB": "EAS",
            "GIH": "SAS",
        }
    else:
        if args.query_mt is None or args.reference_mt is None:
            parser.error("Both --query-mt and --reference-mt required")
        query_mt_path = args.query_mt
        reference_mt_path = args.reference_mt

    # Run with custom reference
    run_with_custom_reference(
        query_mt_path=query_mt_path,
        reference_mt_path=reference_mt_path,
        ancestry_col=args.ancestry_col,
        output_dir=output_dir,
        population_mapping=population_mapping,
        custom_colors=custom_colors,
    )


if __name__ == "__main__":
    main()
