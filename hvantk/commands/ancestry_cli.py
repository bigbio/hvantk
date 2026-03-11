"""
Ancestry Inference CLI - Genetic Ancestry Prediction Command

This module provides the CLI entry point for the ancestry inference pipeline,
which predicts genetic ancestry for samples using PCA and Random Forest
classification against a labeled reference panel.

Example:
    hvantk ancestry-inference \\
        -q /data/cohort.mt \\
        -r /data/1kg_reference.mt \\
        --ancestry-col super_pop \\
        -o /results/ancestry.ht \\
        --generate-report
"""

import logging
from pathlib import Path
from urllib.parse import urlparse

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


@click.command(
    name="ancestry-inference",
    context_settings=CONTEXT_SETTINGS,
)
# Required inputs
@click.option(
    "-q",
    "--query-mt",
    type=str,
    required=True,
    help="Path to query MatrixTable (cohort with unknown ancestry)",
)
@click.option(
    "-r",
    "--reference-mt",
    type=str,
    required=True,
    help="Path to reference MatrixTable with known ancestry labels",
)
@click.option(
    "--ancestry-col",
    type=str,
    default="ancestry",
    help="Column name containing ancestry labels in reference MT [default: ancestry]",
)
# Output options
@click.option(
    "-o",
    "--output-ht",
    type=click.Path(),
    required=True,
    help="Output path for predictions Hail Table",
)
@click.option(
    "--output-dir",
    type=click.Path(),
    default=None,
    help="Directory for additional outputs (TSV, model, report). Defaults to <output-ht parent>/ancestry_results.",
)
# Variant filtering
@click.option(
    "--min-af",
    type=float,
    default=0.01,
    help="Minimum allele frequency for variant inclusion [default: 0.01]",
)
@click.option(
    "--max-af",
    type=float,
    default=0.99,
    help="Maximum allele frequency for variant inclusion [default: 0.99]",
)
@click.option(
    "--min-call-rate",
    type=float,
    default=0.98,
    help="Minimum variant call rate [default: 0.98]",
)
@click.option(
    "--apply-hwe-filter",
    is_flag=True,
    default=False,
    help="Apply Hardy-Weinberg equilibrium filter",
)
@click.option(
    "--hwe-p",
    type=float,
    default=1e-6,
    help="HWE p-value threshold (variants below are removed) [default: 1e-6]",
)
# LD pruning
@click.option(
    "--ld-r2",
    type=float,
    default=0.2,
    help="LD pruning r-squared threshold [default: 0.2]",
)
@click.option(
    "--ld-window",
    type=int,
    default=500000,
    help="LD pruning window size in base pairs [default: 500000]",
)
@click.option(
    "--skip-ld-pruning",
    is_flag=True,
    default=False,
    help="Skip LD pruning step (faster but may affect accuracy)",
)
# PCA options
@click.option(
    "--n-pcs",
    type=int,
    default=20,
    help="Number of principal components to compute [default: 20]",
)
@click.option(
    "--n-pcs-classify",
    type=int,
    default=10,
    help="Number of PCs to use for classification [default: 10]",
)
# Classification options
@click.option(
    "--n-estimators",
    type=int,
    default=100,
    help="Number of trees in Random Forest [default: 100]",
)
@click.option(
    "--min-prob",
    type=float,
    default=0.75,
    help="Minimum probability for ancestry assignment [default: 0.75]",
)
@click.option(
    "--seed",
    type=int,
    default=42,
    help="Random seed for reproducibility [default: 42]",
)
# Validation options
@click.option(
    "--skip-validation",
    is_flag=True,
    default=False,
    help="Skip cross-validation of model",
)
@click.option(
    "--n-cv-folds",
    type=int,
    default=5,
    help="Number of cross-validation folds [default: 5]",
)
# Merge settings
@click.option(
    "--min-shared-variants",
    type=int,
    default=None,
    help="Minimum shared variants required [default: 10000]",
)
# Output formats
@click.option(
    "--generate-report",
    is_flag=True,
    default=False,
    help="Generate HTML report with visualizations (requires visualization extras)",
)
@click.option(
    "--export-tsv",
    is_flag=True,
    default=False,
    help="Export predictions as TSV file",
)
@click.option(
    "--save-model",
    is_flag=True,
    default=False,
    help="Save trained Random Forest model (pickle format)",
)
@click.option(
    "--save-loadings",
    is_flag=True,
    default=False,
    help="Save PCA loadings for projection of new samples",
)
# Checkpointing
@click.option(
    "--checkpoint-path",
    type=click.Path(),
    default=None,
    help="Path for intermediate checkpoints (useful for long-running pipelines)",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite existing outputs and checkpoints",
)
# Logging
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"]),
    default="INFO",
    help="Set logging level [default: INFO]",
)
@click.pass_context
def ancestry_inference_cmd(
    ctx,
    query_mt: str,
    reference_mt: str,
    ancestry_col: str,
    output_ht: str,
    output_dir: str,
    min_af: float,
    max_af: float,
    min_call_rate: float,
    apply_hwe_filter: bool,
    hwe_p: float,
    ld_r2: float,
    ld_window: int,
    skip_ld_pruning: bool,
    n_pcs: int,
    n_pcs_classify: int,
    n_estimators: int,
    min_prob: float,
    seed: int,
    skip_validation: bool,
    n_cv_folds: int,
    min_shared_variants: int,
    generate_report: bool,
    export_tsv: bool,
    save_model: bool,
    save_loadings: bool,
    checkpoint_path: str,
    overwrite: bool,
    log_level: str,
):
    """
    Infer genetic ancestry using PCA and Random Forest classification.

    This command predicts ancestry for samples in a query cohort using a
    labeled reference panel (e.g., 1000 Genomes, HapMap). The workflow
    follows standard population genetics practices:

    \b
    Workflow:
      1. Merge query and reference MatrixTables by shared variants
      2. Filter to high-quality, common, biallelic SNPs
      3. LD prune to remove correlated variants
      4. Compute HWE-normalized PCA
      5. Train Random Forest classifier on reference samples
      6. Predict ancestry for query samples

    \b
    Examples:

      # Basic usage
      hvantk ancestry-inference \\
          -q cohort.mt \\
          -r 1kg_phase3.mt \\
          --ancestry-col super_pop \\
          -o cohort_ancestry.ht

      # With HTML report and TSV export
      hvantk ancestry-inference \\
          -q cohort.mt \\
          -r 1kg_phase3.mt \\
          --ancestry-col super_pop \\
          -o cohort_ancestry.ht \\
          --generate-report \\
          --export-tsv

      # Custom parameters for stricter assignment
      hvantk ancestry-inference \\
          -q cohort.mt \\
          -r 1kg_phase3.mt \\
          --ancestry-col super_pop \\
          -o cohort_ancestry.ht \\
          --min-af 0.05 \\
          --min-prob 0.90 \\
          --n-pcs 30 \\
          --n-pcs-classify 15

      # With checkpointing for large datasets
      hvantk ancestry-inference \\
          -q large_cohort.mt \\
          -r reference.mt \\
          -o ancestry.ht \\
          --checkpoint-path /tmp/ancestry_checkpoints

    \b
    Output Files:
      {output_ht}                          Hail Table with predictions
      {output_dir}/predictions.tsv         TSV export (if --export-tsv)
      {output_dir}/rf_model.pkl            Trained model (if --save-model)
      {output_dir}/pca_loadings.ht         PCA loadings (if --save-loadings)
      {output_dir}/ancestry_report.html    HTML report (if --generate-report)
      {output_dir}/pipeline_stats.json     Pipeline statistics

    \b
    Probability Threshold Guidance:
      0.50  Aggressive assignment; may misclassify admixed individuals
      0.75  Balanced choice (default)
      0.90  Conservative; more "unassigned" but higher confidence
    """

    try:
        # Import Hail and pipeline after logging is configured
        import hail as hl
        from hvantk.core.hail_context import init_hail
        from hvantk.ancestry.pipeline import run_ancestry_inference, PipelineConfig
        from hailtop import fs

        # Initialize Hail (handle case where it's already running externally)
        logger.info("Initializing Hail")
        try:
            init_hail()
        except AssertionError:
            # Hail may already be initialized externally (e.g., in tests)
            logger.debug("Hail appears to be already initialized externally")

        # Validate input paths (cloud-aware)
        if not fs.exists(query_mt):
            raise click.BadParameter(
                f"Query MatrixTable does not exist: {query_mt}", param_hint="--query-mt"
            )
        if not fs.exists(reference_mt):
            raise click.BadParameter(
                f"Reference MatrixTable does not exist: {reference_mt}",
                param_hint="--reference-mt",
            )

        # Validate parameters
        if min_af < 0 or min_af >= 0.5:
            click.echo(
                "Error: --min-af must be between 0 and 0.5 (exclusive)",
                err=True,
            )
            ctx.exit(1)

        if max_af <= 0.5 or max_af > 1:
            click.echo(
                "Error: --max-af must be between 0.5 (exclusive) and 1",
                err=True,
            )
            ctx.exit(1)

        if min_af >= max_af:
            click.echo(
                "Error: --min-af must be less than --max-af",
                err=True,
            )
            ctx.exit(1)

        if min_prob < 0 or min_prob > 1:
            click.echo(
                "Error: --min-prob must be between 0 and 1",
                err=True,
            )
            ctx.exit(1)

        if n_pcs_classify > n_pcs:
            click.echo(
                f"Warning: --n-pcs-classify ({n_pcs_classify}) > --n-pcs ({n_pcs}). "
                f"Will use {n_pcs} PCs for classification.",
                err=True,
            )
            n_pcs_classify = n_pcs

        # Determine output directory
        def _is_cloud_uri(path: str) -> bool:
            parsed = urlparse(path)
            if not parsed.scheme:
                return False
            if parsed.scheme.lower() == "file":
                return False
            if len(path) >= 2 and path[1] == ":" and path[0].isalpha():
                return False
            return True

        output_is_cloud = _is_cloud_uri(output_ht)
        if output_dir is None:
            if output_is_cloud:
                click.echo(
                    "Error: --output-ht points to a cloud URI. "
                    "Provide a local --output-dir for additional outputs.",
                    err=True,
                )
                ctx.exit(1)
            output_dir = str(Path(output_ht).parent / "ancestry_results")

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Load MatrixTables
        click.echo("\n" + "=" * 70)
        click.echo("ANCESTRY INFERENCE PIPELINE")
        click.echo("=" * 70 + "\n")

        logger.info(f"Loading query MatrixTable: {query_mt}")
        click.echo(f"Loading query MatrixTable: {query_mt}")
        query = hl.read_matrix_table(query_mt)

        logger.info(f"Loading reference MatrixTable: {reference_mt}")
        click.echo(f"Loading reference MatrixTable: {reference_mt}")
        reference = hl.read_matrix_table(reference_mt)

        # Build configuration
        config = PipelineConfig(
            min_af=min_af,
            max_af=max_af,
            min_call_rate=min_call_rate,
            apply_hwe_filter=apply_hwe_filter,
            hwe_p_threshold=hwe_p,
            ld_r2=ld_r2,
            ld_window=ld_window,
            skip_ld_pruning=skip_ld_pruning,
            n_pcs=n_pcs,
            n_pcs_classify=n_pcs_classify,
            n_estimators=n_estimators,
            min_prob=min_prob,
            random_seed=seed,
            validate_model=not skip_validation,
            n_cv_folds=n_cv_folds,
            min_shared_variants=min_shared_variants,
            checkpoint_path=checkpoint_path,
            overwrite_checkpoints=overwrite,
        )

        # Run pipeline
        click.echo("\nStarting ancestry inference pipeline...")
        result = run_ancestry_inference(
            query_mt=query,
            reference_mt=reference,
            ancestry_col=ancestry_col,
            config=config,
        )

        # Save main predictions table
        click.echo(f"\nSaving predictions to: {output_ht}")
        result.predictions.write(output_ht, overwrite=overwrite)

        # Generate optional outputs
        if export_tsv or save_model or save_loadings or generate_report:
            click.echo(f"\nGenerating additional outputs in: {output_dir}")

        if export_tsv:
            tsv_path = output_path / "predictions.tsv"
            result.get_predictions_df().to_csv(tsv_path, sep="\t", index=False)
            click.echo(f"  Exported TSV: {tsv_path}")

        if save_model:
            import pickle

            model_path = output_path / "rf_model.pkl"
            with open(model_path, "wb") as f:
                pickle.dump(result.model, f)
            click.echo(f"  Saved model: {model_path}")

        if save_loadings:
            loadings_path = str(output_path / "pca_loadings.ht")
            result.loadings.write(loadings_path, overwrite=overwrite)
            click.echo(f"  Saved loadings: {loadings_path}")

        if generate_report:
            try:
                from hvantk.ancestry.report import generate_ancestry_report

                report_path = output_path / "ancestry_report.html"
                generate_ancestry_report(result, report_path)
                click.echo(f"  Generated report: {report_path}")
            except ImportError as e:
                click.echo(
                    f"  Warning: Could not generate report (missing dependencies): {e}",
                    err=True,
                )

        # Save pipeline statistics
        import json

        stats_path = output_path / "pipeline_stats.json"
        stats_to_save = {
            **result.pipeline_stats,
            "config": result.config.to_dict(),
        }
        with open(stats_path, "w") as f:
            json.dump(stats_to_save, f, indent=2, default=str)
        click.echo(f"  Saved statistics: {stats_path}")

        # Print summary
        predictions_df = result.get_predictions_df()
        from hvantk.ancestry.constants import SOURCE_COL, PREDICTED_ANCESTRY_COL

        query_preds = predictions_df[predictions_df[SOURCE_COL] == "query"]

        click.echo("\n" + "=" * 70)
        click.echo("ANCESTRY INFERENCE COMPLETE")
        click.echo("=" * 70)
        click.echo(f"Query samples:     {len(query_preds)}")
        n_assigned = (query_preds[PREDICTED_ANCESTRY_COL] != "unassigned").sum()
        n_unassigned = (query_preds[PREDICTED_ANCESTRY_COL] == "unassigned").sum()
        click.echo(f"Assigned:          {n_assigned}")
        click.echo(f"Unassigned:        {n_unassigned}")

        if result.get_accuracy() is not None:
            click.echo(f"CV Accuracy:       {result.get_accuracy():.2%}")

        click.echo("\nAncestry distribution:")
        if len(query_preds) == 0:
            click.echo("  No query samples to display")
        else:
            for pop, count in (
                query_preds[PREDICTED_ANCESTRY_COL].value_counts().items()
            ):
                pct = 100 * count / len(query_preds)
                click.echo(f"  {pop}: {count} ({pct:.1f}%)")

        click.echo("\nOutput files:")
        click.echo(f"  Predictions: {output_ht}")
        if export_tsv:
            click.echo(f"  TSV: {output_path / 'predictions.tsv'}")
        if save_model:
            click.echo(f"  Model: {output_path / 'rf_model.pkl'}")
        if save_loadings:
            click.echo(f"  Loadings: {output_path / 'pca_loadings.ht'}")
        if generate_report:
            click.echo(f"  Report: {output_path / 'ancestry_report.html'}")
        click.echo(f"  Stats: {stats_path}")
        click.echo("=" * 70 + "\n")

    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
    except Exception as e:
        logger.exception(f"Pipeline failed: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
