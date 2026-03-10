"""
PSROC CLI - Prediction Score ROC Analysis Command

This module provides the CLI entry point for the PSROC pipeline, which evaluates
variant pathogenicity prediction scores against ClinVar truth labels using
ROC curve analysis.

Example:
    hvantk psroc \\
        --genes BRCA1,BRCA2 \\
        --clinvar-ht /data/clinvar.ht \\
        --dbnsfp-ht /data/dbnsfp.ht \\
        --scores "CADD_phred,REVEL_score,MetaLR_score" \\
        --output-dir /results/psroc
"""

import logging
import click

logger = logging.getLogger(__name__)


def _parse_comma_separated(value: str) -> list:
    """Parse comma-separated string into list of stripped values."""
    if not value:
        return []
    return [v.strip() for v in value.split(",") if v.strip()]


def _display_single_result(result, pipeline) -> None:
    """Display results for a single PSROC run."""
    click.echo("\n" + "=" * 70)
    if pipeline.state.errors:
        click.echo("Pipeline completed with errors:")
        for error in pipeline.state.errors:
            click.echo(f"  - {error}")
        click.echo("=" * 70 + "\n")
    else:
        click.echo("Pipeline completed successfully!")
        click.echo("=" * 70)

        click.echo("\nResults Summary:")
        click.echo(f"  Variants analyzed: {result.n_total}")
        click.echo(f"    Pathogenic: {result.n_pathogenic}")
        click.echo(f"    Benign: {result.n_benign}")
        click.echo(f"    Excluded: {result.n_excluded}")

        click.echo(f"\n  Scores evaluated: {len(result.metrics)}")
        if result.metrics:
            sorted_metrics = sorted(
                result.metrics.items(), key=lambda x: x[1].auc, reverse=True
            )
            for name, roc in sorted_metrics:
                click.echo(
                    f"    {name}: AUC={roc.auc:.3f}, "
                    f"threshold={roc.optimal_threshold:.3f}"
                )

        if result.scores_excluded:
            click.echo(
                f"\n  Scores excluded (high missingness): "
                f"{len(result.scores_excluded)}"
            )
            for name in result.scores_excluded:
                miss = result.missingness[name]
                click.echo(f"    {name}: {miss.missingness_rate:.1%} missing")

        click.echo(f"\nOutput directory: {result.output_dir}")
        click.echo("")


def _display_collection_results(results, output_dir) -> None:
    """Display results for a multi-group PSROC run."""
    click.echo("\n" + "=" * 70)
    click.echo("Gene Set Collection - PSROC Results")
    click.echo("=" * 70)

    for group_name, result in sorted(results.items()):
        click.echo(f"\n  [{group_name}]")
        click.echo(f"    Variants: {result.n_total} "
                    f"(P={result.n_pathogenic}, B={result.n_benign})")
        if result.metrics:
            sorted_metrics = sorted(
                result.metrics.items(), key=lambda x: x[1].auc, reverse=True
            )
            for name, roc in sorted_metrics:
                click.echo(
                    f"    {name}: AUC={roc.auc:.3f}"
                )

    click.echo(f"\n  Groups completed: {len(results)}")
    click.echo(f"  Output directory: {output_dir}")
    click.echo("=" * 70 + "\n")


@click.command(name="psroc")
@click.option(
    "--genes",
    type=str,
    default=None,
    help="Comma-separated gene symbols (e.g., BRCA1,BRCA2,TP53)",
)
@click.option(
    "--genes-file",
    type=click.Path(exists=True),
    default=None,
    help="Path to file containing gene symbols (one per line)",
)
@click.option(
    "--variants",
    type=click.Path(exists=True),
    default=None,
    help="Path to variant list file (chr:pos:ref:alt format, one per line)",
)
@click.option(
    "--gene-sets",
    type=click.Path(exists=True),
    default=None,
    help="Path to a gene set collection file (JSON or GMT). "
    "Runs PSROC independently for each named gene set.",
)
@click.option(
    "--hgnc",
    type=click.Path(exists=True),
    default=None,
    help="Path to HGNC data file (TSV or .ht) for gene alias resolution. "
    "When provided, the gene set is expanded to include known aliases "
    "and previous symbols.",
)
@click.option(
    "--clinvar-ht",
    type=click.Path(exists=True),
    required=True,
    help="Path to pre-built ClinVar Hail Table",
)
@click.option(
    "--dbnsfp-ht",
    type=click.Path(exists=True),
    required=True,
    help="Path to pre-built dbNSFP Hail Table",
)
@click.option(
    "--scores",
    type=str,
    required=True,
    help="Comma-separated dbNSFP score field names to evaluate "
    "(e.g., CADD_phred,REVEL_score,MetaLR_score)",
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    required=True,
    help="Output directory for results",
)
@click.option(
    "--reference-genome",
    type=click.Choice(["GRCh37", "GRCh38"]),
    default="GRCh38",
    help="Reference genome build [default: GRCh38]",
)
@click.option(
    "--min-stars",
    type=int,
    default=1,
    help="Minimum ClinVar review status stars (0-4) [default: 1]",
)
@click.option(
    "--max-missingness",
    type=float,
    default=0.3,
    help="Maximum allowed missingness rate per score (0.0-1.0). "
    "Scores exceeding this threshold are excluded [default: 0.3]",
)
@click.option(
    "--threshold-method",
    type=click.Choice(["youden", "closest_to_corner", "f1"]),
    default="youden",
    help="Method for finding optimal classification threshold [default: youden]",
)
@click.option(
    "--output-prefix",
    type=str,
    default="psroc",
    help="Prefix for output file names [default: psroc]",
)
@click.option(
    "--export-tsv",
    is_flag=True,
    default=False,
    help="Export annotated variants as TSV file",
)
@click.option(
    "--no-plots",
    is_flag=True,
    default=False,
    help="Skip generating visualization plots",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite existing output files",
)
@click.option(
    "--min-variants",
    type=int,
    default=10,
    help="Minimum labeled (P+B) variants required for ROC analysis. "
    "Groups below this threshold are skipped [default: 10]",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Show execution plan without running the pipeline",
)
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"]),
    default="INFO",
    help="Set logging level [default: INFO]",
)
@click.pass_context
def psroc_cmd(
    ctx,
    genes,
    genes_file,
    variants,
    gene_sets,
    hgnc,
    clinvar_ht,
    dbnsfp_ht,
    scores,
    output_dir,
    reference_genome,
    min_stars,
    max_missingness,
    threshold_method,
    output_prefix,
    export_tsv,
    no_plots,
    overwrite,
    min_variants,
    dry_run,
    log_level,
):
    """
    PSROC: Prediction Score ROC Analysis

    Evaluate variant pathogenicity prediction scores against ClinVar truth labels
    using ROC curve analysis. This command computes AUC, optimal thresholds, and
    generates visualization plots for each score.

    \b
    Workflow:
      1. Load ClinVar and dbNSFP Hail Tables
      2. Filter variants by genes or variant list
      3. Assign binary labels (Pathogenic/Benign)
      4. Annotate with dbNSFP prediction scores
      5. Compute missingness statistics per score
      6. Compute ROC metrics for qualifying scores
      7. Generate plots, metrics JSON, and reports

    \b
    Input Sources (provide exactly one):
      --genes         Comma-separated gene symbols
      --genes-file    File with gene symbols (one per line)
      --variants      File with variants (chr:pos:ref:alt format)
      --gene-sets     Gene set collection file (JSON/GMT) for multi-group analysis

    \b
    Examples:

      # Evaluate scores for specific genes
      hvantk psroc \\
          --genes BRCA1,BRCA2 \\
          --clinvar-ht /data/clinvar_grch38.ht \\
          --dbnsfp-ht /data/dbnsfp_grch38.ht \\
          --scores "CADD_phred,REVEL_score,MetaLR_score" \\
          --output-dir /results/brca_psroc

      # Use genes from file with stricter missingness threshold
      hvantk psroc \\
          --genes-file /data/cardiac_genes.txt \\
          --clinvar-ht /data/clinvar_grch38.ht \\
          --dbnsfp-ht /data/dbnsfp_grch38.ht \\
          --scores "CADD_phred,REVEL_score" \\
          --output-dir /results/cardiac \\
          --max-missingness 0.1 \\
          --min-stars 2

      # Evaluate specific variants
      hvantk psroc \\
          --variants /data/my_variants.txt \\
          --clinvar-ht /data/clinvar_grch38.ht \\
          --dbnsfp-ht /data/dbnsfp_grch38.ht \\
          --scores "REVEL_score,ClinPred_score" \\
          --output-dir /results/variants \\
          --export-tsv

      # Dry run to preview execution plan
      hvantk psroc \\
          --genes BRCA1 \\
          --clinvar-ht /data/clinvar.ht \\
          --dbnsfp-ht /data/dbnsfp.ht \\
          --scores "CADD_phred" \\
          --output-dir /results \\
          --dry-run

    \b
    Output Files:
      {prefix}_metrics.json       ROC metrics (AUC, thresholds, sensitivity/specificity)
      {prefix}_missingness.json   Per-score missingness statistics
      {prefix}_annotated.ht       Hail Table with labels and scores
      {prefix}_annotated.tsv      TSV export (if --export-tsv)
      plots/{prefix}_*.png        Visualization plots (unless --no-plots)
    """
    try:
        from hvantk.psroc.pipeline import PSROCConfig, PSROCPipeline

        # Parse comma-separated inputs
        genes_list = _parse_comma_separated(genes) if genes else None
        scores_list = _parse_comma_separated(scores)

        if not scores_list:
            click.echo("Error: --scores must contain at least one score name", err=True)
            ctx.exit(1)

        # Validate input sources
        sources = [genes_list, genes_file, variants, gene_sets]
        provided = [s for s in sources if s]
        if len(provided) == 0:
            click.echo(
                "Error: Must provide exactly one of: "
                "--genes, --genes-file, --variants, or --gene-sets",
                err=True,
            )
            ctx.exit(1)
        if len(provided) > 1:
            click.echo(
                "Error: Cannot provide multiple variant sources. "
                "Use only one of: --genes, --genes-file, --variants, or --gene-sets",
                err=True,
            )
            ctx.exit(1)

        # Load gene set collection if provided
        gene_set_collection = None
        if gene_sets:
            from hvantk.utils.gene_sets import load_gene_sets

            collection = load_gene_sets(gene_sets)
            gene_set_collection = {
                gs.name: gs.genes for gs in collection
            }
            click.echo(
                f"Loaded {len(gene_set_collection)} gene sets from {gene_sets}"
            )

        # Create configuration
        config = PSROCConfig(
            genes=genes_list,
            genes_file=genes_file,
            variants_path=variants,
            gene_set_collection=gene_set_collection,
            clinvar_ht=clinvar_ht,
            dbnsfp_ht=dbnsfp_ht,
            scores=scores_list,
            output_dir=output_dir,
            output_prefix=output_prefix,
            reference_genome=reference_genome,
            min_stars=min_stars,
            max_missingness=max_missingness,
            threshold_method=threshold_method,
            hgnc_path=hgnc,
            export_tsv=export_tsv,
            overwrite=overwrite,
            generate_plots=not no_plots,
            min_variants=min_variants,
        )

        # Validate configuration
        errors = config.validate()
        if errors:
            click.echo("Configuration validation failed:", err=True)
            for error in errors:
                click.echo(f"  - {error}", err=True)
            ctx.exit(1)

        # Create pipeline
        pipeline = PSROCPipeline(config)

        # Show plan if dry-run
        if dry_run:
            pipeline.show_plan()
            return

        # Display starting message
        click.echo("\n" + "=" * 70)
        click.echo("Starting PSROC Pipeline Execution")
        click.echo("=" * 70 + "\n")

        # Run pipeline
        if config.gene_set_collection:
            results = pipeline.run_collection()
            _display_collection_results(results, output_dir)
        else:
            result = pipeline.run()
            _display_single_result(result, pipeline)

    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
    except Exception as e:
        logger.exception(f"Pipeline failed: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
