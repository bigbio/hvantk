"""
HGC QC CLI Commands Module

This module provides command-line interface commands for quality control operations
in the HGC (Hail-based Genotype Combiner) toolkit.

Commands:
- hvantk hgc compute-qc: Compute QC metrics for samples and variants
- hvantk hgc filter-qc: Filter MatrixTable based on QC metrics
- hvantk hgc qc-summary: Generate summary statistics from QC metrics files
- hvantk hgc qc-report: Generate comprehensive HTML QC report with embedded plots
"""

import logging
import click
import os

logger = logging.getLogger(__name__)

# Import HGC QC functionality
from hvantk.algorithms.hgc import (
    check_path_exists_and_readable,
    compute_full_qc,
    compute_sample_qc,
    compute_variant_qc,
    filter_samples_by_qc,
    filter_variants_by_qc,
    save_qc_metrics,
)

# Import utility functions
from .utils import validate_input_files, validate_output_path


def register_qc_commands(group):
    """
    Register QC commands to the HGC command group.

    Args:
        group: Click command group to add commands to
    """
    group.add_command(compute_qc)
    group.add_command(filter_qc)
    group.add_command(qc_summary)
    group.add_command(qc_report)


@click.command(name="compute-qc")
@click.option("--input", "-i", required=True, help="Input MatrixTable path")
@click.option(
    "--output-dir", "-o", required=True, help="Output directory for QC metrics"
)
@click.option(
    "--sample-qc/--no-sample-qc", default=True, help="Compute sample-level QC metrics"
)
@click.option(
    "--variant-qc/--no-variant-qc",
    default=True,
    help="Compute variant-level QC metrics",
)
@click.option(
    "--call-field", default="GT", help="Name of the call field to use (default: GT)"
)
@click.option("--prefix", default="qc_metrics", help="Prefix for output files")
@click.option(
    "--save-mt/--no-save-mt", default=True, help="Save MatrixTable with QC annotations"
)
@click.option(
    "--dry-run", is_flag=True, help="Show what would be done without executing"
)
@click.pass_context
def compute_qc(
    ctx, input, output_dir, sample_qc, variant_qc, call_field, prefix, save_mt, dry_run
):
    """
    Compute comprehensive quality control metrics for samples and variants.

    This command computes QC metrics including call rates, allele frequencies,
    Hardy-Weinberg equilibrium tests, and genotype quality statistics for
    both samples and variants in a MatrixTable.

    Examples:
        hvantk hgc compute-qc -i cohort.mt -o qc_results/
        hvantk hgc compute-qc -i cohort.mt -o qc_results/ --no-sample-qc
        hvantk hgc compute-qc -i cohort.mt -o qc_results/ --call-field LGT
    """
    try:
        logger.info("Starting QC metrics computation")

        # Validate input
        is_valid, errors = validate_input_files([input], "mt")
        if not is_valid:
            click.echo("❌ Input file validation failed:", err=True)
            for error in errors:
                click.echo(f"   • {error}", err=True)
            ctx.exit(1)

        # Validate output directory
        if not validate_output_path(output_dir, create_dirs=True):
            click.echo("❌ Invalid output directory", err=True)
            ctx.exit(1)

        if not sample_qc and not variant_qc:
            click.echo(
                "❌ At least one of --sample-qc or --variant-qc must be enabled",
                err=True,
            )
            ctx.exit(1)

        if dry_run:
            click.echo("🔍 Dry run mode - would execute QC computation with:")
            click.echo(f"   • Input: {input}")
            click.echo(f"   • Output directory: {output_dir}")
            click.echo(f"   • Sample QC: {sample_qc}")
            click.echo(f"   • Variant QC: {variant_qc}")
            click.echo(f"   • Call field: {call_field}")
            click.echo(f"   • Prefix: {prefix}")
            click.echo(f"   • Save MatrixTable: {save_mt}")
            return

        # Import hail_context (init_hail) before hail: it applies the NumPy
        # np.bool compatibility shim before importing Hail.
        from hvantk.core.utils.hail_context import init_hail
        import hail as hl

        init_hail(quiet=True)

        click.echo("🔄 Loading MatrixTable...")
        mt = hl.read_matrix_table(input)

        # Compute QC metrics based on options
        if sample_qc and variant_qc:
            click.echo("🔄 Computing comprehensive QC metrics...")
            qc_results = compute_full_qc(mt, call_field=call_field)
        elif sample_qc:
            click.echo("🔄 Computing sample QC metrics...")
            mt_qc = compute_sample_qc(mt, call_field=call_field)
            from hvantk.algorithms.hgc.qc import QCMetrics

            sample_qc_table = mt_qc.cols().select("sample_qc")
            qc_results = QCMetrics(mt_qc, sample_qc_table, None)
        else:  # variant_qc only
            click.echo("🔄 Computing variant QC metrics...")
            mt_qc = compute_variant_qc(mt, call_field=call_field)
            from hvantk.algorithms.hgc.qc import QCMetrics

            variant_qc_table = mt_qc.rows().select("variant_qc")
            qc_results = QCMetrics(mt_qc, None, variant_qc_table)

        # Save QC metrics
        click.echo("💾 Saving QC metrics...")
        saved_files = save_qc_metrics(qc_results, output_dir, prefix)

        # Remove MatrixTable from saved files if not requested
        if not save_mt and "matrix_table" in saved_files:
            import shutil

            if os.path.exists(saved_files["matrix_table"]):
                shutil.rmtree(saved_files["matrix_table"])
            del saved_files["matrix_table"]

        click.echo("✅ Successfully computed and saved QC metrics:")
        for file_type, file_path in saved_files.items():
            click.echo(f"   • {file_type}: {file_path}")

    except Exception as e:
        logger.exception(f"QC computation failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)


@click.command(name="filter-qc")
@click.option(
    "--input", "-i", required=True, help="Input MatrixTable with QC annotations"
)
@click.option("--output", "-o", required=True, help="Output filtered MatrixTable path")
@click.option(
    "--min-sample-call-rate", type=float, default=0.85, help="Minimum sample call rate"
)
@click.option(
    "--min-variant-call-rate",
    type=float,
    default=0.85,
    help="Minimum variant call rate",
)
@click.option("--min-ac", type=int, default=1, help="Minimum allele count for variants")
@click.option("--max-ac", type=int, help="Maximum allele count for variants")
@click.option("--min-af", type=float, help="Minimum allele frequency for variants")
@click.option("--max-af", type=float, help="Maximum allele frequency for variants")
@click.option(
    "--hwe-threshold",
    type=float,
    default=1e-6,
    help="Hardy-Weinberg equilibrium p-value threshold",
)
@click.option("--min-mean-dp", type=float, help="Minimum mean depth for samples")
@click.option("--max-mean-dp", type=float, help="Maximum mean depth for samples")
@click.option(
    "--min-mean-gq", type=float, help="Minimum mean genotype quality for samples"
)
@click.option(
    "--sample-qc-name", default="sample_qc", help="Name of sample QC annotation"
)
@click.option(
    "--variant-qc-name", default="variant_qc", help="Name of variant QC annotation"
)
@click.option(
    "--overwrite/--no-overwrite", default=False, help="Overwrite output if exists"
)
@click.option(
    "--dry-run", is_flag=True, help="Show what would be done without executing"
)
@click.pass_context
def filter_qc(
    ctx,
    input,
    output,
    min_sample_call_rate,
    min_variant_call_rate,
    min_ac,
    max_ac,
    min_af,
    max_af,
    hwe_threshold,
    min_mean_dp,
    max_mean_dp,
    min_mean_gq,
    sample_qc_name,
    variant_qc_name,
    overwrite,
    dry_run,
):
    """
    Filter MatrixTable based on quality control metrics.

    Apply comprehensive QC filters to remove low-quality samples and variants
    based on call rates, allele frequencies, Hardy-Weinberg equilibrium,
    and other quality metrics.

    Examples:
        hvantk hgc filter-qc -i cohort_qc.mt -o cohort_filtered.mt
        hvantk hgc filter-qc -i cohort_qc.mt -o cohort_filtered.mt --min-ac 5 --min-af 0.01
        hvantk hgc filter-qc -i cohort_qc.mt -o cohort_filtered.mt --min-sample-call-rate 0.9
    """
    try:
        logger.info("Starting QC-based filtering")

        # Validate input
        is_valid, errors = validate_input_files([input], "mt")
        if not is_valid:
            click.echo("❌ Input file validation failed:", err=True)
            for error in errors:
                click.echo(f"   • {error}", err=True)
            ctx.exit(1)

        # Validate output path
        if not validate_output_path(output, create_dirs=True):
            click.echo("❌ Invalid output path", err=True)
            ctx.exit(1)

        if dry_run:
            click.echo("🔍 Dry run mode - would execute QC filtering with:")
            click.echo(f"   • Input: {input}")
            click.echo(f"   • Output: {output}")
            click.echo(f"   • Sample filters:")
            click.echo(f"     - Min call rate: {min_sample_call_rate}")
            if min_mean_dp:
                click.echo(f"     - Min mean depth: {min_mean_dp}")
            if max_mean_dp:
                click.echo(f"     - Max mean depth: {max_mean_dp}")
            if min_mean_gq:
                click.echo(f"     - Min mean GQ: {min_mean_gq}")
            click.echo(f"   • Variant filters:")
            click.echo(f"     - Min call rate: {min_variant_call_rate}")
            click.echo(f"     - Min AC: {min_ac}")
            if max_ac:
                click.echo(f"     - Max AC: {max_ac}")
            if min_af:
                click.echo(f"     - Min AF: {min_af}")
            if max_af:
                click.echo(f"     - Max AF: {max_af}")
            click.echo(f"     - HWE threshold: {hwe_threshold}")
            return

        # Import hail_context (init_hail) before hail: it applies the NumPy
        # np.bool compatibility shim before importing Hail.
        from hvantk.core.utils.hail_context import init_hail
        import hail as hl

        init_hail(quiet=True)

        click.echo("🔄 Loading MatrixTable...")
        mt = hl.read_matrix_table(input)

        # Get initial counts
        n_samples_initial = mt.count_cols()
        n_variants_initial = mt.count_rows()

        # Apply sample filters
        click.echo("🔄 Filtering samples based on QC metrics...")
        try:
            mt_filtered = filter_samples_by_qc(
                mt,
                min_call_rate=min_sample_call_rate,
                min_mean_dp=min_mean_dp,
                max_mean_dp=max_mean_dp,
                min_mean_gq=min_mean_gq,
                sample_qc_name=sample_qc_name,
            )
        except ValueError as e:
            if "Sample QC annotation" in str(e):
                click.echo(f"⚠️  Skipping sample filtering: {e}")
                mt_filtered = mt
            else:
                raise

        # Apply variant filters
        click.echo("🔄 Filtering variants based on QC metrics...")
        try:
            mt_filtered = filter_variants_by_qc(
                mt_filtered,
                min_call_rate=min_variant_call_rate,
                min_ac=min_ac,
                max_ac=max_ac,
                min_af=min_af,
                max_af=max_af,
                hwe_threshold=hwe_threshold,
                variant_qc_name=variant_qc_name,
            )
        except ValueError as e:
            if "Variant QC annotation" in str(e):
                click.echo(f"⚠️  Skipping variant filtering: {e}")
            else:
                raise

        # Get final counts
        n_samples_final = mt_filtered.count_cols()
        n_variants_final = mt_filtered.count_rows()

        # Save filtered MatrixTable
        click.echo("💾 Saving filtered MatrixTable...")
        mt_filtered.write(output, overwrite=overwrite)

        # Report filtering results
        samples_removed = n_samples_initial - n_samples_final
        variants_removed = n_variants_initial - n_variants_final

        click.echo("✅ Successfully applied QC filters:")
        click.echo(
            f"   • Samples: {n_samples_initial} → {n_samples_final} "
            f"({samples_removed} removed, {samples_removed/n_samples_initial*100:.1f}%)"
        )
        click.echo(
            f"   • Variants: {n_variants_initial} → {n_variants_final} "
            f"({variants_removed} removed, {variants_removed/n_variants_initial*100:.1f}%)"
        )
        click.echo(f"   • Output: {output}")

    except Exception as e:
        logger.exception(f"QC filtering failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)


@click.command(name="qc-summary")
@click.option(
    "--qc-dir", "-d", required=True, help="Directory containing QC metrics files"
)
@click.option("--sample-file", help="Specific sample QC file (optional)")
@click.option("--variant-file", help="Specific variant QC file (optional)")
@click.option("--output", "-o", help="Output file for summary report")
@click.option(
    "--format",
    type=click.Choice(["csv", "json", "markdown"]),
    default="markdown",
    help="Output format for summary",
)
@click.pass_context
def qc_summary(ctx, qc_dir, sample_file, variant_file, output, format):
    """
    Generate summary statistics from QC metrics files.

    Create comprehensive summary reports of quality control metrics
    including distributions, outliers, and recommended filtering thresholds.

    Examples:
        hvantk hgc qc-summary -d qc_results/
        hvantk hgc qc-summary -d qc_results/ -o summary.md --format markdown
        hvantk hgc qc-summary --sample-file sample_qc.csv --variant-file variant_qc.csv
    """
    try:
        logger.info("Generating QC summary statistics")

        import pandas as pd
        from pathlib import Path

        qc_dir = Path(qc_dir)

        # Find QC files
        if sample_file:
            sample_path = Path(sample_file)
        else:
            sample_files = list(qc_dir.glob("*sample_qc*.csv"))
            if not sample_files:
                click.echo("⚠️  No sample QC files found")
                sample_path = None
            else:
                sample_path = sample_files[0]
                if len(sample_files) > 1:
                    click.echo(
                        f"ℹ️  Multiple sample QC files found, using {sample_path.name}"
                    )

        if variant_file:
            variant_path = Path(variant_file)
        else:
            variant_files = list(qc_dir.glob("*variant_qc*.csv"))
            if not variant_files:
                click.echo("⚠️  No variant QC files found")
                variant_path = None
            else:
                variant_path = variant_files[0]
                if len(variant_files) > 1:
                    click.echo(
                        f"ℹ️  Multiple variant QC files found, using {variant_path.name}"
                    )

        if not sample_path and not variant_path:
            click.echo("❌ No QC files found", err=True)
            ctx.exit(1)

        # Generate summary
        summary_data = {}

        if sample_path and sample_path.exists():
            click.echo(f"📊 Processing sample QC metrics from {sample_path}")
            sample_df = pd.read_csv(sample_path)
            from hvantk.algorithms.hgc.qc import get_qc_summary_stats

            sample_summary = get_qc_summary_stats(sample_df)
            summary_data["sample_qc"] = {
                "file": str(sample_path),
                "n_samples": len(sample_df),
                "summary_stats": (
                    sample_summary.to_dict() if not sample_summary.empty else {}
                ),
                "columns": sample_df.columns.tolist(),
            }

        if variant_path and variant_path.exists():
            click.echo(f"📊 Processing variant QC metrics from {variant_path}")
            variant_df = pd.read_csv(variant_path)
            from hvantk.algorithms.hgc.qc import get_qc_summary_stats

            variant_summary = get_qc_summary_stats(variant_df)
            summary_data["variant_qc"] = {
                "file": str(variant_path),
                "n_variants": len(variant_df),
                "summary_stats": (
                    variant_summary.to_dict() if not variant_summary.empty else {}
                ),
                "columns": variant_df.columns.tolist(),
            }

        # Output summary
        if output:
            output_path = Path(output)
            if format == "csv":
                # Combine summaries into CSV
                all_stats = []
                for qc_type, data in summary_data.items():
                    for metric, stats in data["summary_stats"].items():
                        for stat_name, value in stats.items():
                            all_stats.append(
                                {
                                    "qc_type": qc_type,
                                    "metric": metric,
                                    "statistic": stat_name,
                                    "value": value,
                                }
                            )
                summary_df = pd.DataFrame(all_stats)
                summary_df.to_csv(output_path, index=False)

            elif format == "json":
                import json

                with open(output_path, "w") as f:
                    json.dump(summary_data, f, indent=2)

            elif format == "markdown":
                from hvantk.algorithms.hgc.qc_report import render_qc_summary_markdown

                md_content = render_qc_summary_markdown(summary_data)

                with open(output_path, "w") as f:
                    f.write(md_content)

            click.echo(f"💾 Summary saved to {output_path}")
        else:
            # Print to console
            click.echo("\n📊 QC Summary Statistics:")
            for qc_type, data in summary_data.items():
                click.echo(f"\n{qc_type.replace('_', ' ').title()}:")
                click.echo(f"  File: {data['file']}")
                click.echo(
                    f"  Count: {data.get('n_samples', data.get('n_variants', 0)):,}"
                )
                click.echo(f"  Metrics: {len(data['columns'])}")

        click.echo("✅ QC summary completed successfully")

    except Exception as e:
        logger.exception(f"QC summary generation failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)


@click.command(name="qc-report")
@click.option(
    "--input", "-i", required=True, help="Input MatrixTable with QC annotations"
)
@click.option("--output", "-o", required=True, help="Output HTML report file")
@click.option("--title", default="Quality Control Report", help="Report title")
@click.option(
    "--include-plots",
    multiple=True,
    type=click.Choice(
        [
            "sample_overview",
            "variant_overview",
            "sample_call_rates",
            "variant_call_rates",
            "allele_frequencies",
            "hwe",
            "titv",
        ]
    ),
    help="Specific plots to include (can be used multiple times)",
)
@click.option(
    "--style",
    type=click.Choice(["default", "publication"]),
    default="default",
    help="Plot style for embedded plots",
)
@click.option(
    "--dry-run", is_flag=True, help="Show what would be done without executing"
)
@click.pass_context
def qc_report(ctx, input, output, title, include_plots, style, dry_run):
    """
    Generate comprehensive HTML QC report with embedded plots.

    Create a self-contained HTML report with quality control metrics, plots,
    summary tables, and recommendations for genomic variant data.

    Examples:
        hvantk hgc qc-report -i cohort_qc.mt -o qc_report.html
        hvantk hgc qc-report -i cohort_qc.mt -o report.html --title "Cohort QC Analysis"
        hvantk hgc qc-report -i cohort_qc.mt -o report.html --include-plots sample_overview variant_overview
        hvantk hgc qc-report -i cohort_qc.mt -o report.html --style publication --dry-run
    """
    try:
        # hail_context (init_hail) first: it applies the NumPy np.bool
        # compatibility shim before Hail is imported.
        from hvantk.core.utils.hail_context import init_hail
        import hail as hl
        from pathlib import Path

        # Validate inputs
        if not check_path_exists_and_readable(input):
            click.echo(f"❌ Input file not found or not readable: {input}", err=True)
            ctx.exit(1)

        output_path = Path(output)

        if dry_run:
            click.echo("🔍 Dry run mode - showing what would be done:")
            click.echo(f"   • Input MatrixTable: {input}")
            click.echo(f"   • Output HTML report: {output_path.absolute()}")
            click.echo(f"   • Report title: {title}")
            click.echo(f"   • Plot style: {style}")
            if include_plots:
                click.echo(f"   • Included plots: {', '.join(include_plots)}")
            else:
                click.echo("   • Included plots: all (default)")
            return

        # Initialize Hail
        init_hail(quiet=True)

        # Load MatrixTable
        click.echo(f"📥 Loading MatrixTable from {input}")
        mt = hl.read_matrix_table(input)

        # Check for QC annotations
        if "sample_qc" not in mt.col and "variant_qc" not in mt.row:
            click.echo(
                "❌ No QC annotations found in MatrixTable. Run compute-qc first.",
                err=True,
            )
            ctx.exit(1)

        # Extract QC metrics and create QCMetrics object
        from hvantk.algorithms.hgc.qc import QCMetrics

        sample_qc = mt.cols().select("sample_qc") if "sample_qc" in mt.col else None
        variant_qc = mt.rows().select("variant_qc") if "variant_qc" in mt.row else None

        qc_results = QCMetrics(mt=mt, sample_qc=sample_qc, variant_qc=variant_qc)

        click.echo(
            f"📊 Found QC data: Sample QC: {qc_results.has_sample_qc}, Variant QC: {qc_results.has_variant_qc}"
        )

        # Prepare report parameters
        report_kwargs = {"title": title}

        if include_plots:
            report_kwargs["include_plots"] = list(include_plots)

        # Generate HTML report
        click.echo("📝 Generating comprehensive HTML QC report...")

        report_path = qc_results.generate_html_report(
            output_path=output_path, **report_kwargs
        )

        # Get file size
        file_size = report_path.stat().st_size / 1024  # KB

        click.echo(f"📊 HTML QC report created:")
        click.echo(f"   • File: {report_path.absolute()}")
        click.echo(f"   • Size: {file_size:.1f} KB")
        click.echo(f"   • Title: {title}")

        sample_df = (
            qc_results.get_sample_metrics_df() if qc_results.has_sample_qc else None
        )
        variant_df = (
            qc_results.get_variant_metrics_df() if qc_results.has_variant_qc else None
        )

        if sample_df is not None:
            click.echo(f"   • Samples: {len(sample_df):,}")
        if variant_df is not None:
            click.echo(f"   • Variants: {len(variant_df):,}")

        click.echo(f"\n🌐 Open the report in your web browser:")
        click.echo(f"   file://{report_path.absolute()}")
        click.echo("✅ QC HTML report generated successfully")

    except Exception as e:
        logger.exception(f"QC HTML report generation failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)
