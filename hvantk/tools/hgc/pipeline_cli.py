"""
HGC CLI - Pipeline Command

End-to-end pipeline command for complete gVCF processing workflow.
"""

import logging
import click

logger = logging.getLogger(__name__)


def register_pipeline_command(group):
    """Register pipeline command with the HGC command group."""
    group.add_command(pipeline)


@click.command(name="pipeline")
@click.option(
    "-i",
    "--input-dir",
    type=click.Path(exists=True),
    required=True,
    help="Path to directory containing input gVCF files",
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    required=True,
    help="Path to output directory",
)
# Stage control flags
@click.option(
    "--import-interval-size",
    type=click.IntRange(min=1),
    default=None,
    help=(
        "Size (bp) of the even genomic intervals used to partition gVCF import in "
        "stage 1. Hail derives ONE PARTITION PER INTERVAL, so this caps the combine "
        "stage's parallelism. Default (Hail's genome default) is 1.2 Mb, which yields "
        "few partitions for a single chromosome (chr20 -> 54, chr1 -> 208), leaving "
        "extra cores idle. Lower it so partitions comfortably exceed your core count "
        "(~2-4x is a good target). Mutually exclusive with "
        "--use-exome-default-intervals."
    ),
)
@click.option(
    "--use-exome-default-intervals",
    is_flag=True,
    default=False,
    help="Partition gVCF import with Hail's exome default interval size (60 Mb).",
)
@click.option(
    "--gvcf-batch-size",
    type=click.IntRange(min=1),
    default=None,
    help="Number of gVCFs to combine per tree-merge batch (Hail default: 50).",
)
@click.option(
    "--branch-factor",
    type=click.IntRange(min=2),
    default=None,
    help="Branch factor of the combiner's hierarchical merge (Hail default: 100).",
)
@click.option(
    "--skip-combine-gvcfs",
    is_flag=True,
    default=False,
    help="Skip combining gVCF files (use existing VDS)",
)
@click.option(
    "--skip-vds-to-mt",
    is_flag=True,
    default=False,
    help="Skip VDS to MatrixTable conversion (use existing MT)",
)
@click.option(
    "--skip-compute-sample-qc",
    is_flag=True,
    default=False,
    help="Skip computing sample QC metrics",
)
@click.option(
    "--skip-compute-variant-qc",
    is_flag=True,
    default=False,
    help="Skip computing variant QC metrics",
)
@click.option(
    "--skip-export-pvcf",
    is_flag=True,
    default=False,
    help="Skip exporting the cohort (project) VCF",
)
@click.option(
    "--skip-validation",
    is_flag=True,
    default=False,
    help=(
        "Skip the biallelic audit and genotype repair during VDS -> MatrixTable conversion. "
        "The audit runs on the sparse variant data and is cheap; skip it only for a trusted, "
        "already-validated VDS."
    ),
)
# Path overrides
@click.option(
    "--vds-path",
    type=click.Path(),
    default=None,
    help="Path to existing VDS (required if --skip-combine-gvcfs)",
)
@click.option(
    "--mt-path",
    type=click.Path(),
    default=None,
    help="Path to existing MatrixTable (required if --skip-vds-to-mt)",
)
# Processing configuration
@click.option(
    "--tmp-dir",
    type=click.Path(),
    default=None,
    help="Path to temporary directory for intermediate files",
)
@click.option(
    "--reference-genome",
    type=click.Choice(["GRCh37", "GRCh38"]),
    default="GRCh38",
    help="Reference genome build",
)
@click.option(
    "--n-partitions",
    type=int,
    default=None,
    help=(
        "Coalesce the dense MatrixTable to this many partitions in the VDS -> MT stage. "
        "Reduces only. Default keeps the VDS's own layout, which is reference-block-derived "
        "and saturates as sample count grows, leaving partitions too thin to amortise task "
        "overhead (see #207). Size this from the dense matrix. Does not affect the gVCF "
        "combiner -- use --import-interval-size / --gvcf-batch-size / --branch-factor / "
        "--use-exome-default-intervals for stage 1."
    ),
)
@click.option(
    "--overwrite", is_flag=True, default=False, help="Overwrite existing output files"
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Show what would be done without executing",
)
# QC configuration
@click.option(
    "--min-sample-call-rate",
    type=float,
    default=0.85,
    help="Minimum sample call rate for QC filtering",
)
@click.option(
    "--min-variant-call-rate",
    type=float,
    default=0.85,
    help="Minimum variant call rate for QC filtering",
)
@click.option(
    "--apply-qc-filters",
    is_flag=True,
    default=False,
    help="Apply QC filters before exporting pVCF",
)
# Output options
@click.option(
    "--keep-intermediates",
    is_flag=True,
    default=True,
    help="Keep intermediate files (VDS, MT)",
)
@click.option(
    "--generate-qc-report",
    is_flag=True,
    default=False,
    help="Generate HTML QC report after computing QC metrics",
)
@click.option(
    "--output-prefix", type=str, default="cohort", help="Prefix for output files"
)
@click.pass_context
def pipeline(
    ctx,
    input_dir,
    output_dir,
    import_interval_size,
    use_exome_default_intervals,
    gvcf_batch_size,
    branch_factor,
    skip_combine_gvcfs,
    skip_vds_to_mt,
    skip_compute_sample_qc,
    skip_compute_variant_qc,
    skip_export_pvcf,
    skip_validation,
    vds_path,
    mt_path,
    tmp_dir,
    reference_genome,
    n_partitions,
    overwrite,
    dry_run,
    min_sample_call_rate,
    min_variant_call_rate,
    apply_qc_filters,
    keep_intermediates,
    generate_qc_report,
    output_prefix,
):
    """
    Run end-to-end gVCF processing pipeline.

    Orchestrates the complete workflow from gVCF files to cohort VCF:

    \b
    Stages:
      1. Combine gVCFs into VDS (Variant Dataset)
      2. Convert VDS to MatrixTable
      3. Compute sample QC metrics
      4. Compute variant QC metrics
      5. Export to project VCF (pVCF)

    \b
    Examples:
      # Run full pipeline
      hvantk hgc pipeline -i /path/to/gvcfs -o /path/to/output

      # Skip sample QC
      hvantk hgc pipeline -i /path/to/gvcfs -o /path/to/output --skip-compute-sample-qc

      # Start from existing VDS
      hvantk hgc pipeline -i /path/to/gvcfs -o /path/to/output \\
          --skip-combine-gvcfs --vds-path /path/to/existing.vds

      # Apply QC filters before export
      hvantk hgc pipeline -i /path/to/gvcfs -o /path/to/output \\
          --apply-qc-filters --min-sample-call-rate 0.9
    """
    try:
        from hvantk.algorithms.hgc.pipeline import PipelineConfig, PipelineRunner

        # Create configuration
        config = PipelineConfig(
            input_dir=input_dir,
            output_dir=output_dir,
            tmp_dir=tmp_dir,
            reference_genome=reference_genome,
            n_partitions=n_partitions,
            overwrite=overwrite,
            output_prefix=output_prefix,
            import_interval_size=import_interval_size,
            use_exome_default_intervals=use_exome_default_intervals,
            gvcf_batch_size=gvcf_batch_size,
            branch_factor=branch_factor,
            skip_combine_gvcfs=skip_combine_gvcfs,
            skip_vds_to_mt=skip_vds_to_mt,
            skip_compute_sample_qc=skip_compute_sample_qc,
            skip_compute_variant_qc=skip_compute_variant_qc,
            skip_export_pvcf=skip_export_pvcf,
            skip_validation=skip_validation,
            vds_path=vds_path,
            mt_path=mt_path,
            min_sample_call_rate=min_sample_call_rate,
            min_variant_call_rate=min_variant_call_rate,
            apply_qc_filters=apply_qc_filters,
            keep_intermediates=keep_intermediates,
            generate_qc_report=generate_qc_report,
        )

        # Validate configuration
        errors = config.validate()
        if errors:
            click.echo("❌ Configuration validation failed:", err=True)
            for error in errors:
                click.echo(f"   • {error}", err=True)
            ctx.exit(1)

        # Create pipeline runner
        runner = PipelineRunner(config)

        # Show plan if dry-run
        if dry_run:
            runner.show_plan()
            return

        # Display starting message
        click.echo("\n" + "=" * 70)
        click.echo("🚀 Starting HGC Pipeline Execution")
        click.echo("=" * 70 + "\n")

        # Run pipeline
        state = runner.run()

        # Display results
        click.echo("\n" + "=" * 70)
        if state.errors:
            click.echo("❌ Pipeline completed with errors:")
            for error in state.errors:
                click.echo(f"   • {error}")
            click.echo("=" * 70 + "\n")
            ctx.exit(1)
        else:
            click.echo("✅ Pipeline completed successfully!")
            click.echo("=" * 70)

            # Show outputs
            click.echo("\n📦 Output files:")
            for stage, output_path in state.outputs.items():
                click.echo(f"   • {stage}: {output_path}")

            click.echo(f"\n📊 Completed stages: {len(state.completed_stages)}")
            click.echo(f"⏱️  Duration: {state.start_time} → {state.end_time}")
            click.echo("")

    except Exception as e:
        logger.exception(f"Pipeline failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)
