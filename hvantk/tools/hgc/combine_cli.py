"""
HGC CLI - Combine Commands

Commands for combining GVCF files and VDS datasets.
"""

import logging
import os
import click

from hvantk.algorithms.hgc import combine_gvcfs, combine_vdses
from .utils import validate_output_path, DEFAULT_TEMP_DIR

logger = logging.getLogger(__name__)


def register_combine_commands(group):
    """Register combine commands with the HGC command group."""
    group.add_command(gvcf_combine)
    group.add_command(vds_combine)


@click.command(name="gvcf-combine")
@click.option("--gvcf-dir", "-g", help="Directory containing GVCF files")
@click.option(
    "--vds-paths", "-v", multiple=True, help="VDS paths to combine with GVCFs"
)
@click.option("--output", "-o", required=True, help="Output VDS path")
@click.option(
    "--temp-dir",
    "--tmp",
    default=DEFAULT_TEMP_DIR,
    help="Temporary directory for intermediate files",
)
@click.option("--save-path", "-s", help="Path to save the combiner plan")
@click.option(
    "--import-interval-size",
    type=int,
    default=None,
    help=(
        "Size (bp) of the even genomic intervals used to partition GVCF import. "
        "Hail derives ONE PARTITION PER INTERVAL, so this caps the combiner's "
        "parallelism. Default (Hail's genome default) is 1.2 Mb. Lower it so that "
        "partitions comfortably exceed your core count (~2-4x is a good target). "
        "Mutually exclusive with --use-genome-default-intervals / "
        "--use-exome-default-intervals."
    ),
)
@click.option(
    "--use-genome-default-intervals",
    is_flag=True,
    default=False,
    help="Partition with Hail's genome default interval size (1.2 Mb). This is the default behaviour.",
)
@click.option(
    "--use-exome-default-intervals",
    is_flag=True,
    default=False,
    help="Partition with Hail's exome default interval size (60 Mb).",
)
@click.option(
    "--gvcf-batch-size",
    type=int,
    default=None,
    help="Number of GVCFs to combine per tree-merge batch (Hail default: 50).",
)
@click.option(
    "--branch-factor",
    type=int,
    default=None,
    help="Branch factor of the combiner's hierarchical merge (Hail default: 100).",
)
@click.option(
    "--dry-run", is_flag=True, help="Show what would be done without executing"
)
@click.pass_context
def gvcf_combine(
    ctx,
    gvcf_dir,
    vds_paths,
    output,
    temp_dir,
    save_path,
    import_interval_size,
    use_genome_default_intervals,
    use_exome_default_intervals,
    gvcf_batch_size,
    branch_factor,
    dry_run,
):
    """
    Combine GVCF files for joint genotyping.

    This command takes GVCF files from a directory and/or existing VDS datasets
    and combines them into a single joint-called variant dataset using Hail's optimized combiner.

    Tuning parallelism: the combiner creates one partition per genomic interval, so the
    interval size sets the ceiling on how many cores can do useful work. With the 1.2 Mb
    genome default, a single chromosome yields relatively few partitions (chr20 -> 54,
    chr1 -> 208) and any cores beyond that idle. Use --import-interval-size to raise the
    partition count above your core count.

    Examples:
        hvantk hgc gvcf-combine -g /path/to/gvcfs -o combined.vds
        hvantk hgc gvcf-combine -g /path/to/gvcfs -v existing.vds -o combined.vds
        hvantk hgc gvcf-combine -g /path/to/gvcfs -o combined.vds --import-interval-size 600000
    """
    try:
        logger.info("Starting GVCF combination workflow")

        if not gvcf_dir and not vds_paths:
            click.echo("❌ Either --gvcf-dir or --vds-paths must be provided", err=True)
            ctx.exit(1)

        # Validate output path
        if not validate_output_path(output, create_dirs=True):
            click.echo("❌ Invalid output path", err=True)
            ctx.exit(1)

        # Exactly one partitioning strategy may be chosen (Hail only warns and silently
        # picks one when they collide, so reject it here).
        partition_opts = [
            import_interval_size is not None,
            use_genome_default_intervals,
            use_exome_default_intervals,
        ]
        if sum(partition_opts) > 1:
            click.echo(
                "❌ --import-interval-size, --use-genome-default-intervals and "
                "--use-exome-default-intervals are mutually exclusive",
                err=True,
            )
            ctx.exit(1)

        # Only forward keys the user actually set; an empty dict preserves the previous
        # behaviour (Hail's 1.2 Mb genome default).
        combiner_kwargs = {}
        if import_interval_size is not None:
            combiner_kwargs["import_interval_size"] = import_interval_size
        if use_genome_default_intervals:
            combiner_kwargs["use_genome_default_intervals"] = True
        if use_exome_default_intervals:
            combiner_kwargs["use_exome_default_intervals"] = True
        if gvcf_batch_size is not None:
            combiner_kwargs["gvcf_batch_size"] = gvcf_batch_size
        if branch_factor is not None:
            combiner_kwargs["branch_factor"] = branch_factor

        if dry_run:
            click.echo("🔍 Dry run mode - would execute GVCF combination with:")
            click.echo(f"   • GVCF directory: {gvcf_dir or 'None'}")
            click.echo(f"   • VDS paths: {list(vds_paths) if vds_paths else 'None'}")
            click.echo(f"   • Output: {output}")
            click.echo(f"   • Temp directory: {temp_dir}")
            click.echo(f"   • Combiner options: {combiner_kwargs or 'defaults'}")
            return

        # Execute combination
        click.echo("🔄 Starting GVCF combination...")
        combine_gvcfs(
            gvcf_dir=gvcf_dir,
            vds_output_path=output,
            tmp_path=temp_dir,
            save_path=save_path or f"{output}.plan",
            vdses=list(vds_paths) if vds_paths else [],
            kwargs=combiner_kwargs,
        )

        click.echo(f"✅ Successfully combined GVCFs to {output}")

    except Exception as e:
        logger.exception(f"GVCF combination failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)


@click.command(name="vds-combine")
@click.option(
    "--input-dir", "-i", required=True, help="Directory containing VDS datasets"
)
@click.option("--output", "-o", required=True, help="Output path for combined VDS")
@click.option(
    "--validate/--no-validate", default=True, help="Validate the combined VDS"
)
@click.option(
    "--overwrite/--no-overwrite", default=False, help="Overwrite output if exists"
)
@click.option(
    "--dry-run", is_flag=True, help="Show what would be done without executing"
)
@click.pass_context
def vds_combine(ctx, input_dir, output, validate, overwrite, dry_run):
    """
    Combine Variant DataSets (VDS) for joint analysis.

    VDS is Hail's optimized format for large-scale variant data storage
    and processing with efficient compression and query capabilities.

    Examples:
        hvantk hgc vds-combine -i /path/to/vds_datasets -o combined.vds
        hvantk hgc vds-combine -i /path/to/vds_datasets -o combined.vds --no-validate
    """
    try:
        logger.info("Starting VDS combination workflow")

        # Validate input directory
        if not os.path.isdir(input_dir):
            click.echo(f"❌ Input directory not found: {input_dir}", err=True)
            ctx.exit(1)

        # Validate output path
        if not validate_output_path(output, create_dirs=True):
            click.echo("❌ Invalid output path", err=True)
            ctx.exit(1)

        if dry_run:
            click.echo("🔍 Dry run mode - would execute VDS combination with:")
            click.echo(f"   • Input directory: {input_dir}")
            click.echo(f"   • Output: {output}")
            click.echo(f"   • Validate: {validate}")
            click.echo(f"   • Overwrite: {overwrite}")
            return

        # Execute combination
        click.echo("🔄 Starting VDS combination...")
        combine_vdses(
            vdses_dir=input_dir,
            output_path=output,
            validate=validate,
            overwrite=overwrite,
        )

        click.echo(f"✅ Successfully combined VDS datasets to {output}")

    except Exception as e:
        logger.exception(f"VDS combination failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)
