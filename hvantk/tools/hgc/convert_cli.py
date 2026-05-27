"""
HGC CLI - Convert Commands

Commands for converting between different variant data formats.
"""

import logging
import click

from hvantk.algorithms.hgc import convert_vds_to_mt, convert_mt_to_multi_sample_vcf
from .utils import validate_input_files, validate_output_path

logger = logging.getLogger(__name__)


def register_convert_commands(group):
    """Register convert commands with the HGC command group."""
    group.add_command(vds2mt)
    group.add_command(mt2vcf)


@click.command(name="vds2mt")
@click.option("--input", "-i", required=True, help="Input VDS dataset")
@click.option("--output", "-o", required=True, help="Output MatrixTable path")
@click.option(
    "--adjust-genotypes/--no-adjust-genotypes",
    default=True,
    help="Annotate with adjusted genotypes",
)
@click.option(
    "--skip-split-multi", is_flag=True, help="Skip splitting multi-allelic variants"
)
@click.option(
    "--skip-validation",
    is_flag=True,
    help="Skip biallelic validation (faster, use only if confident)",
)
@click.option(
    "--skip-keying-by-cols", is_flag=True, help="Skip keying MatrixTable by columns"
)
@click.option(
    "--overwrite/--no-overwrite", default=False, help="Overwrite output if exists"
)
@click.option(
    "--dry-run", is_flag=True, help="Show what would be done without executing"
)
@click.pass_context
def vds2mt(
    ctx,
    input,
    output,
    adjust_genotypes,
    skip_split_multi,
    skip_validation,
    skip_keying_by_cols,
    overwrite,
    dry_run,
):
    """
    Convert Variant DataSet (VDS) to MatrixTable format.

    VDS is optimized for storage while MatrixTable is better for analysis.
    This conversion creates a dense matrix which is recommended for most analyses.

    Examples:
        hvantk hgc vds2mt -i dataset.vds -o dataset.mt
        hvantk hgc vds2mt -i dataset.vds -o dataset.mt --skip-validation
    """
    try:
        logger.info("Starting VDS to MatrixTable conversion")

        # Validate input
        is_valid, errors = validate_input_files([input], "vds")
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
            click.echo(
                "🔍 Dry run mode - would execute VDS to MatrixTable conversion with:"
            )
            click.echo(f"   • Input: {input}")
            click.echo(f"   • Output: {output}")
            click.echo(f"   • Adjust genotypes: {adjust_genotypes}")
            click.echo(f"   • Skip split multi: {skip_split_multi}")
            click.echo(f"   • Skip validation: {skip_validation}")
            return

        # Execute conversion
        click.echo("🔄 Starting VDS to MatrixTable conversion...")
        convert_vds_to_mt(
            vds_path=input,
            output_path=output,
            adjust_genotypes=adjust_genotypes,
            skip_split_multi=skip_split_multi,
            skip_validation=skip_validation,
            skip_keying_by_cols=skip_keying_by_cols,
            overwrite=overwrite,
        )

        click.echo(f"✅ Successfully converted {input} to MatrixTable at {output}")

    except Exception as e:
        logger.exception(f"VDS to MatrixTable conversion failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)


@click.command(name="mt2vcf")
@click.option("--input", "-i", required=True, help="Input MatrixTable")
@click.option("--output", "-o", required=True, help="Output VCF file path")
@click.option(
    "--filter-adj/--no-filter-adj",
    default=True,
    help="Filter to adjusted genotypes (recommended)",
)
@click.option("--min-ac", default=1, type=int, help="Minimum alternate allele count")
@click.option(
    "--split-multi/--no-split-multi", default=True, help="Split multi-allelic variants"
)
@click.option(
    "--dry-run", is_flag=True, help="Show what would be done without executing"
)
@click.pass_context
def mt2vcf(ctx, input, output, filter_adj, min_ac, split_multi, dry_run):
    """
    Convert MatrixTable to VCF format.

    Exports processed genomic data back to standard VCF format for
    compatibility with other tools and pipelines.

    Examples:
        hvantk hgc mt2vcf -i analysis.mt -o results.vcf.bgz
        hvantk hgc mt2vcf -i analysis.mt -o results.vcf.bgz --min-ac 2
    """
    try:
        logger.info("Starting MatrixTable to VCF conversion")

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
            click.echo(
                "🔍 Dry run mode - would execute MatrixTable to VCF conversion with:"
            )
            click.echo(f"   • Input: {input}")
            click.echo(f"   • Output: {output}")
            click.echo(f"   • Filter adjusted genotypes: {filter_adj}")
            click.echo(f"   • Minimum AC: {min_ac}")
            click.echo(f"   • Split multi: {split_multi}")
            return

        # Execute conversion
        click.echo("🔄 Starting MatrixTable to VCF conversion...")
        convert_mt_to_multi_sample_vcf(
            mt_path=input,
            vcf_path=output,
            filter_adj_genotypes=filter_adj,
            min_ac=min_ac,
            split_multi=split_multi,
        )

        click.echo(f"✅ Successfully converted {input} to VCF at {output}")

    except Exception as e:
        logger.exception(f"MatrixTable to VCF conversion failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)
