"""
HGC CLI Commands Module

This module provides command-line interface commands for the HGC (Hail-based Genotype Combiner) toolkit.
Integrates with hvantk's CLI architecture to provide joint genotyping capabilities.

Commands:
- hvantk hgc gvcf-combine: Combine GVCF files for joint genotyping
- hvantk hgc vds-combine: Combine VDS datasets
- hvantk hgc vds2mt: Convert VDS to MatrixTable
- hvantk hgc mt2vcf: Convert MatrixTable to VCF
"""

import logging
import click
import glob
import os

logger = logging.getLogger(__name__)

# Import HGC functionality - use the actual functions that exist
from hvantk.hgc import (
    combine_gvcfs,
    combine_vdses,
    convert_vds_to_mt,
    convert_mt_to_multi_sample_vcf,
    check_path_exists_and_readable,
    validate_vds_paths,
)


# Utility functions for CLI
def setup_logging_for_hgc(log_level='INFO'):
    """Set up logging for HGC operations."""
    numeric_level = getattr(logging, log_level.upper(), logging.INFO)
    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )


def expand_file_patterns(patterns):
    """Expand file patterns/wildcards into actual file paths."""
    expanded = []
    for pattern in patterns:
        matches = glob.glob(pattern)
        expanded.extend(matches)
    return expanded


def validate_input_files(file_paths, file_type='gvcf'):
    """Validate input files and return (is_valid, errors) tuple."""
    errors = []
    try:
        if file_type == 'gvcf':
            for path in file_paths:
                check_path_exists_and_readable(path)
        elif file_type == 'vds':
            validate_vds_paths(file_paths)
        else:
            for path in file_paths:
                check_path_exists_and_readable(path)
        return (True, [])
    except Exception as e:
        errors.append(str(e))
        return (False, errors)


def validate_output_path(output_path, create_dirs=False):
    """Validate output path and optionally create parent directories."""
    try:
        parent_dir = os.path.dirname(output_path)
        if parent_dir and not os.path.exists(parent_dir):
            if create_dirs:
                os.makedirs(parent_dir, exist_ok=True)
            else:
                return False
        return True
    except Exception:
        return False


def estimate_resource_requirements(file_paths, operation='combine'):
    """Estimate resource requirements for an operation."""
    total_size = 0
    for path in file_paths:
        try:
            if os.path.isfile(path):
                total_size += os.path.getsize(path)
            elif os.path.isdir(path):
                for root, dirs, files in os.walk(path):
                    for file in files:
                        total_size += os.path.getsize(os.path.join(root, file))
        except Exception:
            pass

    total_size_gb = total_size / (1024 ** 3)

    # Simple heuristic estimates
    memory = f"{max(4, int(total_size_gb * 2))}g"
    partitions = max(100, int(total_size_gb * 10))
    estimated_runtime = max(5, int(total_size_gb * 2))

    return {
        'memory': memory,
        'partitions': partitions,
        'estimated_runtime_minutes': estimated_runtime,
        'total_size_gb': round(total_size_gb, 2)
    }


@click.group(name="hgc", help="HGC (Hail-based Genotype Combiner) commands for joint genotyping workflows")
@click.option('--log-level', default='INFO',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              help='Set logging level')
@click.pass_context
def hgc_group(ctx, log_level):
    """
    HGC command group for joint genotyping operations.
    """
    # Ensure that ctx.obj exists and is a dict
    ctx.ensure_object(dict)
    ctx.obj['log_level'] = log_level

    # Set up logging for HGC operations
    setup_logging_for_hgc(log_level)
    logger.info(f"Starting HGC command with log level: {log_level}")


@hgc_group.command(name="gvcf-combine")
@click.option('--gvcf-dir', '-g', help='Directory containing GVCF files')
@click.option('--vds-paths', '-v', multiple=True, help='VDS paths to combine with GVCFs')
@click.option('--output', '-o', required=True, help='Output VDS path')
@click.option('--temp-dir', '--tmp', default='/tmp/hgc', help='Temporary directory for intermediate files')
@click.option('--save-path', '-s', help='Path to save the combiner plan')
@click.option('--dry-run', is_flag=True, help='Show what would be done without executing')
@click.pass_context
def gvcf_combine(ctx, gvcf_dir, vds_paths, output, temp_dir, save_path, dry_run):
    """
    Combine GVCF files for joint genotyping.

    This command takes GVCF files from a directory and/or existing VDS datasets
    and combines them into a single joint-called variant dataset using Hail's optimized combiner.

    Examples:
        hvantk hgc gvcf-combine -g /path/to/gvcfs -o combined.vds
        hvantk hgc gvcf-combine -g /path/to/gvcfs -v existing.vds -o combined.vds
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

        if dry_run:
            click.echo("🔍 Dry run mode - would execute GVCF combination with:")
            click.echo(f"   • GVCF directory: {gvcf_dir or 'None'}")
            click.echo(f"   • VDS paths: {list(vds_paths) if vds_paths else 'None'}")
            click.echo(f"   • Output: {output}")
            click.echo(f"   • Temp directory: {temp_dir}")
            return

        # Execute combination
        click.echo("🔄 Starting GVCF combination...")
        combine_gvcfs(
            gvcf_dir=gvcf_dir,
            vds_output_path=output,
            tmp_path=temp_dir,
            save_path=save_path or f"{output}.plan",
            vdses=list(vds_paths) if vds_paths else [],
            kwargs={}
        )

        click.echo(f"✅ Successfully combined GVCFs to {output}")

    except Exception as e:
        logger.error(f"GVCF combination failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)


@hgc_group.command(name="vds-combine")
@click.option('--input-dir', '-i', required=True, help='Directory containing VDS datasets')
@click.option('--output', '-o', required=True, help='Output path for combined VDS')
@click.option('--validate/--no-validate', default=True, help='Validate the combined VDS')
@click.option('--overwrite/--no-overwrite', default=False, help='Overwrite output if exists')
@click.option('--dry-run', is_flag=True, help='Show what would be done without executing')
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
            overwrite=overwrite
        )

        click.echo(f"✅ Successfully combined VDS datasets to {output}")

    except Exception as e:
        logger.error(f"VDS combination failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)


@hgc_group.command(name="vds2mt")
@click.option('--input', '-i', required=True, help='Input VDS dataset')
@click.option('--output', '-o', required=True, help='Output MatrixTable path')
@click.option('--adjust-genotypes/--no-adjust-genotypes', default=True,
              help='Annotate with adjusted genotypes')
@click.option('--skip-split-multi', is_flag=True, help='Skip splitting multi-allelic variants')
@click.option('--convert-lgt-to-gt/--no-convert-lgt-to-gt', default=True,
              help='Convert LGT to GT (recommended)')
@click.option('--skip-keying-by-cols', is_flag=True, help='Skip keying MatrixTable by columns')
@click.option('--overwrite/--no-overwrite', default=False, help='Overwrite output if exists')
@click.option('--dry-run', is_flag=True, help='Show what would be done without executing')
@click.pass_context
def vds2mt(ctx, input, output, adjust_genotypes, skip_split_multi, convert_lgt_to_gt,
           skip_keying_by_cols, overwrite, dry_run):
    """
    Convert Variant DataSet (VDS) to MatrixTable format.

    VDS is optimized for storage while MatrixTable is better for analysis.
    This conversion creates a dense matrix which is recommended for most analyses.

    Examples:
        hvantk hgc vds2mt -i dataset.vds -o dataset.mt
        hvantk hgc vds2mt -i dataset.vds -o dataset.mt --skip-split-multi
    """
    try:
        logger.info("Starting VDS to MatrixTable conversion")

        # Validate input
        is_valid, errors = validate_input_files([input], 'vds')
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
            click.echo("🔍 Dry run mode - would execute VDS to MatrixTable conversion with:")
            click.echo(f"   • Input: {input}")
            click.echo(f"   • Output: {output}")
            click.echo(f"   • Adjust genotypes: {adjust_genotypes}")
            click.echo(f"   • Skip split multi: {skip_split_multi}")
            click.echo(f"   • Convert LGT to GT: {convert_lgt_to_gt}")
            return

        # Execute conversion
        click.echo("🔄 Starting VDS to MatrixTable conversion...")
        convert_vds_to_mt(
            vds_path=input,
            output_path=output,
            adjust_genotypes=adjust_genotypes,
            skip_split_multi=skip_split_multi,
            convert_lgt_to_gt=convert_lgt_to_gt,
            skip_keying_by_cols=skip_keying_by_cols,
            overwrite=overwrite
        )

        click.echo(f"✅ Successfully converted {input} to MatrixTable at {output}")

    except Exception as e:
        logger.error(f"VDS to MatrixTable conversion failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)


@hgc_group.command(name="mt2vcf")
@click.option('--input', '-i', required=True, help='Input MatrixTable')
@click.option('--output', '-o', required=True, help='Output VCF file path')
@click.option('--filter-adj/--no-filter-adj', default=True,
              help='Filter to adjusted genotypes (recommended)')
@click.option('--min-ac', default=1, type=int, help='Minimum alternate allele count')
@click.option('--split-multi/--no-split-multi', default=True, help='Split multi-allelic variants')
@click.option('--dry-run', is_flag=True, help='Show what would be done without executing')
@click.pass_context
def mt2vcf(ctx, input, output, filter_adj, min_ac, split_multi, dry_run):
    """
    Convert MatrixTable to VCF format.

    Exports processed genomic data back to standard VCF format for
    compatibility with other tools and pipelines.

    Examples:
        hvantk hgc mt2vcf -i analysis.mt -o results.vcf.gz
        hvantk hgc mt2vcf -i analysis.mt -o results.vcf.gz --min-ac 2
    """
    try:
        logger.info("Starting MatrixTable to VCF conversion")

        # Validate input
        is_valid, errors = validate_input_files([input], 'mt')
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
            click.echo("🔍 Dry run mode - would execute MatrixTable to VCF conversion with:")
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
            split_multi=split_multi
        )

        click.echo(f"✅ Successfully converted {input} to VCF at {output}")

    except Exception as e:
        logger.error(f"MatrixTable to VCF conversion failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)

