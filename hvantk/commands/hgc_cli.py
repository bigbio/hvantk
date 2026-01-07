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
import tempfile

logger = logging.getLogger(__name__)

# Import HGC functionality - use the actual functions that exist
from hvantk.hgc import (
    combine_gvcfs,
    combine_vdses,
    convert_vds_to_mt,
    convert_mt_to_multi_sample_vcf,
    check_path_exists_and_readable,
    validate_vds_paths,
    compute_full_qc,
    compute_sample_qc,
    compute_variant_qc,
    filter_samples_by_qc,
    filter_variants_by_qc,
    save_qc_metrics
)

DEFAULT_TEMP_DIR = os.environ.get('HGC_TEMP_DIR', tempfile.gettempdir())


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


def estimate_resource_requirements(file_paths):
    """Estimate resource requirements for an operation."""
    logger = logging.getLogger(__name__)
    total_size = 0
    for path in file_paths:
        try:
            if os.path.isfile(path):
                total_size += os.path.getsize(path)
            elif os.path.isdir(path):
                for root, dirs, files in os.walk(path):
                    for file in files:
                        file_path = os.path.join(root, file)
                        try:
                            total_size += os.path.getsize(file_path)
                        except OSError as e:
                            logger.warning(
                                f"Failed to get size of file '{file_path}': {e}"
                            )
                            continue
        except OSError as e:
            logger.warning(
                f"Failed to access path '{path}': {e}"
            )
            continue
        except Exception as e:
            logger.error(
                f"Unexpected error accessing path '{path}': {e}"
            )
            raise

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
@click.option('--temp-dir', '--tmp', default=DEFAULT_TEMP_DIR, help='Temporary directory for intermediate files')
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
        logger.exception(f"GVCF combination failed: {e}")
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
        logger.exception(f"VDS combination failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)


@hgc_group.command(name="vds2mt")
@click.option('--input', '-i', required=True, help='Input VDS dataset')
@click.option('--output', '-o', required=True, help='Output MatrixTable path')
@click.option('--adjust-genotypes/--no-adjust-genotypes', default=True,
              help='Annotate with adjusted genotypes')
@click.option('--skip-split-multi', is_flag=True, help='Skip splitting multi-allelic variants')
@click.option('--skip-validation', is_flag=True,
              help='Skip biallelic validation (faster, use only if confident)')
@click.option('--skip-keying-by-cols', is_flag=True, help='Skip keying MatrixTable by columns')
@click.option('--overwrite/--no-overwrite', default=False, help='Overwrite output if exists')
@click.option('--dry-run', is_flag=True, help='Show what would be done without executing')
@click.pass_context
def vds2mt(ctx, input, output, adjust_genotypes, skip_split_multi, skip_validation,
           skip_keying_by_cols, overwrite, dry_run):
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
            overwrite=overwrite
        )

        click.echo(f"✅ Successfully converted {input} to MatrixTable at {output}")

    except Exception as e:
        logger.exception(f"VDS to MatrixTable conversion failed: {e}")
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
        hvantk hgc mt2vcf -i analysis.mt -o results.vcf.bgz
        hvantk hgc mt2vcf -i analysis.mt -o results.vcf.bgz --min-ac 2
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
        logger.exception(f"MatrixTable to VCF conversion failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)

@hgc_group.command(name="compute-qc")
@click.option('--input', '-i', required=True, help='Input MatrixTable path')
@click.option('--output-dir', '-o', required=True, help='Output directory for QC metrics')
@click.option('--sample-qc/--no-sample-qc', default=True, help='Compute sample-level QC metrics')
@click.option('--variant-qc/--no-variant-qc', default=True, help='Compute variant-level QC metrics')
@click.option('--call-field', default='GT', help='Name of the call field to use (default: GT)')
@click.option('--prefix', default='qc_metrics', help='Prefix for output files')
@click.option('--save-mt/--no-save-mt', default=True, help='Save MatrixTable with QC annotations')
@click.option('--dry-run', is_flag=True, help='Show what would be done without executing')
@click.pass_context
def compute_qc(ctx, input, output_dir, sample_qc, variant_qc, call_field, prefix, save_mt, dry_run):
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
        is_valid, errors = validate_input_files([input], 'mt')
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
            click.echo("❌ At least one of --sample-qc or --variant-qc must be enabled", err=True)
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

        # Import hail and read MatrixTable
        import hail as hl
        hl.init(quiet=True)

        click.echo("🔄 Loading MatrixTable...")
        mt = hl.read_matrix_table(input)

        # Compute QC metrics based on options
        if sample_qc and variant_qc:
            click.echo("🔄 Computing comprehensive QC metrics...")
            qc_results = compute_full_qc(mt, call_field=call_field)
        elif sample_qc:
            click.echo("🔄 Computing sample QC metrics...")
            mt_qc = compute_sample_qc(mt, call_field=call_field)
            from hvantk.hgc.qc import QCMetrics
            sample_qc_table = mt_qc.cols().select('sample_qc')
            qc_results = QCMetrics(mt_qc, sample_qc_table, None)
        else:  # variant_qc only
            click.echo("🔄 Computing variant QC metrics...")
            mt_qc = compute_variant_qc(mt, call_field=call_field)
            from hvantk.hgc.qc import QCMetrics
            variant_qc_table = mt_qc.rows().select('variant_qc')
            qc_results = QCMetrics(mt_qc, None, variant_qc_table)

        # Save QC metrics
        click.echo("💾 Saving QC metrics...")
        saved_files = save_qc_metrics(qc_results, output_dir, prefix)

        # Remove MatrixTable from saved files if not requested
        if not save_mt and 'matrix_table' in saved_files:
            import os
            import shutil
            if os.path.exists(saved_files['matrix_table']):
                shutil.rmtree(saved_files['matrix_table'])
            del saved_files['matrix_table']

        click.echo("✅ Successfully computed and saved QC metrics:")
        for file_type, file_path in saved_files.items():
            click.echo(f"   • {file_type}: {file_path}")

    except Exception as e:
        logger.exception(f"QC computation failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)


@hgc_group.command(name="filter-qc")
@click.option('--input', '-i', required=True, help='Input MatrixTable with QC annotations')
@click.option('--output', '-o', required=True, help='Output filtered MatrixTable path')
@click.option('--min-sample-call-rate', type=float, default=0.85, help='Minimum sample call rate')
@click.option('--min-variant-call-rate', type=float, default=0.85, help='Minimum variant call rate')
@click.option('--min-ac', type=int, default=1, help='Minimum allele count for variants')
@click.option('--max-ac', type=int, help='Maximum allele count for variants')
@click.option('--min-af', type=float, help='Minimum allele frequency for variants')
@click.option('--max-af', type=float, help='Maximum allele frequency for variants')
@click.option('--hwe-threshold', type=float, default=1e-6, help='Hardy-Weinberg equilibrium p-value threshold')
@click.option('--min-mean-dp', type=float, help='Minimum mean depth for samples')
@click.option('--max-mean-dp', type=float, help='Maximum mean depth for samples')
@click.option('--min-mean-gq', type=float, help='Minimum mean genotype quality for samples')
@click.option('--sample-qc-name', default='sample_qc', help='Name of sample QC annotation')
@click.option('--variant-qc-name', default='variant_qc', help='Name of variant QC annotation')
@click.option('--overwrite/--no-overwrite', default=False, help='Overwrite output if exists')
@click.option('--dry-run', is_flag=True, help='Show what would be done without executing')
@click.pass_context
def filter_qc(ctx, input, output, min_sample_call_rate, min_variant_call_rate,
              min_ac, max_ac, min_af, max_af, hwe_threshold, min_mean_dp, max_mean_dp,
              min_mean_gq, sample_qc_name, variant_qc_name, overwrite, dry_run):
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
            click.echo("🔍 Dry run mode - would execute QC filtering with:")
            click.echo(f"   • Input: {input}")
            click.echo(f"   • Output: {output}")
            click.echo(f"   • Sample filters:")
            click.echo(f"     - Min call rate: {min_sample_call_rate}")
            if min_mean_dp: click.echo(f"     - Min mean depth: {min_mean_dp}")
            if max_mean_dp: click.echo(f"     - Max mean depth: {max_mean_dp}")
            if min_mean_gq: click.echo(f"     - Min mean GQ: {min_mean_gq}")
            click.echo(f"   • Variant filters:")
            click.echo(f"     - Min call rate: {min_variant_call_rate}")
            click.echo(f"     - Min AC: {min_ac}")
            if max_ac: click.echo(f"     - Max AC: {max_ac}")
            if min_af: click.echo(f"     - Min AF: {min_af}")
            if max_af: click.echo(f"     - Max AF: {max_af}")
            click.echo(f"     - HWE threshold: {hwe_threshold}")
            return

        # Import hail and read MatrixTable
        import hail as hl
        hl.init(quiet=True)

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
                sample_qc_name=sample_qc_name
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
                variant_qc_name=variant_qc_name
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
        click.echo(f"   • Samples: {n_samples_initial} → {n_samples_final} "
                  f"({samples_removed} removed, {samples_removed/n_samples_initial*100:.1f}%)")
        click.echo(f"   • Variants: {n_variants_initial} → {n_variants_final} "
                  f"({variants_removed} removed, {variants_removed/n_variants_initial*100:.1f}%)")
        click.echo(f"   • Output: {output}")

    except Exception as e:
        logger.exception(f"QC filtering failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)


@hgc_group.command(name="qc-summary")
@click.option('--qc-dir', '-d', required=True, help='Directory containing QC metrics files')
@click.option('--sample-file', help='Specific sample QC file (optional)')
@click.option('--variant-file', help='Specific variant QC file (optional)')
@click.option('--output', '-o', help='Output file for summary report')
@click.option('--format', type=click.Choice(['csv', 'json', 'markdown']), default='markdown',
              help='Output format for summary')
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
                    click.echo(f"ℹ️  Multiple sample QC files found, using {sample_path.name}")

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
                    click.echo(f"ℹ️  Multiple variant QC files found, using {variant_path.name}")

        if not sample_path and not variant_path:
            click.echo("❌ No QC files found", err=True)
            ctx.exit(1)

        # Generate summary
        summary_data = {}

        if sample_path and sample_path.exists():
            click.echo(f"📊 Processing sample QC metrics from {sample_path}")
            sample_df = pd.read_csv(sample_path)
            from hvantk.hgc.qc import get_qc_summary_stats
            sample_summary = get_qc_summary_stats(sample_df)
            summary_data['sample_qc'] = {
                'file': str(sample_path),
                'n_samples': len(sample_df),
                'summary_stats': sample_summary.to_dict() if not sample_summary.empty else {},
                'columns': sample_df.columns.tolist()
            }

        if variant_path and variant_path.exists():
            click.echo(f"📊 Processing variant QC metrics from {variant_path}")
            variant_df = pd.read_csv(variant_path)
            from hvantk.hgc.qc import get_qc_summary_stats
            variant_summary = get_qc_summary_stats(variant_df)
            summary_data['variant_qc'] = {
                'file': str(variant_path),
                'n_variants': len(variant_df),
                'summary_stats': variant_summary.to_dict() if not variant_summary.empty else {},
                'columns': variant_df.columns.tolist()
            }

        # Output summary
        if output:
            output_path = Path(output)
            if format == 'csv':
                # Combine summaries into CSV
                all_stats = []
                for qc_type, data in summary_data.items():
                    for metric, stats in data['summary_stats'].items():
                        for stat_name, value in stats.items():
                            all_stats.append({
                                'qc_type': qc_type,
                                'metric': metric,
                                'statistic': stat_name,
                                'value': value
                            })
                summary_df = pd.DataFrame(all_stats)
                summary_df.to_csv(output_path, index=False)

            elif format == 'json':
                import json
                with open(output_path, 'w') as f:
                    json.dump(summary_data, f, indent=2)

            elif format == 'markdown':
                # Generate markdown report
                md_content = "# Quality Control Summary Report\n\n"

                for qc_type, data in summary_data.items():
                    md_content += f"## {qc_type.replace('_', ' ').title()}\n\n"
                    md_content += f"- **File**: {data['file']}\n"
                    md_content += f"- **Count**: {data.get('n_samples', data.get('n_variants', 0)):,}\n"
                    md_content += f"- **Metrics**: {', '.join(data['columns'])}\n\n"

                    if data['summary_stats']:
                        md_content += "### Summary Statistics\n\n"
                        md_content += "| Metric | Count | Mean | Std | Min | 25% | 50% | 75% | Max |\n"
                        md_content += "|--------|--------|--------|--------|--------|--------|--------|--------|--------|\n"

                        for metric, stats in data['summary_stats'].items():
                            if isinstance(stats, dict):
                                row = f"| {metric} |"
                                for stat in ['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']:
                                    value = stats.get(stat, 'N/A')
                                    if isinstance(value, (int, float)) and stat != 'count':
                                        value = f"{value:.3f}" if abs(value) < 1000 else f"{value:.2e}"
                                    row += f" {value} |"
                                md_content += row + "\n"
                        md_content += "\n"

                with open(output_path, 'w') as f:
                    f.write(md_content)

            click.echo(f"💾 Summary saved to {output_path}")
        else:
            # Print to console
            click.echo("\n📊 QC Summary Statistics:")
            for qc_type, data in summary_data.items():
                click.echo(f"\n{qc_type.replace('_', ' ').title()}:")
                click.echo(f"  File: {data['file']}")
                click.echo(f"  Count: {data.get('n_samples', data.get('n_variants', 0)):,}")
                click.echo(f"  Metrics: {len(data['columns'])}")

        click.echo("✅ QC summary completed successfully")

    except Exception as e:
        logger.exception(f"QC summary generation failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)


@hgc_group.command(name="plot-qc")
@click.option('--input', '-i', required=True, help='Input MatrixTable with QC annotations')
@click.option('--output-dir', '-o', required=True, help='Output directory for plots')
@click.option('--plot-type',
              type=click.Choice(['overview', 'individual', 'dashboard', 'all']),
              default='overview',
              help='Type of plots to generate')
@click.option('--format', 'output_format',
              type=click.Choice(['png', 'pdf', 'svg']),
              default='png',
              help='Output format for plots')
@click.option('--style',
              type=click.Choice(['default', 'publication']),
              default='default',
              help='Plot style')
@click.option('--figsize', default='10,6', help='Figure size as width,height (e.g., 10,6)')
@click.option('--dpi', default=300, help='Resolution for output plots')
@click.option('--interactive', is_flag=True, help='Generate interactive plots (requires plotly)')
@click.option('--dry-run', is_flag=True, help='Show what would be done without executing')
@click.pass_context
def plot_qc(ctx, input, output_dir, plot_type, output_format, style, figsize, dpi, interactive, dry_run):
    """
    Generate QC plots from MatrixTable with QC annotations.

    Create individual plots, overviews, or comprehensive dashboards from quality control
    metrics computed on genomic variant data.

    Examples:
        hvantk hgc plot-qc -i cohort_qc.mt -o plots/ --plot-type overview
        hvantk hgc plot-qc -i cohort_qc.mt -o plots/ --plot-type dashboard --style publication
        hvantk hgc plot-qc -i cohort_qc.mt -o plots/ --plot-type all --format pdf
        hvantk hgc plot-qc -i cohort_qc.mt -o plots/ --interactive --dry-run
    """
    try:
        import hail as hl
        from pathlib import Path

        # Validate inputs
        if not check_path_exists_and_readable(input):
            click.echo(f"❌ Input file not found or not readable: {input}", err=True)
            ctx.exit(1)

        # Parse figsize
        try:
            width, height = map(float, figsize.split(','))
            figsize_tuple = (width, height)
        except ValueError:
            click.echo(f"❌ Invalid figsize format: {figsize}. Use format: width,height", err=True)
            ctx.exit(1)

        output_path = Path(output_dir)

        if dry_run:
            click.echo("🔍 Dry run mode - showing what would be done:")
            click.echo(f"   • Input MatrixTable: {input}")
            click.echo(f"   • Output directory: {output_path.absolute()}")
            click.echo(f"   • Plot type: {plot_type}")
            click.echo(f"   • Output format: {output_format}")
            click.echo(f"   • Style: {style}")
            click.echo(f"   • Figure size: {figsize_tuple}")
            click.echo(f"   • DPI: {dpi}")
            click.echo(f"   • Interactive: {interactive}")
            return

        # Create output directory
        output_path.mkdir(parents=True, exist_ok=True)

        # Initialize Hail
        hl.init(quiet=True)

        # Load MatrixTable
        click.echo(f"📥 Loading MatrixTable from {input}")
        mt = hl.read_matrix_table(input)

        # Check for QC annotations
        if 'sample_qc' not in mt.col and 'variant_qc' not in mt.row:
            click.echo("❌ No QC annotations found in MatrixTable. Run compute-qc first.", err=True)
            ctx.exit(1)

        # Extract QC metrics
        from hvantk.hgc.qc import extract_qc_metrics, QCMetrics

        # Create QCMetrics object
        sample_qc = mt.cols().select('sample_qc') if 'sample_qc' in mt.col else None
        variant_qc = mt.rows().select('variant_qc') if 'variant_qc' in mt.row else None

        qc_results = QCMetrics(
            mt=mt,
            sample_qc=sample_qc,
            variant_qc=variant_qc
        )

        click.echo(f"📊 Found QC data: Sample QC: {qc_results.has_sample_qc}, Variant QC: {qc_results.has_variant_qc}")

        # Check for interactive plotting
        if interactive:
            try:
                from hvantk.visualization.interactive_qc import check_plotly_available
                check_plotly_available()
                click.echo("🎨 Using interactive plotly plots")
                use_interactive = True
            except ImportError as e:
                click.echo(f"⚠️  Interactive plotting not available: {e}")
                click.echo("    Using standard matplotlib plots instead")
                use_interactive = False
        else:
            use_interactive = False

        plot_kwargs = {
            'style': style,
            'figsize': figsize_tuple,
            'dpi': dpi
        }

        # Generate plots based on type
        created_files = []

        if plot_type == 'overview' or plot_type == 'all':
            if qc_results.has_sample_qc:
                click.echo("🎨 Creating sample QC overview...")
                if use_interactive:
                    # For interactive, we'll create a dashboard instead of overview
                    fig = qc_results.plot_interactive_dashboard()
                    from hvantk.visualization.interactive_qc import save_interactive_plot
                    save_interactive_plot(
                        fig,
                        output_path / 'interactive_sample_overview.html',
                        format='html'
                    )
                    created_files.append('interactive_sample_overview.html')
                else:
                    fig = qc_results.plot_sample_overview(
                        save_path=output_path / f'sample_overview.{output_format}',
                        **plot_kwargs
                    )
                    created_files.append(f'sample_overview.{output_format}')

            if qc_results.has_variant_qc and not use_interactive:
                click.echo("🎨 Creating variant QC overview...")
                fig = qc_results.plot_variant_overview(
                    save_path=output_path / f'variant_overview.{output_format}',
                    **plot_kwargs
                )
                created_files.append(f'variant_overview.{output_format}')

        if plot_type == 'individual' or plot_type == 'all':
            if qc_results.has_sample_qc:
                click.echo("🎨 Creating individual sample plots...")

                if use_interactive:
                    from hvantk.visualization.interactive_qc import save_interactive_plot

                    # Sample call rates (interactive)
                    fig = qc_results.plot_interactive_sample_call_rates()
                    save_interactive_plot(fig, output_path / 'interactive_sample_call_rates.html')
                    created_files.append('interactive_sample_call_rates.html')

                    # Sample Ti/Tv (interactive)
                    try:
                        fig = qc_results.plot_interactive_sample_titv()
                        save_interactive_plot(fig, output_path / 'interactive_sample_titv.html')
                        created_files.append('interactive_sample_titv.html')
                    except (ValueError, KeyError, AttributeError, TypeError) as e:
                        click.echo(f"⚠️  Interactive Ti/Tv plot skipped: {type(e).__name__}: {e}")

                    # Sample scatter plot (bonus interactive feature)
                    try:
                        fig = qc_results.plot_interactive_sample_scatter()
                        save_interactive_plot(fig, output_path / 'interactive_sample_scatter.html')
                        created_files.append('interactive_sample_scatter.html')
                    except (ValueError, KeyError, AttributeError, TypeError) as e:
                        click.echo(f"⚠️  Interactive scatter plot skipped: {type(e).__name__}: {e}")
                else:
                    # Standard matplotlib plots
                    fig = qc_results.plot_sample_call_rates(
                        save_path=output_path / f'sample_call_rates.{output_format}',
                        **plot_kwargs
                    )
                    created_files.append(f'sample_call_rates.{output_format}')

                    # Sample Ti/Tv (if available)
                    try:
                        fig = qc_results.plot_sample_titv(
                            save_path=output_path / f'sample_titv.{output_format}',
                            **plot_kwargs
                        )
                        created_files.append(f'sample_titv.{output_format}')
                    except (ValueError, KeyError, AttributeError, TypeError) as e:
                        click.echo(f"⚠️  Ti/Tv plot skipped: {type(e).__name__}: {e}")

            if qc_results.has_variant_qc:
                click.echo("🎨 Creating individual variant plots...")

                if use_interactive:
                    from hvantk.visualization.interactive_qc import save_interactive_plot

                    # Variant call rates (interactive)
                    fig = qc_results.plot_interactive_variant_call_rates()
                    save_interactive_plot(fig, output_path / 'interactive_variant_call_rates.html')
                    created_files.append('interactive_variant_call_rates.html')

                    # Allele frequencies (interactive)
                    fig = qc_results.plot_interactive_allele_frequencies()
                    save_interactive_plot(fig, output_path / 'interactive_allele_frequencies.html')
                    created_files.append('interactive_allele_frequencies.html')

                    # HWE p-values (interactive)
                    fig = qc_results.plot_interactive_hwe_pvalues()
                    save_interactive_plot(fig, output_path / 'interactive_hwe_pvalues.html')
                    created_files.append('interactive_hwe_pvalues.html')
                else:
                    # Standard matplotlib plots
                    fig = qc_results.plot_variant_call_rates(
                        save_path=output_path / f'variant_call_rates.{output_format}',
                        **plot_kwargs
                    )
                    created_files.append(f'variant_call_rates.{output_format}')

                    # Allele frequencies
                    fig = qc_results.plot_allele_frequencies(
                        save_path=output_path / f'allele_frequencies.{output_format}',
                        **plot_kwargs
                    )
                    created_files.append(f'allele_frequencies.{output_format}')

                    # HWE p-values
                    fig = qc_results.plot_hwe_pvalues(
                        save_path=output_path / f'hwe_pvalues.{output_format}',
                        **plot_kwargs
                    )
                    created_files.append(f'hwe_pvalues.{output_format}')

        if plot_type == 'dashboard' or plot_type == 'all':
            click.echo("🎨 Creating comprehensive QC dashboard...")

            if use_interactive:
                from hvantk.visualization.interactive_qc import save_interactive_plot
                fig = qc_results.plot_interactive_dashboard()
                save_interactive_plot(fig, output_path / 'interactive_qc_dashboard.html')
                created_files.append('interactive_qc_dashboard.html')
            else:
                fig = qc_results.plot_dashboard(
                    save_path=output_path / f'qc_dashboard.{output_format}',
                    figsize=(20, 12),
                    **{k: v for k, v in plot_kwargs.items() if k != 'figsize'}
                )
                created_files.append(f'qc_dashboard.{output_format}')

        # Summary
        click.echo(f"\n📊 Successfully created {len(created_files)} QC plots:")
        for i, filename in enumerate(created_files, 1):
            file_path = output_path / filename
            file_size = file_path.stat().st_size / 1024  # KB
            click.echo(f"   {i:2d}. {filename} ({file_size:.1f} KB)")

        click.echo(f"\n📁 All plots saved in: {output_path.absolute()}")
        click.echo("✅ QC plotting completed successfully")

    except Exception as e:
        logger.exception(f"QC plotting failed: {e}")
        click.echo(f"❌ Error: {e}", err=True)
        ctx.exit(1)


@hgc_group.command(name="qc-report")
@click.option('--input', '-i', required=True, help='Input MatrixTable with QC annotations')
@click.option('--output', '-o', required=True, help='Output HTML report file')
@click.option('--title', default='Quality Control Report', help='Report title')
@click.option('--include-plots',
              multiple=True,
              type=click.Choice(['sample_overview', 'variant_overview', 'sample_call_rates',
                               'variant_call_rates', 'allele_frequencies', 'hwe', 'titv']),
              help='Specific plots to include (can be used multiple times)')
@click.option('--style',
              type=click.Choice(['default', 'publication']),
              default='default',
              help='Plot style for embedded plots')
@click.option('--dry-run', is_flag=True, help='Show what would be done without executing')
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
        hl.init(quiet=True)

        # Load MatrixTable
        click.echo(f"📥 Loading MatrixTable from {input}")
        mt = hl.read_matrix_table(input)

        # Check for QC annotations
        if 'sample_qc' not in mt.col and 'variant_qc' not in mt.row:
            click.echo("❌ No QC annotations found in MatrixTable. Run compute-qc first.", err=True)
            ctx.exit(1)

        # Extract QC metrics and create QCMetrics object
        from hvantk.hgc.qc import QCMetrics

        sample_qc = mt.cols().select('sample_qc') if 'sample_qc' in mt.col else None
        variant_qc = mt.rows().select('variant_qc') if 'variant_qc' in mt.row else None

        qc_results = QCMetrics(
            mt=mt,
            sample_qc=sample_qc,
            variant_qc=variant_qc
        )

        click.echo(f"📊 Found QC data: Sample QC: {qc_results.has_sample_qc}, Variant QC: {qc_results.has_variant_qc}")

        # Prepare report parameters
        report_kwargs = {
            'title': title
        }

        if include_plots:
            report_kwargs['include_plots'] = list(include_plots)

        # Generate HTML report
        click.echo("📝 Generating comprehensive HTML QC report...")

        report_path = qc_results.generate_html_report(
            output_path=output_path,
            **report_kwargs
        )

        # Get file size
        file_size = report_path.stat().st_size / 1024  # KB

        click.echo(f"📊 HTML QC report created:")
        click.echo(f"   • File: {report_path.absolute()}")
        click.echo(f"   • Size: {file_size:.1f} KB")
        click.echo(f"   • Title: {title}")

        sample_df = qc_results.get_sample_metrics_df() if qc_results.has_sample_qc else None
        variant_df = qc_results.get_variant_metrics_df() if qc_results.has_variant_qc else None

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


