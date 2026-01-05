#!/usr/bin/env python3
"""
HGC Scalability Benchmark - Python Workflow Runner

This script runs the complete HGC workflow for a given set of GVCF files:
1. Combine GVCFs → VDS
2. Convert VDS → MatrixTable
3. Compute QC metrics and export QC tables
4. Export clean cohort VCF

Usage:
    python hgc_scalability_benchmark.py \
        --gvcf-list samples_100.txt \
        --output-dir ./run_100 \
        --sample-size 100 \
        --reference GRCh38
"""

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Dict, List

import hail as hl

# Import HGC functions
from hvantk.hgc import (
    combine_gvcfs,
    convert_vds_to_mt,
    compute_full_qc,
    convert_mt_to_multi_sample_vcf
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


def setup_hail(tmp_dir: Path, log_file: Path):
    """Initialize Hail with appropriate settings."""
    logger.info("Initializing Hail...")

    # Ensure directories exist and have correct permissions
    tmp_dir.mkdir(parents=True, exist_ok=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    try:
        hl.init(
            tmp_dir=str(tmp_dir),
            log=str(log_file),
            quiet=False,
            append=False,
            min_block_size=0,  # Important for combiner
            default_reference='GRCh38'
        )
        logger.info(f"✓ Hail initialized successfully")
        logger.info(f"  Hail version: {hl.__version__}")
        logger.info(f"  Temp directory: {tmp_dir}")
        logger.info(f"  Log file: {log_file}")
    except Exception as e:
        logger.error(f"Failed to initialize Hail: {e}")
        raise


def read_gvcf_list(gvcf_list_file: Path) -> List[str]:
    """Read list of GVCF files from a text file."""
    with open(gvcf_list_file, 'r') as f:
        gvcf_files = [line.strip() for line in f if line.strip()]
    logger.info(f"Read {len(gvcf_files)} GVCF files from {gvcf_list_file}")
    return gvcf_files


def run_hgc_workflow(
    gvcf_files: List[str],
    output_dir: Path,
    sample_size: int,
    reference_genome: str = 'GRCh38'
) -> Dict[str, float]:
    """
    Run complete end-to-end HGC workflow: GVCF → VDS → MT → QC Tables + Cohort VCF

    Args:
        gvcf_files: List of GVCF file paths
        output_dir: Output directory for results
        sample_size: Number of samples (for logging)
        reference_genome: Reference genome (GRCh38 or GRCh37)

    Returns:
        Dictionary with timing breakdown for each step
    """
    timings = {}

    # Define paths
    vds_path = str(output_dir / f"combined_{sample_size}.vds")
    mt_path = str(output_dir / f"analysis_{sample_size}.mt")
    sample_qc_path = str(output_dir / f"sample_qc_{sample_size}.ht")
    variant_qc_path = str(output_dir / f"variant_qc_{sample_size}.ht")
    vcf_path = str(output_dir / f"cohort_{sample_size}.vcf.gz")
    tmp_path = str(output_dir / "tmp")
    combiner_plan = str(output_dir / f"combiner_plan_{sample_size}.json")

    logger.info("="*80)
    logger.info(f"Starting HGC workflow for {sample_size} samples")
    logger.info("="*80)

    # STEP 1: Combine GVCFs → VDS
    logger.info(f"[{sample_size}] Step 1/4: Combining {len(gvcf_files)} GVCFs to VDS...")
    logger.info(f"[{sample_size}]   GVCF files:")
    for i, gvcf in enumerate(gvcf_files[:5], 1):  # Show first 5
        logger.info(f"[{sample_size}]     {i}. {gvcf}")
    if len(gvcf_files) > 5:
        logger.info(f"[{sample_size}]     ... and {len(gvcf_files) - 5} more")

    start = time.time()
    try:
        # Create a temporary directory with symlinks to the GVCF files
        # This allows combine_gvcfs to find them via directory scanning
        gvcf_links_dir = output_dir / "gvcf_links"
        gvcf_links_dir.mkdir(exist_ok=True)

        logger.info(f"[{sample_size}]   Creating symbolic links to GVCFs and their indexes...")
        created_links = 0
        for gvcf_path in gvcf_files:
            # Create symlink for GVCF file
            link_name = gvcf_links_dir / os.path.basename(gvcf_path)
            if os.path.lexists(link_name):
                # Remove broken symlink if exists
                if not link_name.exists():
                    os.remove(link_name)
                    os.symlink(gvcf_path, link_name)
                    created_links += 1
            else:
                os.symlink(gvcf_path, link_name)
                created_links += 1

            # Create symlink for index file (.tbi)
            tbi_path = gvcf_path + '.tbi'
            tbi_link_name = gvcf_links_dir / (os.path.basename(gvcf_path) + '.tbi')
            if os.path.exists(tbi_path):
                if os.path.lexists(tbi_link_name):
                    # Remove broken symlink if exists
                    if not tbi_link_name.exists():
                        os.remove(tbi_link_name)
                        os.symlink(tbi_path, tbi_link_name)
                        created_links += 1
                else:
                    os.symlink(tbi_path, tbi_link_name)
                    created_links += 1
            else:
                logger.warning(f"[{sample_size}]   Index file not found: {tbi_path}")

        logger.info(f"[{sample_size}]   Created {created_links} symbolic links")

        # Verify symlinks are valid
        link_files = list(gvcf_links_dir.glob("*.g.vcf.gz"))
        logger.info(f"[{sample_size}]   Found {len(link_files)} GVCF files in symlink directory")
        if len(link_files) != len(gvcf_files):
            logger.warning(f"[{sample_size}]   Expected {len(gvcf_files)} but found {len(link_files)} GVCF files")

        logger.info(f"[{sample_size}]   Calling hvantk combine_gvcfs...")
        combine_gvcfs(
            gvcf_dir=str(gvcf_links_dir),
            vds_output_path=vds_path,
            tmp_path=tmp_path,
            save_path=combiner_plan,
            vdses=[],
            kwargs={},
            reference_genome=reference_genome
        )
        timings['gvcf_combine'] = time.time() - start
        logger.info(f"[{sample_size}] Step 1 completed in {timings['gvcf_combine']:.1f}s")
        logger.info(f"[{sample_size}]   → VDS written to: {vds_path}")
    except Exception as e:
        logger.error(f"[{sample_size}] Step 1 FAILED: {e}")
        raise

    # STEP 2: Convert VDS → MatrixTable
    logger.info(f"[{sample_size}] Step 2/4: Converting VDS to MatrixTable...")
    start = time.time()
    try:
        convert_vds_to_mt(
            vds_path=vds_path,
            output_path=mt_path,
            adjust_genotypes=True,
            skip_split_multi=False,
            convert_lgt_to_gt=True,
            overwrite=True
        )
        timings['vds_to_mt'] = time.time() - start
        logger.info(f"[{sample_size}] Step 2 completed in {timings['vds_to_mt']:.1f}s")
        logger.info(f"[{sample_size}]   → MatrixTable written to: {mt_path}")
    except Exception as e:
        logger.error(f"[{sample_size}] Step 2 FAILED: {e}")
        raise

    # STEP 3: Compute QC metrics and export QC tables
    logger.info(f"[{sample_size}] Step 3/4: Computing QC metrics and exporting QC tables...")
    start = time.time()
    try:
        mt = hl.read_matrix_table(mt_path)
        qc_metrics = compute_full_qc(mt)

        # Export QC tables separately (NOT in the VCF)
        qc_metrics.sample_qc.write(sample_qc_path, overwrite=True)
        qc_metrics.variant_qc.write(variant_qc_path, overwrite=True)

        timings['compute_qc'] = time.time() - start
        logger.info(f"[{sample_size}] Step 3 completed in {timings['compute_qc']:.1f}s")
        logger.info(f"[{sample_size}]   → Sample QC table: {sample_qc_path}")
        logger.info(f"[{sample_size}]   → Variant QC table: {variant_qc_path}")
    except Exception as e:
        logger.error(f"[{sample_size}] Step 3 FAILED: {e}")
        raise

    # STEP 4: Export clean MatrixTable to cohort VCF
    logger.info(f"[{sample_size}] Step 4/4: Exporting clean cohort VCF...")
    start = time.time()
    try:
        convert_mt_to_multi_sample_vcf(
            mt_path=mt_path,
            vcf_path=vcf_path,
            filter_adj_genotypes=True,
            min_ac=1,
            split_multi=True
        )
        timings['mt_to_vcf'] = time.time() - start
        logger.info(f"[{sample_size}] Step 4 completed in {timings['mt_to_vcf']:.1f}s")
        logger.info(f"[{sample_size}]   → Cohort VCF: {vcf_path}")
    except Exception as e:
        logger.error(f"[{sample_size}] Step 4 FAILED: {e}")
        raise

    # Summary
    timings['total'] = sum(timings.values())
    logger.info("="*80)
    logger.info(f"[{sample_size}] ✅ Complete workflow finished in {timings['total']:.1f}s")
    logger.info(f"[{sample_size}] Timing breakdown:")
    logger.info(f"[{sample_size}]   - GVCF → VDS:    {timings['gvcf_combine']:.1f}s ({timings['gvcf_combine']/timings['total']*100:.1f}%)")
    logger.info(f"[{sample_size}]   - VDS → MT:      {timings['vds_to_mt']:.1f}s ({timings['vds_to_mt']/timings['total']*100:.1f}%)")
    logger.info(f"[{sample_size}]   - Compute QC:    {timings['compute_qc']:.1f}s ({timings['compute_qc']/timings['total']*100:.1f}%)")
    logger.info(f"[{sample_size}]   - MT → VCF:      {timings['mt_to_vcf']:.1f}s ({timings['mt_to_vcf']/timings['total']*100:.1f}%)")
    logger.info(f"[{sample_size}] Final outputs:")
    logger.info(f"[{sample_size}]   - Cohort VCF: {vcf_path}")
    logger.info(f"[{sample_size}]   - Sample QC: {sample_qc_path}")
    logger.info(f"[{sample_size}]   - Variant QC: {variant_qc_path}")
    logger.info("="*80)

    return timings


def main():
    parser = argparse.ArgumentParser(
        description='Run HGC scalability benchmark for a given sample set',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--gvcf-list',
        type=Path,
        required=True,
        help='Text file with list of GVCF file paths (one per line)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        required=True,
        help='Output directory for results'
    )
    parser.add_argument(
        '--sample-size',
        type=int,
        required=True,
        help='Number of samples (for logging and naming)'
    )
    parser.add_argument(
        '--reference',
        type=str,
        default='GRCh38',
        choices=['GRCh38', 'GRCh37'],
        help='Reference genome (default: GRCh38)'
    )

    args = parser.parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Add file handler for logging
    log_file = args.output_dir / f"workflow_{args.sample_size}.log"
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    ))
    logger.addHandler(file_handler)

    logger.info(f"HGC Scalability Benchmark - Sample Size: {args.sample_size}")
    logger.info(f"GVCF list: {args.gvcf_list}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Reference genome: {args.reference}")

    # Initialize Hail
    tmp_dir = args.output_dir / "hail_tmp"
    tmp_dir.mkdir(exist_ok=True)
    hail_log = args.output_dir / f"hail_{args.sample_size}.log"
    setup_hail(tmp_dir, hail_log)

    # Read GVCF list
    gvcf_files = read_gvcf_list(args.gvcf_list)

    if len(gvcf_files) != args.sample_size:
        logger.warning(
            f"Expected {args.sample_size} samples but got {len(gvcf_files)} from list. "
            f"Proceeding with {len(gvcf_files)} samples."
        )

    # Run workflow
    try:
        timings = run_hgc_workflow(
            gvcf_files=gvcf_files,
            output_dir=args.output_dir,
            sample_size=args.sample_size,
            reference_genome=args.reference
        )

        # Save timing results to JSON
        timing_file = args.output_dir / f"timing_{args.sample_size}.json"
        with open(timing_file, 'w') as f:
            json.dump(timings, f, indent=2)
        logger.info(f"Timing results saved to: {timing_file}")

        # Return success
        logger.info("Workflow completed successfully!")
        sys.exit(0)

    except Exception as e:
        logger.error(f"Workflow failed with error: {e}")
        sys.exit(1)
    finally:
        # Stop Hail
        hl.stop()


if __name__ == '__main__':
    main()

