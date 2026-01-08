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


def setup_hail(
    tmp_dir: Path,
    log_file: Path,
    *,
    driver_memory_gb: int = 128,
    max_partition_bytes_mb: int = 128,
    local_dir: Path | None = None,
    aqe_enabled: bool = False,
) -> None:
    """
    Initialize Hail for sample scalability benchmarking (weak scaling).

    Benchmarking goal
    -----------------
    You want runtime to scale linearly with sample count: T(n) ≈ a*n for large n.
    The amount of work grows with sample count, so partition count should also grow.

    Key benchmarking choices (and why)
    ----------------------------------
    1) Keep partition SIZE constant (max_partition_bytes), not COUNT:
       - Each partition should have ~128MB of data (Spark standard)
       - As samples increase, data increases, so partition count increases naturally
       - This reflects real-world scaling: more data = more partitions
       - Opposite of CPU scaling, where partition COUNT is fixed!

    2) Use a large, FIXED driver heap:
       - Genomic workflows need memory for sample/variant metadata
       - Fixed heap prevents GC variance across runs
       - Keep constant across all sample sizes

    3) AQE off for reproducibility:
       - Adaptive Query Execution changes execution plans at runtime
       - Makes runs non-comparable
       - For scalability measurement: OFF
       - For production: ON (for best performance)

    4) local[144] CPU allocation (fixed):
       - Constant CPU count isolates sample scaling
       - Measures algorithm complexity: O(n), O(n²), etc.
       - Not measuring CPU parallelism

    Practical benchmarking tips
    ---------------------------
    - Do one warm-up run (discard it) to stabilize JIT and filesystem cache.
    - Run each sample size once (not replicates) since you're measuring scaling, not variance.
    - Plot runtime vs sample size and fit to y = a*n + b to check linearity.
    - Expected R² > 0.98 for linear scaling

    Args:
        tmp_dir: Hail temporary directory (should be on fast local storage).
        log_file: Hail log file path.
        driver_memory_gb: Fixed driver heap size across all runs (default: 128GB).
        max_partition_bytes_mb: Max partition size in MB for input reading (default: 128MB).
                                Partition COUNT will scale with data size.
        local_dir: Optional Spark local directory for shuffle spill (recommend NVMe).
        aqe_enabled: Whether to enable Adaptive Query Execution (recommend False, default: False).

    Returns:
        None (initializes Hail global context).
    """
    logger.info("Initializing Hail for sample scalability benchmarking (weak scaling)")
    logger.info("Benchmark invariants (kept constant across sample size runs):")
    logger.info(f"  driver_memory_gb     = {driver_memory_gb}g")
    logger.info(f"  max_partition_bytes  = {max_partition_bytes_mb} MB")
    logger.info(f"  AQE enabled          = {aqe_enabled}")
    if local_dir is not None:
        logger.info(f"  spark.local.dir      = {local_dir}")
    logger.info("Note: Partition COUNT will scale with data size (samples)")
    logger.info("      This is CORRECT for weak scaling (unlike CPU scaling)")

    # Ensure directories exist
    tmp_dir.mkdir(parents=True, exist_ok=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    if local_dir is not None:
        local_dir.mkdir(parents=True, exist_ok=True)

    spark_conf = {
        # --- Memory (fixed across sample size runs) ---
        "spark.driver.memory": f"{driver_memory_gb}g",

        # --- Input file splitting (fixed to control work distribution) ---
        # Partition SIZE is fixed; partition COUNT scales with data (correct for weak scaling)
        "spark.sql.files.maxPartitionBytes": str(max_partition_bytes_mb * 1024 * 1024),

        # --- Shuffles (adaptive: let Spark decide based on data size) ---
        # Do not fix shuffle partitions! Let them scale with data.
        # Use Spark defaults or calculate from sample count.
        "spark.sql.shuffle.partitions": "200",  # Conservative default for small starts

        # --- AQE (off for reproducible weak scaling measurement) ---
        "spark.sql.adaptive.enabled": "true" if aqe_enabled else "false",
        "spark.sql.adaptive.coalescePartitions.enabled": "true" if aqe_enabled else "false",
    }

    # Optional: place shuffle spill & temp files on fast local disk
    if local_dir is not None:
        spark_conf["spark.local.dir"] = str(local_dir)

    # Initialize Hail context
    try:
        # Stop previous context if exists
        if hl.current_backend() is not None:
            try:
                hl.stop()
            except Exception:
                pass

        hl.init(
            tmp_dir=str(tmp_dir),
            log=str(log_file),
            quiet=False,
            append=False,
            min_block_size=0,          # Important for combiner workflows
            default_reference="GRCh38",
            master="local[144]",        # Fixed CPU allocation for weak scaling
            spark_conf=spark_conf,
        )

        logger.info("✓ Hail initialized successfully")
        logger.info(f"  Hail version: {hl.__version__}")
        logger.info(f"  Spark master: local[144]")
        logger.info(f"  tmp_dir: {tmp_dir}")
        logger.info(f"  log_file: {log_file}")

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
    vcf_path = str(output_dir / f"cohort_{sample_size}.vcf.bgz")
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
            skip_validation=False,  # Keep validation for benchmarking to catch issues
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
        '--driver-memory',
        type=int,
        default=128,
        help='Driver memory in GB (default: 128, kept FIXED across runs)'
    )
    parser.add_argument(
        '--max-partition-bytes',
        type=int,
        default=128,
        help='Max partition bytes in MB (default: 128, kept FIXED - partition count will scale with data)'
    )
    parser.add_argument(
        '--local-dir',
        type=Path,
        default=None,
        help='Spark local directory for shuffle spill (recommend fast NVMe)'
    )
    parser.add_argument(
        '--aqe-enabled',
        action='store_true',
        help='Enable Adaptive Query Execution (default: False for reproducible weak scaling)'
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

    logger.info(f"HGC Scalability Benchmark - Weak Scaling Test")
    logger.info(f"  Sample Size: {args.sample_size} (VARIABLE)")
    logger.info(f"  Driver Memory: {args.driver_memory}GB (FIXED)")
    logger.info(f"  Max Partition Bytes: {args.max_partition_bytes}MB (FIXED)")
    logger.info(f"  AQE Enabled: {args.aqe_enabled}")
    logger.info(f"  GVCF list: {args.gvcf_list}")
    logger.info(f"  Output directory: {args.output_dir}")
    logger.info(f"  Reference genome: {args.reference}")

    # Initialize Hail
    tmp_dir = args.output_dir / "hail_tmp"
    tmp_dir.mkdir(exist_ok=True)
    hail_log = args.output_dir / f"hail_{args.sample_size}.log"
    setup_hail(
        tmp_dir=tmp_dir,
        log_file=hail_log,
        driver_memory_gb=args.driver_memory,
        max_partition_bytes_mb=args.max_partition_bytes,
        local_dir=args.local_dir,
        aqe_enabled=args.aqe_enabled,
    )

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

