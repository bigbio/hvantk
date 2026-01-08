#!/usr/bin/env python3
"""
HGC CPU Scaling Benchmark - Fixed Cohort, Variable CPU Count

This script runs the complete HGC workflow with a FIXED cohort size (e.g., 500 samples)
but varying the number of CPUs allocated to the Hail/Spark cluster.

Purpose: Measure strong scaling - how runtime changes with CPU count for fixed workload.

Usage:
    python hgc_cpu_scaling_benchmark.py \
        --gvcf-list samples_500.txt \
        --output-dir ./cpu_scaling_500 \
        --sample-size 500 \
        --num-cpus 32 \
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
    num_cpus: int,
    *,
    driver_memory_gb: int = 128,
    shuffle_partitions: int = 2048,
    max_partition_bytes_mb: int = 128,
    local_dir: Path | None = None,
    aqe_enabled: bool = False,
) -> None:
    """
    Initialize Hail/Spark on a *single node* for a CPU-scaling benchmark (local[N]).

    Benchmarking goal
    -----------------
    You want runtime to vary primarily as a function of CPU cores (N), *not* because Spark/Hail
    changes the amount of work. The biggest way benchmarks get "weird" is when partition counts
    or execution plans vary with N.

    Key benchmarking choices (and why)
    ----------------------------------
    1) master=local[N] is the CPU knob:
       - In local mode, Spark uses a thread pool of size N. This is the cleanest way to vary CPU.
       - Many spark.executor.* settings are not meaningful in local mode and can be misleading.

    2) Keep shuffle_partitions FIXED across runs:
       - If you scale shuffle partitions with num_cpus, you change the number of tasks/shuffle files,
         which changes overhead and can dominate runtime (especially for smaller workloads).
       - For 1000 WGS samples on chr20, 1024–4096 is usually reasonable; default here is 2048.

    3) Use a large, FIXED driver heap:
       - In local mode, the driver JVM is the main process. A stable heap avoids GC/spill swings
         across runs. Keep this constant for fair comparisons.

    4) AQE off for "pure scaling":
       - Adaptive Query Execution can change partitioning at runtime, muddying comparability.
       - You can run a separate experiment with AQE on for "best practical runtime".

    Practical benchmarking tips
    ---------------------------
    - Do one warm-up run (discard it) to stabilize JIT and filesystem cache.
    - Run 3 replicates per N and report median/min/max.
    - Keep tmp_dir and (optionally) spark.local.dir on fast local NVMe with lots of free space.

    Args:
        tmp_dir: Hail temporary directory (should be fast local storage for shuffle spill).
        log_file: Hail log file path.
        num_cpus: Number of local threads to use (e.g., 32, 64, 96, 128).
        driver_memory_gb: Fixed driver heap size across benchmark runs (default: 128).
        shuffle_partitions: Fixed shuffle partition count across benchmark runs (default: 2048).
        max_partition_bytes_mb: Input file split size for Spark SQL readers in MB (default: 128).
        local_dir: Optional Spark local directory (recommended: fast NVMe). If None, Spark default.
        aqe_enabled: Whether to enable Spark AQE (recommend False for pure CPU scaling, default: False).

    Returns:
        None (initializes Hail global context).
    """
    logger.info(f"Initializing Hail for CPU benchmark: local[{num_cpus}]")
    logger.info("Benchmark invariants (keep fixed across runs except num_cpus):")
    logger.info(f"  driver_memory_gb     = {driver_memory_gb}g")
    logger.info(f"  shuffle_partitions   = {shuffle_partitions}")
    logger.info(f"  max_partition_bytes  = {max_partition_bytes_mb} MB")
    logger.info(f"  AQE enabled          = {aqe_enabled}")
    if local_dir is not None:
        logger.info(f"  spark.local.dir      = {local_dir}")

    # Ensure directories exist
    tmp_dir.mkdir(parents=True, exist_ok=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    if local_dir is not None:
        local_dir.mkdir(parents=True, exist_ok=True)

    spark_conf = {
        # --- CPU scaling knob ---
        # local[N] controls the parallel thread pool; this is what you vary per benchmark point.
        # Keep everything else fixed for comparability.
        #
        # --- Memory (local mode) ---
        "spark.driver.memory": f"{driver_memory_gb}g",

        # --- Partitioning (fixed across CPU runs) ---
        "spark.sql.shuffle.partitions": str(shuffle_partitions),
        "spark.default.parallelism": str(shuffle_partitions),

        # Controls how Spark splits input files for reading (keep fixed across runs).
        "spark.sql.files.maxPartitionBytes": str(max_partition_bytes_mb * 1024 * 1024),

        # --- AQE (off for pure scaling, on for best-real runtime in a separate experiment) ---
        "spark.sql.adaptive.enabled": "true" if aqe_enabled else "false",
        "spark.sql.adaptive.coalescePartitions.enabled": "true" if aqe_enabled else "false",
    }

    # Optional: place shuffle spill & temp files on fast local disk
    if local_dir is not None:
        spark_conf["spark.local.dir"] = str(local_dir)

    # Initialize / reset Hail context
    try:
        # If re-running in the same Python process, stop the previous context cleanly.
        # (Hail raises if you init twice without stopping.)
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
            min_block_size=0,          # important for some combiner workflows
            default_reference="GRCh38",
            master=f"local[{num_cpus}]",
            spark_conf=spark_conf,
        )

        logger.info("✓ Hail initialized successfully")
        logger.info(f"  Hail version: {hl.__version__}")
        logger.info(f"  Spark master: local[{num_cpus}]")
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
    num_cpus: int,
    reference_genome: str = 'GRCh38'
) -> Dict[str, float]:
    """
    Run complete end-to-end HGC workflow: GVCF → VDS → MT → QC Tables + Cohort VCF

    Args:
        gvcf_files: List of GVCF file paths
        output_dir: Output directory for results
        sample_size: Number of samples (for logging)
        num_cpus: Number of CPUs allocated
        reference_genome: Reference genome (GRCh38 or GRCh37)

    Returns:
        Dictionary with timing breakdown for each step
    """
    timings = {}

    # Define paths - include CPU count in naming to avoid conflicts
    vds_path = str(output_dir / f"combined_{sample_size}_cpu{num_cpus}.vds")
    mt_path = str(output_dir / f"analysis_{sample_size}_cpu{num_cpus}.mt")
    sample_qc_path = str(output_dir / f"sample_qc_{sample_size}_cpu{num_cpus}.ht")
    variant_qc_path = str(output_dir / f"variant_qc_{sample_size}_cpu{num_cpus}.ht")
    vcf_path = str(output_dir / f"cohort_{sample_size}_cpu{num_cpus}.vcf.bgz")
    tmp_path = str(output_dir / "tmp")
    combiner_plan = str(output_dir / f"combiner_plan_{sample_size}_cpu{num_cpus}.json")

    logger.info("="*80)
    logger.info(f"Starting HGC workflow: {sample_size} samples @ {num_cpus} CPUs")
    logger.info("="*80)

    # STEP 1: Combine GVCFs → VDS
    logger.info(f"[{num_cpus} CPUs] Step 1/4: Combining {len(gvcf_files)} GVCFs to VDS...")
    logger.info(f"[{num_cpus} CPUs]   GVCF files:")
    for i, gvcf in enumerate(gvcf_files[:5], 1):  # Show first 5
        logger.info(f"[{num_cpus} CPUs]     {i}. {gvcf}")
    if len(gvcf_files) > 5:
        logger.info(f"[{num_cpus} CPUs]     ... and {len(gvcf_files) - 5} more")

    start = time.time()
    try:
        # Create a temporary directory with symlinks to the GVCF files
        gvcf_links_dir = output_dir / f"gvcf_links_cpu{num_cpus}"
        gvcf_links_dir.mkdir(exist_ok=True)

        logger.info(f"[{num_cpus} CPUs]   Creating symbolic links to GVCFs and their indexes...")
        created_links = 0
        for gvcf_path in gvcf_files:
            # Create symlink for GVCF file
            link_name = gvcf_links_dir / os.path.basename(gvcf_path)
            if os.path.lexists(link_name):
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
                    if not tbi_link_name.exists():
                        os.remove(tbi_link_name)
                        os.symlink(tbi_path, tbi_link_name)
                        created_links += 1
                else:
                    os.symlink(tbi_path, tbi_link_name)
                    created_links += 1
            else:
                logger.warning(f"[{num_cpus} CPUs]   Index file not found: {tbi_path}")

        logger.info(f"[{num_cpus} CPUs]   Created {created_links} symbolic links")

        logger.info(f"[{num_cpus} CPUs]   Calling hvantk combine_gvcfs...")
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
        logger.info(f"[{num_cpus} CPUs] Step 1 completed in {timings['gvcf_combine']:.1f}s")
        logger.info(f"[{num_cpus} CPUs]   → VDS written to: {vds_path}")
    except Exception as e:
        logger.error(f"[{num_cpus} CPUs] Step 1 FAILED: {e}")
        raise

    # STEP 2: Convert VDS → MatrixTable
    logger.info(f"[{num_cpus} CPUs] Step 2/4: Converting VDS to MatrixTable...")
    start = time.time()
    try:
        convert_vds_to_mt(
            vds_path=vds_path,
            output_path=mt_path,
            adjust_genotypes=True,
            skip_split_multi=False,
            skip_validation=False,
            overwrite=True
        )
        timings['vds_to_mt'] = time.time() - start
        logger.info(f"[{num_cpus} CPUs] Step 2 completed in {timings['vds_to_mt']:.1f}s")
        logger.info(f"[{num_cpus} CPUs]   → MatrixTable written to: {mt_path}")
    except Exception as e:
        logger.error(f"[{num_cpus} CPUs] Step 2 FAILED: {e}")
        raise

    # STEP 3: Compute QC metrics and export QC tables
    logger.info(f"[{num_cpus} CPUs] Step 3/4: Computing QC metrics and exporting QC tables...")
    start = time.time()
    try:
        mt = hl.read_matrix_table(mt_path)
        qc_metrics = compute_full_qc(mt)

        # Export QC tables separately
        qc_metrics.sample_qc.write(sample_qc_path, overwrite=True)
        qc_metrics.variant_qc.write(variant_qc_path, overwrite=True)

        timings['compute_qc'] = time.time() - start
        logger.info(f"[{num_cpus} CPUs] Step 3 completed in {timings['compute_qc']:.1f}s")
        logger.info(f"[{num_cpus} CPUs]   → Sample QC table: {sample_qc_path}")
        logger.info(f"[{num_cpus} CPUs]   → Variant QC table: {variant_qc_path}")
    except Exception as e:
        logger.error(f"[{num_cpus} CPUs] Step 3 FAILED: {e}")
        raise

    # STEP 4: Export clean MatrixTable to cohort VCF
    logger.info(f"[{num_cpus} CPUs] Step 4/4: Exporting clean cohort VCF...")
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
        logger.info(f"[{num_cpus} CPUs] Step 4 completed in {timings['mt_to_vcf']:.1f}s")
        logger.info(f"[{num_cpus} CPUs]   → Cohort VCF: {vcf_path}")
    except Exception as e:
        logger.error(f"[{num_cpus} CPUs] Step 4 FAILED: {e}")
        raise

    # Summary
    timings['total'] = sum(timings.values())
    timings['num_cpus'] = num_cpus
    timings['sample_size'] = sample_size

    logger.info("="*80)
    logger.info(f"[{num_cpus} CPUs] ✅ Complete workflow finished in {timings['total']:.1f}s")
    logger.info(f"[{num_cpus} CPUs] Timing breakdown:")
    logger.info(f"[{num_cpus} CPUs]   - GVCF → VDS:    {timings['gvcf_combine']:.1f}s ({timings['gvcf_combine']/timings['total']*100:.1f}%)")
    logger.info(f"[{num_cpus} CPUs]   - VDS → MT:      {timings['vds_to_mt']:.1f}s ({timings['vds_to_mt']/timings['total']*100:.1f}%)")
    logger.info(f"[{num_cpus} CPUs]   - Compute QC:    {timings['compute_qc']:.1f}s ({timings['compute_qc']/timings['total']*100:.1f}%)")
    logger.info(f"[{num_cpus} CPUs]   - MT → VCF:      {timings['mt_to_vcf']:.1f}s ({timings['mt_to_vcf']/timings['total']*100:.1f}%)")
    logger.info("="*80)

    return timings


def main():
    parser = argparse.ArgumentParser(
        description='Run HGC CPU scaling benchmark for fixed cohort size',
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
        '--num-cpus',
        type=int,
        required=True,
        help='Number of CPU cores to allocate to Hail/Spark'
    )
    parser.add_argument(
        '--driver-memory',
        type=int,
        default=128,
        help='Driver memory in GB (default: 128, kept FIXED across runs for fair comparison)'
    )
    parser.add_argument(
        '--shuffle-partitions',
        type=int,
        default=2048,
        help='Number of shuffle partitions (default: 2048, kept FIXED across runs)'
    )
    parser.add_argument(
        '--max-partition-bytes',
        type=int,
        default=128,
        help='Max partition bytes in MB (default: 128)'
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
        help='Enable Adaptive Query Execution (default: False for pure CPU scaling)'
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
    log_file = args.output_dir / f"workflow_cpu{args.num_cpus}.log"
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    ))
    logger.addHandler(file_handler)

    logger.info(f"HGC CPU Scaling Benchmark")
    logger.info(f"  Sample Size: {args.sample_size} (FIXED)")
    logger.info(f"  CPU Count: {args.num_cpus} (VARIABLE)")
    logger.info(f"  Driver Memory: {args.driver_memory}GB (FIXED)")
    logger.info(f"  Shuffle Partitions: {args.shuffle_partitions} (FIXED)")
    logger.info(f"  AQE Enabled: {args.aqe_enabled}")
    logger.info(f"  GVCF list: {args.gvcf_list}")
    logger.info(f"  Output directory: {args.output_dir}")
    logger.info(f"  Reference genome: {args.reference}")

    # Initialize Hail with specific CPU count
    tmp_dir = args.output_dir / "hail_tmp"
    tmp_dir.mkdir(exist_ok=True)
    hail_log = args.output_dir / f"hail_cpu{args.num_cpus}.log"
    setup_hail(
        tmp_dir=tmp_dir,
        log_file=hail_log,
        num_cpus=args.num_cpus,
        driver_memory_gb=args.driver_memory,
        shuffle_partitions=args.shuffle_partitions,
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
            num_cpus=args.num_cpus,
            reference_genome=args.reference
        )

        # Save timing results to JSON
        timing_file = args.output_dir / f"timing_cpu{args.num_cpus}.json"
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

