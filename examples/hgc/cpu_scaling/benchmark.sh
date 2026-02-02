#!/bin/bash
#
# HGC CPU Scaling Benchmark - Orchestration Script
#
# Runs the HGC workflow with a FIXED cohort size but VARYING CPU counts
# to measure strong scaling (speedup vs. CPU count)
#
# Key benchmarking principle: Keep all parameters FIXED except CPU count
# to ensure runtime varies primarily as a function of CPU cores.
#
# Usage:
#   bash hgc_cpu_scaling_benchmark.sh [OPTIONS]
#
# Options:
#   --gvcf-list FILE           File containing list of GVCF paths (default: samples_500.txt)
#   --output-dir DIR           Output directory for results (default: ./cpu_scaling_results)
#   --sample-size N            Fixed sample size (default: 500)
#   --cpu-counts COUNTS        Comma-separated CPU counts (default: 16,24,32,48,64,96,144)
#   --driver-memory GB         Driver memory in GB, FIXED across runs (default: 128)
#   --shuffle-partitions N     Shuffle partitions, FIXED across runs (default: 2048)
#   --max-partition-bytes MB   Max partition bytes in MB (default: 128)
#   --local-dir DIR            Spark local dir for spill (recommend fast NVMe)
#   --aqe-enabled              Enable Adaptive Query Execution (default: off for pure scaling)
#   --reference REF            Reference genome (default: GRCh38)
#   --conda-env ENV            Conda environment name (default: auto-detect or use hvantk)
#   --resume                   Resume from previous run (skip completed runs)
#   --help                     Show this help message

set -euo pipefail

# Default configuration
GVCF_LIST="samples_500.txt"
OUTPUT_DIR="./cpu_scaling_results"
SAMPLE_SIZE=500
CPU_COUNTS="16,24,32,48,64,96,144"
DRIVER_MEMORY=128
SHUFFLE_PARTITIONS=2048
MAX_PARTITION_BYTES=128
LOCAL_DIR=""
AQE_ENABLED=false
REFERENCE="GRCh38"
RESUME=false
CONDA_ENV=""

# =============================================================================
# Proxy Configuration for Hail (Critical for HPC/Corporate Environments)
# =============================================================================
echo "Configuring proxy settings for Hail backend..."

# Save original proxy settings
ORIGINAL_HTTP_PROXY="${HTTP_PROXY:-}"
ORIGINAL_HTTPS_PROXY="${HTTPS_PROXY:-}"

# Bypass proxy for localhost traffic (Hail/Spark communication)
export NO_PROXY="localhost,127.0.0.1,0.0.0.0,::1"
export no_proxy="localhost,127.0.0.1,0.0.0.0,::1"
unset HTTP_PROXY HTTPS_PROXY http_proxy https_proxy

echo "  ✓ Proxy disabled for localhost (Hail backend communication)"
echo "  NO_PROXY: $NO_PROXY"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --gvcf-list)
            GVCF_LIST="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --sample-size)
            SAMPLE_SIZE="$2"
            shift 2
            ;;
        --cpu-counts)
            CPU_COUNTS="$2"
            shift 2
            ;;
        --driver-memory)
            DRIVER_MEMORY="$2"
            shift 2
            ;;
        --shuffle-partitions)
            SHUFFLE_PARTITIONS="$2"
            shift 2
            ;;
        --max-partition-bytes)
            MAX_PARTITION_BYTES="$2"
            shift 2
            ;;
        --local-dir)
            LOCAL_DIR="$2"
            shift 2
            ;;
        --aqe-enabled)
            AQE_ENABLED=true
            shift
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        --resume)
            RESUME=true
            shift
            ;;
        --help)
            grep "^#" "$0" | sed 's/^# *//'
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Run with --help for usage information"
            exit 1
            ;;
    esac
done

# Convert comma-separated CPU counts to array
IFS=',' read -ra CPUS <<< "$CPU_COUNTS"

# Detect and activate conda environment
if [ -z "$CONDA_ENV" ]; then
    if [ -n "${CONDA_DEFAULT_ENV:-}" ]; then
        CONDA_ENV="$CONDA_DEFAULT_ENV"
        echo "Using current conda environment: $CONDA_ENV"
    else
        CONDA_ENV="hvantk"
        echo "No conda environment specified, using default: $CONDA_ENV"
    fi
fi

# Activate conda environment
echo "Activating conda environment: $CONDA_ENV"
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV" || {
    echo "ERROR: Failed to activate conda environment '$CONDA_ENV'"
    echo "Please create the environment or specify a different one with --conda-env"
    exit 1
}

# Verify Python environment
echo ""
echo "Environment verification:"
python --version
echo "Python location: $(which python)"
echo "Hail installation:"
python -c "import hail; print(f'  Hail version: {hail.__version__}')" || {
    echo "ERROR: Hail not found in environment"
    exit 1
}
echo ""

# Get script directory for finding benchmark.py and plot_results.py
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Verify GVCF list exists
if [ ! -f "$GVCF_LIST" ]; then
    echo "ERROR: GVCF list file not found: $GVCF_LIST"
    exit 1
fi

# Count GVCFs in list
GVCF_COUNT=$(wc -l < "$GVCF_LIST" | tr -d ' ')
echo "Found $GVCF_COUNT GVCF files in $GVCF_LIST"

if [ "$GVCF_COUNT" -ne "$SAMPLE_SIZE" ]; then
    echo "WARNING: GVCF count ($GVCF_COUNT) != requested sample size ($SAMPLE_SIZE)"
    echo "Proceeding with $GVCF_COUNT samples..."
fi

# Print benchmark configuration
echo ""
echo "============================================================================"
echo "HGC CPU Scaling Benchmark Configuration"
echo "============================================================================"
echo "Fixed Parameters (kept constant across runs):"
echo "  Cohort Size:         $SAMPLE_SIZE samples"
echo "  Driver Memory:       ${DRIVER_MEMORY}GB"
echo "  Shuffle Partitions:  $SHUFFLE_PARTITIONS"
echo "  Max Partition Bytes: ${MAX_PARTITION_BYTES}MB"
echo "  AQE Enabled:         $AQE_ENABLED"
if [ -n "$LOCAL_DIR" ]; then
    echo "  Spark Local Dir:     $LOCAL_DIR"
fi
echo ""
echo "Variable Parameter:"
echo "  CPU Counts:          ${CPUS[*]}"
echo ""
echo "Other Settings:"
echo "  GVCF List:           $GVCF_LIST"
echo "  Output Directory:    $OUTPUT_DIR"
echo "  Reference Genome:    $REFERENCE"
echo "Resume Mode:           $RESUME"
echo "============================================================================"
echo ""

# Create summary CSV
SUMMARY_CSV="$OUTPUT_DIR/cpu_scaling_summary.csv"
if [ ! -f "$SUMMARY_CSV" ]; then
    echo "num_cpus,sample_size,gvcf_combine_sec,vds_to_mt_sec,compute_qc_sec,mt_to_vcf_sec,total_sec,speedup,efficiency" > "$SUMMARY_CSV"
fi

# Baseline timing (for speedup calculation)
BASELINE_TIME=""

# Run benchmark for each CPU count
TOTAL_RUNS=${#CPUS[@]}
CURRENT_RUN=0

for NUM_CPUS in "${CPUS[@]}"; do
    CURRENT_RUN=$((CURRENT_RUN + 1))

    echo ""
    echo "========================================================================"
    echo "Run $CURRENT_RUN/$TOTAL_RUNS: Testing with $NUM_CPUS CPUs"
    echo "========================================================================"

    # Check if already completed
    TIMING_FILE="$OUTPUT_DIR/timing_cpu${NUM_CPUS}.json"
    if [ "$RESUME" = true ] && [ -f "$TIMING_FILE" ]; then
        echo "✓ Already completed (found $TIMING_FILE). Skipping..."

        # Read timing for baseline calculation
        if [ -z "$BASELINE_TIME" ]; then
            BASELINE_TIME=$(python3 -c "import json; print(json.load(open('$TIMING_FILE'))['total'])")
            echo "  Using as baseline time: ${BASELINE_TIME}s"
        fi

        continue
    fi

    # Run benchmark
    echo "Starting workflow with $NUM_CPUS CPUs..."
    START_TIME=$(date +%s)

    # Build command with required arguments
    CMD="python \"$SCRIPT_DIR/benchmark.py\" \
        --gvcf-list \"$GVCF_LIST\" \
        --output-dir \"$OUTPUT_DIR\" \
        --sample-size $SAMPLE_SIZE \
        --num-cpus $NUM_CPUS \
        --driver-memory $DRIVER_MEMORY \
        --shuffle-partitions $SHUFFLE_PARTITIONS \
        --max-partition-bytes $MAX_PARTITION_BYTES \
        --reference $REFERENCE"

    # Add optional arguments
    if [ -n "$LOCAL_DIR" ]; then
        CMD="$CMD --local-dir \"$LOCAL_DIR\""
    fi

    if [ "$AQE_ENABLED" = true ]; then
        CMD="$CMD --aqe-enabled"
    fi

    # Execute
    eval $CMD || {
        echo "ERROR: Workflow failed for $NUM_CPUS CPUs"
        echo "Check logs in: $OUTPUT_DIR/workflow_cpu${NUM_CPUS}.log"
        exit 1
    }

    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))

    echo "✓ Completed in ${ELAPSED}s"

    # Extract timings from JSON
    if [ -f "$TIMING_FILE" ]; then
        TIMINGS=$(python3 -c "
import json
with open('$TIMING_FILE') as f:
    t = json.load(f)
print(f\"{t['gvcf_combine']:.2f},{t['vds_to_mt']:.2f},{t['compute_qc']:.2f},{t['mt_to_vcf']:.2f},{t['total']:.2f}\")
")

        TOTAL_TIME=$(echo "$TIMINGS" | cut -d',' -f5)

        # Set baseline if this is the first run
        if [ -z "$BASELINE_TIME" ]; then
            BASELINE_TIME="$TOTAL_TIME"
            echo "  Set baseline time: ${BASELINE_TIME}s (${NUM_CPUS} CPUs)"
        fi

        # Calculate speedup and efficiency
        SPEEDUP=$(python3 -c "print(f'{$BASELINE_TIME / $TOTAL_TIME:.3f}')")
        EFFICIENCY=$(python3 -c "print(f'{($BASELINE_TIME / $TOTAL_TIME) / ($NUM_CPUS / ${CPUS[0]}) * 100:.1f}')")

        # Append to summary CSV
        echo "$NUM_CPUS,$SAMPLE_SIZE,$TIMINGS,$SPEEDUP,$EFFICIENCY" >> "$SUMMARY_CSV"

        echo "  Speedup vs baseline: ${SPEEDUP}x"
        echo "  Parallel efficiency: ${EFFICIENCY}%"
    else
        echo "WARNING: Timing file not found: $TIMING_FILE"
    fi
done

echo ""
echo "============================================================================"
echo "CPU Scaling Benchmark Complete!"
echo "============================================================================"
echo "Results saved to: $OUTPUT_DIR"
echo "Summary CSV: $SUMMARY_CSV"
echo ""

# Generate plots
echo "Generating plots..."
if python "$SCRIPT_DIR/plot_results.py" --results-dir "$OUTPUT_DIR"; then
    echo "✓ Plots generated successfully"
else
    echo "WARNING: Plot generation failed (check if plot_cpu_scaling_results.py exists)"
fi

echo ""
echo "To view results:"
echo "  cat $SUMMARY_CSV"
echo "  ls $OUTPUT_DIR/*.png"
echo ""

