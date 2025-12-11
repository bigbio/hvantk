#!/bin/bash
#
# HGC Scalability Benchmark - Main Orchestration Script
#
# This script runs the HGC workflow at different sample sizes to measure scalability.
# It samples GVCFs, runs the workflow, captures timing and memory metrics.
#
# Usage:
#   bash hgc_scalability_benchmark.sh [OPTIONS]
#
# Options:
#   --gvcf-dir DIR        Directory containing GVCF files (default: /mnt/nfs/KOL_UOL/projects/CHD_1000WGS/variant_calling/split_vcfs/chr20)
#   --output-dir DIR      Output directory for results (default: ./scalability_results)
#   --sample-sizes SIZES  Comma-separated sample sizes (default: 20,50,100,250,500,750,1000)
#   --reference REF       Reference genome (default: GRCh38)
#   --seed SEED          Random seed for sampling (default: 42)
#   --conda-env ENV      Conda environment name (default: auto-detect or use pyvatk)
#   --resume             Resume from previous run (skip completed runs)
#   --help               Show this help message

set -euo pipefail

# Default configuration
GVCF_DIR="/mnt/nfs/KOL_UOL/projects/CHD_1000WGS/variant_calling/split_vcfs/chr20"
OUTPUT_DIR="./scalability_results"
SAMPLE_SIZES="20,50,100,250,500,750,1000"
REFERENCE="GRCh38"
SEED=42
RESUME=false
CONDA_ENV=""  # Auto-detect or specify conda environment name

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --gvcf-dir)
            GVCF_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --sample-sizes)
            SAMPLE_SIZES="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --seed)
            SEED="$2"
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

# Convert comma-separated sizes to array
IFS=',' read -ra SIZES <<< "$SAMPLE_SIZES"

# Detect and activate conda environment
if [ -z "$CONDA_ENV" ]; then
    # Try to auto-detect conda environment
    # Use parameter expansion to avoid unbound variable error
    if [ -n "${CONDA_DEFAULT_ENV:-}" ]; then
        CONDA_ENV="$CONDA_DEFAULT_ENV"
        echo "Auto-detected conda environment: $CONDA_ENV"
    else
        # Default to pyvatk if no environment detected
        CONDA_ENV="pyvatk"
        echo "No conda environment detected, will try to use: $CONDA_ENV"
    fi
fi

# Initialize conda for bash if not already done
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
elif [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
    source "/opt/conda/etc/profile.d/conda.sh"
elif command -v conda &> /dev/null; then
    # Try to initialize conda if command exists
    eval "$(conda shell.bash hook 2>/dev/null)" || true
fi

# Activate the conda environment
if command -v conda &> /dev/null; then
    echo "Activating conda environment: $CONDA_ENV"
    conda activate "$CONDA_ENV" 2>/dev/null || {
        echo "WARNING: Failed to activate conda environment '$CONDA_ENV'"
        echo "Please activate it manually before running this script:"
        echo "  conda activate $CONDA_ENV"
        echo "Or specify a different environment with --conda-env"
        exit 1
    }
    echo "Conda environment activated: $(conda info --envs | grep '*' | awk '{print $1}')"
else
    echo "WARNING: conda not found. Make sure the correct Python environment is active."
    echo "Current Python: $(which python)"
fi

echo ""

echo "========================================================================"
echo "HGC Scalability Benchmark"
echo "========================================================================"
echo "GVCF directory:     $GVCF_DIR"
echo "Output directory:   $OUTPUT_DIR"
echo "Sample sizes:       ${SIZES[*]}"
echo "Reference genome:   $REFERENCE"
echo "Random seed:        $SEED"
echo "Conda environment:  $CONDA_ENV"
echo "Resume mode:        $RESUME"
echo "========================================================================"
echo ""

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="$SCRIPT_DIR/hgc_scalability_benchmark.py"

# Check if Python script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "ERROR: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

# Check if GVCF directory exists
if [ ! -d "$GVCF_DIR" ]; then
    echo "ERROR: GVCF directory not found: $GVCF_DIR"
    exit 1
fi

# Create output directory structure
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/sample_sets"
mkdir -p "$OUTPUT_DIR/run_logs"

# Create list of all GVCF files
echo "Scanning GVCF directory..."
GVCF_LIST="$OUTPUT_DIR/all_gvcfs.txt"
find "$GVCF_DIR" -name "*.g.vcf.gz" -o -name "*.gvcf.gz" -o -name "*.g.vcf" -o -name "*.gvcf" > "$GVCF_LIST"
TOTAL_GVCFS=$(wc -l < "$GVCF_LIST")
echo "Found $TOTAL_GVCFS GVCF files"
echo ""

if [ "$TOTAL_GVCFS" -eq 0 ]; then
    echo "ERROR: No GVCF files found in $GVCF_DIR"
    exit 1
fi

# Initialize CSV files
TIMING_CSV="$OUTPUT_DIR/timings.csv"
MEMORY_CSV="$OUTPUT_DIR/memory_usage.csv"

if [ ! -f "$TIMING_CSV" ] || [ "$RESUME" = false ]; then
    echo "sample_size,gvcf_combine_sec,vds_to_mt_sec,compute_qc_sec,mt_to_vcf_sec,total_sec" > "$TIMING_CSV"
fi

if [ ! -f "$MEMORY_CSV" ] || [ "$RESUME" = false ]; then
    echo "sample_size,peak_memory_mb,run_time_sec" > "$MEMORY_CSV"
fi

# Function to sample GVCFs
sample_gvcfs() {
    local n=$1
    local output_file=$2

    if [ "$n" -gt "$TOTAL_GVCFS" ]; then
        echo "WARNING: Requested $n samples but only $TOTAL_GVCFS available. Using all available GVCFs."
        n=$TOTAL_GVCFS
    fi

    # Use shuf with seed for reproducible sampling
    # On macOS, use gshuf if available, otherwise use sort -R with seed
    if command -v gshuf &> /dev/null; then
        gshuf --random-source=<(yes $SEED) -n "$n" "$GVCF_LIST" > "$output_file"
    elif command -v shuf &> /dev/null; then
        shuf --random-source=<(yes $SEED) -n "$n" "$GVCF_LIST" > "$output_file"
    else
        # Fallback for macOS without gshuf
        export RANDOM=$SEED
        sort -R "$GVCF_LIST" | head -n "$n" > "$output_file"
    fi

    echo "Sampled $n GVCFs to: $output_file"
}

# Function to parse memory from time output
parse_memory() {
    local log_file=$1
    local os_type=$(uname -s)

    if [[ "$os_type" == "Darwin" ]]; then
        # macOS: "maximum resident set size" in bytes
        local mem_bytes=$(grep "maximum resident set size" "$log_file" | awk '{print $1}')
        if [ -n "$mem_bytes" ]; then
            echo "scale=2; $mem_bytes / 1024 / 1024" | bc
        else
            echo "N/A"
        fi
    else
        # Linux: "Maximum resident set size (kbytes)"
        local mem_kb=$(grep "Maximum resident set size" "$log_file" | awk '{print $NF}')
        if [ -n "$mem_kb" ]; then
            echo "scale=2; $mem_kb / 1024" | bc
        else
            echo "N/A"
        fi
    fi
}

# Main benchmark loop
echo "Starting benchmark runs..."
echo ""

for SIZE in "${SIZES[@]}"; do
    echo "========================================================================"
    echo "Processing sample size: $SIZE"
    echo "========================================================================"

    # Check if already completed
    RUN_DIR="$OUTPUT_DIR/run_${SIZE}"
    TIMING_FILE="$RUN_DIR/timing_${SIZE}.json"

    if [ "$RESUME" = true ] && [ -f "$TIMING_FILE" ]; then
        echo "SKIPPING: Run for $SIZE samples already completed (found $TIMING_FILE)"
        echo ""
        continue
    fi

    # Create run directory
    mkdir -p "$RUN_DIR"

    # Sample GVCFs
    SAMPLE_LIST="$OUTPUT_DIR/sample_sets/samples_${SIZE}.txt"
    if [ ! -f "$SAMPLE_LIST" ] || [ "$RESUME" = false ]; then
        echo "Sampling $SIZE GVCFs..."
        sample_gvcfs "$SIZE" "$SAMPLE_LIST"
    else
        echo "Using existing sample list: $SAMPLE_LIST"
    fi

    # Run workflow with timing and memory tracking
    LOG_FILE="$OUTPUT_DIR/run_logs/run_${SIZE}.log"
    TIME_LOG="$OUTPUT_DIR/run_logs/time_${SIZE}.txt"

    echo "Running HGC workflow for $SIZE samples..."
    echo "Log file: $LOG_FILE"

    # Detect OS and use appropriate time command flags
    START_TIME=$(date +%s)
    OS_TYPE=$(uname -s)

    # Run the workflow with appropriate time command
    if [[ "$OS_TYPE" == "Darwin" ]]; then
        # macOS: Use /usr/bin/time -l for detailed stats
        echo "Detected macOS - using 'time -l'"
        /usr/bin/time -l python "$PYTHON_SCRIPT" \
            --gvcf-list "$SAMPLE_LIST" \
            --output-dir "$RUN_DIR" \
            --sample-size "$SIZE" \
            --reference "$REFERENCE" \
            > >(tee "$LOG_FILE") 2> >(tee "$TIME_LOG" >&2)
    else
        # Linux: Try /usr/bin/time -v for detailed stats
        echo "Detected Linux - using 'time -v'"
        # First check if -v flag is supported
        if /usr/bin/time -v echo test >/dev/null 2>&1; then
            /usr/bin/time -v python "$PYTHON_SCRIPT" \
                --gvcf-list "$SAMPLE_LIST" \
                --output-dir "$RUN_DIR" \
                --sample-size "$SIZE" \
                --reference "$REFERENCE" \
                > >(tee "$LOG_FILE") 2> >(tee "$TIME_LOG" >&2)
        else
            # Fallback: Use simple time without detailed stats
            echo "Warning: /usr/bin/time -v not supported, using basic timing"
            python "$PYTHON_SCRIPT" \
                --gvcf-list "$SAMPLE_LIST" \
                --output-dir "$RUN_DIR" \
                --sample-size "$SIZE" \
                --reference "$REFERENCE" \
                2>&1 | tee "$LOG_FILE"
            echo "No detailed memory stats available" > "$TIME_LOG"
        fi
    fi

    EXIT_CODE=${PIPESTATUS[0]}
    END_TIME=$(date +%s)
    WALL_TIME=$((END_TIME - START_TIME))

    if [ $EXIT_CODE -ne 0 ]; then
        echo "ERROR: Workflow failed for $SIZE samples (exit code: $EXIT_CODE)"
        echo "Check log file: $LOG_FILE"
        continue
    fi

    echo "Workflow completed successfully in ${WALL_TIME}s"

    # Extract timing from JSON
    if [ -f "$TIMING_FILE" ]; then
        GVCF_COMBINE=$(python -c "import json; print(json.load(open('$TIMING_FILE'))['gvcf_combine'])" 2>/dev/null || echo "N/A")
        VDS_TO_MT=$(python -c "import json; print(json.load(open('$TIMING_FILE'))['vds_to_mt'])" 2>/dev/null || echo "N/A")
        COMPUTE_QC=$(python -c "import json; print(json.load(open('$TIMING_FILE'))['compute_qc'])" 2>/dev/null || echo "N/A")
        MT_TO_VCF=$(python -c "import json; print(json.load(open('$TIMING_FILE'))['mt_to_vcf'])" 2>/dev/null || echo "N/A")
        TOTAL=$(python -c "import json; print(json.load(open('$TIMING_FILE'))['total'])" 2>/dev/null || echo "N/A")

        # Append to timing CSV
        echo "$SIZE,$GVCF_COMBINE,$VDS_TO_MT,$COMPUTE_QC,$MT_TO_VCF,$TOTAL" >> "$TIMING_CSV"
        echo "Timings saved to: $TIMING_CSV"
    else
        echo "WARNING: Timing file not found: $TIMING_FILE"
    fi

    # Extract memory usage
    if [ -f "$TIME_LOG" ]; then
        PEAK_MEMORY=$(parse_memory "$TIME_LOG")
        echo "$SIZE,$PEAK_MEMORY,$WALL_TIME" >> "$MEMORY_CSV"
        echo "Memory usage saved to: $MEMORY_CSV"
        echo "Peak memory: ${PEAK_MEMORY} MB"
    else
        echo "WARNING: Time log not found: $TIME_LOG"
    fi

    echo ""
done

echo "========================================================================"
echo "Benchmark completed!"
echo "========================================================================"
echo "Results saved to: $OUTPUT_DIR"
echo "  - Timings: $TIMING_CSV"
echo "  - Memory: $MEMORY_CSV"
echo ""
echo "To generate plots, run:"
echo "  python $SCRIPT_DIR/plot_scalability_results.py --results-dir $OUTPUT_DIR"
echo "========================================================================"

