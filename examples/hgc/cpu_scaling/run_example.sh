#!/bin/bash
#
# Example: Run HGC CPU Scaling Benchmark
#
# This script demonstrates how to run a complete CPU scaling analysis
# for a fixed cohort of 500 samples.

set -euo pipefail

# =============================================================================
# Configuration - MODIFY THESE FOR YOUR SYSTEM
# =============================================================================

# Path to directory containing GVCF files
GVCF_DIR="/path/to/your/gvcf_directory"

# Fixed cohort size for CPU scaling test
SAMPLE_SIZE=500

# Output directory
OUTPUT_DIR="./cpu_scaling_results_500"

# CPU counts to test (adjust based on your available resources)
# Current system has 144 total CPUs, so testing from 16 to 144
CPU_COUNTS="16,24,32,48,64,96,144"

# Memory per CPU core (in GB)
MEMORY_PER_CORE=4

# Reference genome
REFERENCE="GRCh38"

# =============================================================================
# Execution
# =============================================================================

echo "============================================================================"
echo "HGC CPU Scaling Benchmark Example"
echo "============================================================================"
echo ""
echo "This will:"
echo "  1. Generate a fixed sample list of $SAMPLE_SIZE GVCFs"
echo "  2. Run the HGC workflow at different CPU counts: $CPU_COUNTS"
echo "  3. Generate scaling analysis plots and summary"
echo ""
echo "Configuration:"
echo "  GVCF Directory:    $GVCF_DIR"
echo "  Sample Size:       $SAMPLE_SIZE (FIXED)"
echo "  CPU Counts:        $CPU_COUNTS (VARIABLE)"
echo "  Memory per Core:   ${MEMORY_PER_CORE}GB"
echo "  Output Directory:  $OUTPUT_DIR"
echo "  Reference:         $REFERENCE"
echo "============================================================================"
echo ""

# Check if GVCF directory exists
if [ ! -d "$GVCF_DIR" ]; then
    echo "ERROR: GVCF directory not found: $GVCF_DIR"
    echo ""
    echo "Please edit this script and set GVCF_DIR to your actual GVCF directory."
    echo ""
    exit 1
fi

# =============================================================================
# Step 1: Generate Sample List
# =============================================================================

SAMPLE_LIST="${OUTPUT_DIR}/samples_${SAMPLE_SIZE}.txt"

echo "Step 1: Generating sample list..."
echo ""

mkdir -p "$OUTPUT_DIR"

if [ -f "$SAMPLE_LIST" ]; then
    echo "Sample list already exists: $SAMPLE_LIST"
    echo "Using existing list (delete to regenerate)"
else
    echo "Generating new sample list: $SAMPLE_LIST"
    bash "$(dirname "$0")/../common/generate_sample_list.sh" "$GVCF_DIR" "$SAMPLE_SIZE" "$SAMPLE_LIST"
fi

echo ""

# =============================================================================
# Step 2: Run CPU Scaling Benchmark
# =============================================================================

echo "Step 2: Running CPU scaling benchmark..."
echo ""
echo "This will run the complete HGC workflow at each CPU count:"
IFS=',' read -ra CPUS <<< "$CPU_COUNTS"
for cpu in "${CPUS[@]}"; do
    echo "  - $cpu CPUs"
done
echo ""
echo "Estimated total time: varies based on system (30min - 3hrs)"
echo ""

read -p "Press Enter to start benchmark (or Ctrl+C to cancel)..."

bash "$(dirname "$0")/benchmark.sh" \
    --gvcf-list "$SAMPLE_LIST" \
    --output-dir "$OUTPUT_DIR" \
    --sample-size "$SAMPLE_SIZE" \
    --cpu-counts "$CPU_COUNTS" \
    --memory-per-core "$MEMORY_PER_CORE" \
    --reference "$REFERENCE"

echo ""
echo "============================================================================"
echo "Benchmark Complete!"
echo "============================================================================"
echo ""

# =============================================================================
# Step 3: Display Results
# =============================================================================

echo "Results Summary:"
echo ""

if [ -f "$OUTPUT_DIR/cpu_scaling_summary.txt" ]; then
    cat "$OUTPUT_DIR/cpu_scaling_summary.txt"
else
    echo "Summary file not found. Check logs for errors."
fi

echo ""
echo "Generated Files:"
echo ""
echo "Summary:"
ls -lh "$OUTPUT_DIR"/*.{csv,txt} 2>/dev/null || echo "  (none)"
echo ""
echo "Plots:"
ls -lh "$OUTPUT_DIR"/*.png 2>/dev/null || echo "  (none)"
echo ""
echo "============================================================================"
echo ""
echo "Next Steps:"
echo ""
echo "1. View comprehensive plot:"
echo "   open $OUTPUT_DIR/cpu_scaling_comprehensive.png"
echo ""
echo "2. View individual plots:"
echo "   open $OUTPUT_DIR/cpu_scaling_speedup.png"
echo "   open $OUTPUT_DIR/cpu_scaling_efficiency.png"
echo ""
echo "3. Read detailed summary:"
echo "   cat $OUTPUT_DIR/cpu_scaling_summary.txt"
echo ""
echo "4. Analyze results to determine optimal CPU count for your workload"
echo ""
echo "============================================================================"

