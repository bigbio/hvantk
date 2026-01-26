#!/bin/bash
#
# Generate GVCF Sample List for CPU Scaling Benchmark
#
# This script creates a sample list file with exactly N GVCF paths
# for use with the CPU scaling benchmark.
#
# Usage:
#   bash generate_sample_list.sh <gvcf_dir> <num_samples> <output_file>
#
# Example:
#   bash generate_sample_list.sh /data/gvcfs 500 samples_500.txt

set -euo pipefail

# Check arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <gvcf_dir> <num_samples> <output_file>"
    echo ""
    echo "Arguments:"
    echo "  gvcf_dir     Directory containing GVCF files (*.g.vcf.gz)"
    echo "  num_samples  Number of samples to select"
    echo "  output_file  Output file to write sample paths"
    echo ""
    echo "Example:"
    echo "  $0 /data/gvcfs 500 samples_500.txt"
    exit 1
fi

GVCF_DIR="$1"
NUM_SAMPLES="$2"
OUTPUT_FILE="$3"

# Validate GVCF directory
if [ ! -d "$GVCF_DIR" ]; then
    echo "ERROR: GVCF directory not found: $GVCF_DIR"
    exit 1
fi

# Count available GVCFs
echo "Searching for GVCF files in: $GVCF_DIR"
TOTAL_GVCFS=$(find "$GVCF_DIR" -name "*.g.vcf.gz" -type f | wc -l | tr -d ' ')
echo "Found $TOTAL_GVCFS GVCF files"

if [ "$TOTAL_GVCFS" -eq 0 ]; then
    echo "ERROR: No GVCF files found in $GVCF_DIR"
    exit 1
fi

if [ "$TOTAL_GVCFS" -lt "$NUM_SAMPLES" ]; then
    echo "ERROR: Not enough GVCF files. Requested $NUM_SAMPLES but only found $TOTAL_GVCFS"
    exit 1
fi

# Create sample list
echo "Selecting $NUM_SAMPLES samples randomly..."

# Use shuf if available (Linux), otherwise use sort -R (macOS/BSD)
if command -v shuf &> /dev/null; then
    # Linux: use shuf
    find "$GVCF_DIR" -name "*.g.vcf.gz" -type f | \
        shuf | \
        head -n "$NUM_SAMPLES" > "$OUTPUT_FILE"
elif command -v gshuf &> /dev/null; then
    # macOS with GNU coreutils installed (brew install coreutils)
    find "$GVCF_DIR" -name "*.g.vcf.gz" -type f | \
        gshuf | \
        head -n "$NUM_SAMPLES" > "$OUTPUT_FILE"
else
    # macOS/BSD: use sort -R (random sort)
    find "$GVCF_DIR" -name "*.g.vcf.gz" -type f | \
        sort -R | \
        head -n "$NUM_SAMPLES" > "$OUTPUT_FILE"
fi

# Verify
SELECTED=$(wc -l < "$OUTPUT_FILE" | tr -d ' ')

if [ "$SELECTED" -ne "$NUM_SAMPLES" ]; then
    echo "ERROR: Expected $NUM_SAMPLES but got $SELECTED in output"
    exit 1
fi

echo "✓ Successfully created sample list: $OUTPUT_FILE"
echo "  Selected: $SELECTED / $TOTAL_GVCFS GVCF files"
echo ""
echo "Sample list preview:"
head -n 5 "$OUTPUT_FILE"
if [ "$SELECTED" -gt 5 ]; then
    echo "  ... and $((SELECTED - 5)) more"
fi

echo ""
echo "To run CPU scaling benchmark with this list:"
echo "  bash hgc_cpu_scaling_benchmark.sh \\"
echo "    --gvcf-list $OUTPUT_FILE \\"
echo "    --output-dir ./cpu_scaling_$NUM_SAMPLES \\"
echo "    --sample-size $NUM_SAMPLES \\"
echo "    --cpu-counts 16,24,32,48,64,96,144"

