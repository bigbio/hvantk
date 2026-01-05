#!/bin/bash
#
# Quick Start Example - Run HGC Scalability Benchmark
#
# This script provides a simple way to run the benchmark with default settings
# for the CHD_1000WGS chr20 data.

set -euo pipefail

# Configuration
# Set the GVCF directory path. You can set the GVCF_DIR environment variable before running this script,
# or edit the placeholder path below to your actual GVCF directory.
GVCF_DIR="${GVCF_DIR:-/path/to/your/gvcf_directory}"
if [ "$GVCF_DIR" = "/path/to/your/gvcf_directory" ]; then
    echo "WARNING: GVCF_DIR is set to a placeholder path. Please set the GVCF_DIR environment variable or edit the script to point to your actual GVCF directory."
fi
OUTPUT_DIR="./scalability_results_chr20_$(date +%Y%m%d)"
SAMPLE_SIZES="20,50,100,250,500,750,1000"
REFERENCE="GRCh38"

# Detect conda environment
if [ -n "${CONDA_DEFAULT_ENV:-}" ]; then
    CONDA_ENV="$CONDA_DEFAULT_ENV"
else
    CONDA_ENV="hvantk"
fi

echo "========================================================================"
echo "HGC Scalability Benchmark - Quick Start"
echo "========================================================================"
echo ""
echo "This will run the complete scalability benchmark with the following settings:"
echo "  - GVCF directory: $GVCF_DIR"
echo "  - Output directory: $OUTPUT_DIR"
echo "  - Sample sizes: $SAMPLE_SIZES"
echo "  - Reference genome: $REFERENCE"
echo "  - Conda environment: $CONDA_ENV"
echo ""
echo "Expected runtime: ~100+ hours for all 7 sample sizes"
echo ""
read -p "Do you want to continue? (y/N) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 0
fi

echo ""
echo "Starting benchmark..."
echo ""

# Run the benchmark
bash "$(dirname "$0")/hgc_scalability_benchmark.sh" \
    --gvcf-dir "$GVCF_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --sample-sizes "$SAMPLE_SIZES" \
    --reference "$REFERENCE" \
    --conda-env "$CONDA_ENV"

echo ""
echo "========================================================================"
echo "Benchmark completed!"
echo "========================================================================"
echo ""
echo "Results saved to: $OUTPUT_DIR"
echo ""
echo "To generate plots, run:"
echo "  python $(dirname "$0")/plot_scalability_results.py --results-dir $OUTPUT_DIR"
echo ""

