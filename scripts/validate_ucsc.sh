#!/bin/bash
# UCSC Dataset Validation Convenience Script

set -e

echo "🔬 UCSC Dataset Automatic Validation"
echo "====================================="

# Default parameters
WORK_DIR="./dataset_validation"
SAMPLE_LINES=100
BATCH_SIZE=5
MAX_DATASETS=""
DRY_RUN=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN="--dry-run"
            shift
            ;;
        --test)
            MAX_DATASETS="--max-datasets 10"
            echo "🧪 Test mode: validating first 10 datasets only"
            shift
            ;;
        --small-batch)
            BATCH_SIZE=3
            echo "📦 Small batch mode: 3 datasets per batch"
            shift
            ;;
        --large-batch)
            BATCH_SIZE=10
            echo "📦 Large batch mode: 10 datasets per batch"
            shift
            ;;
        --quick)
            SAMPLE_LINES=50
            echo "⚡ Quick mode: 50 sample lines"
            shift
            ;;
        --thorough)
            SAMPLE_LINES=200
            echo "🔍 Thorough mode: 200 sample lines"
            shift
            ;;
        --work-dir)
            WORK_DIR="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --dry-run        Show what would be validated without running"
            echo "  --test           Validate only first 10 datasets (for testing)"
            echo "  --small-batch    Use batch size of 3 (safer for limited resources)"
            echo "  --large-batch    Use batch size of 10 (faster with good resources)"
            echo "  --quick          Use 50 sample lines (faster validation)"
            echo "  --thorough       Use 200 sample lines (more thorough validation)"
            echo "  --work-dir DIR   Specify working directory"
            echo "  --help           Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0 --dry-run                    # See what would be validated"
            echo "  $0 --test                       # Test with first 10 datasets"
            echo "  $0 --quick --small-batch        # Fast validation with small batches"
            echo "  $0 --thorough --large-batch     # Thorough validation with large batches"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Show configuration
echo ""
echo "Configuration:"
echo "  Work Directory: $WORK_DIR"
echo "  Sample Lines: $SAMPLE_LINES"
echo "  Batch Size: $BATCH_SIZE"
if [[ -n "$MAX_DATASETS" ]]; then
    echo "  Max Datasets: 10 (test mode)"
fi
if [[ -n "$DRY_RUN" ]]; then
    echo "  Mode: DRY RUN"
fi
echo ""

# Run the validation script
python scripts/validate_all_ucsc_datasets.py \
    --work-dir "$WORK_DIR" \
    --sample-lines $SAMPLE_LINES \
    --batch-size $BATCH_SIZE \
    $MAX_DATASETS \
    $DRY_RUN \
    --continue-on-error

echo ""
echo "🎉 Validation script completed!"
