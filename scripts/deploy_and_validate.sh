#!/bin/bash
#
# Deployment and Validation Script for HVANTK Dataset Validation Framework
#
# This script sets up and runs the validation framework on a new host with conda/Hail/Spark
#

set -e  # Exit on any error

echo "=== HVANTK Dataset Validation Framework Deployment ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo

# Configuration
WORK_DIR="${WORK_DIR:-./dataset_validation}"
MAX_UCSC_DATASETS="${MAX_UCSC_DATASETS:-10}"
MAX_ATLAS_DATASETS="${MAX_ATLAS_DATASETS:-3}"
VERBOSE="${VERBOSE:-true}"

# Function to run validation with proper error handling
run_validation() {
    local cmd="$1"
    local description="$2"

    echo ">>> $description"
    echo "Command: $cmd"
    echo

    if eval "$cmd"; then
        echo "✅ SUCCESS: $description"
    else
        echo "❌ FAILED: $description"
        return 1
    fi
    echo
}

# Main deployment and validation workflow
main() {
    echo "Step 1: Environment Check"
    echo "------------------------"

    # Check if conda environment is active
    if [[ -z "$CONDA_DEFAULT_ENV" ]]; then
        echo "❌ ERROR: No conda environment active"
        echo "Please run: conda activate pyvatk"
        exit 1
    fi

    echo "✅ Conda environment: $CONDA_DEFAULT_ENV"

    # Check if hvantk is installed
    if ! python -c "import hvantk" 2>/dev/null; then
        echo "❌ ERROR: hvantk not installed"
        echo "Please run: pip install -e ."
        exit 1
    fi

    echo "✅ HVANTK package installed"
    echo

    echo "Step 2: List Available Datasets"
    echo "-------------------------------"

    run_validation \
        "python -m hvantk.commands.dataset_validation_cli list-datasets --source both" \
        "Listing available datasets"

    echo "Step 3: Validate Expression Atlas Datasets (Priority)"
    echo "---------------------------------------------------"

    # Test the known working Expression Atlas datasets first
    run_validation \
        "python -m hvantk.commands.dataset_validation_cli validate-datasets --expression-atlas-datasets E-MTAB-6798 E-MTAB-6814 E-MTAB-6769 --work-dir '$WORK_DIR' $([ '$VERBOSE' = 'true' ] && echo '--verbose')" \
        "Validating known working Expression Atlas datasets"

    echo "Step 4: Validate UCSC Datasets (Sample)"
    echo "---------------------------------------"

    # Test a few UCSC datasets
    run_validation \
        "python -m hvantk.commands.dataset_validation_cli validate-datasets --ucsc-datasets cortex-dev zeisel2015 h1-esc-diff --work-dir '$WORK_DIR' $([ '$VERBOSE' = 'true' ] && echo '--verbose')" \
        "Validating sample UCSC datasets"

    echo "Step 5: Batch Validation (Optional - Background)"
    echo "-----------------------------------------------"

    echo "To run full batch validation in background:"
    echo "nohup python -m hvantk.commands.dataset_validation_cli validate-datasets \\"
    echo "    --ucsc-datasets all \\"
    echo "    --expression-atlas-datasets all \\"
    echo "    --max-datasets $MAX_UCSC_DATASETS \\"
    echo "    --work-dir '$WORK_DIR' \\"
    echo "    --report-file validation_report.txt \\"
    echo "    --verbose > validation.log 2>&1 &"
    echo

    echo "Step 6: Generate Validation Report"
    echo "----------------------------------"

    run_validation \
        "python -m hvantk.commands.dataset_validation_cli report --output validation_summary.txt" \
        "Generating validation report"

    echo "Step 7: Show Status Summary"
    echo "--------------------------"

    run_validation \
        "python -m hvantk.commands.dataset_validation_cli status --show-successful" \
        "Showing validation status"

    echo "=== DEPLOYMENT COMPLETE ==="
    echo "✅ Validation framework is ready for use"
    echo "📊 Check validation_summary.txt for detailed results"
    echo "📁 Working directory: $WORK_DIR"
    echo

    # Show final summary
    echo "Summary of validated datasets:"
    if [[ -f validation_summary.txt ]]; then
        cat validation_summary.txt
    fi
}

# Run with error handling
if main "$@"; then
    echo "🎉 Deployment and validation completed successfully!"
    exit 0
else
    echo "💥 Deployment failed. Check the errors above."
    exit 1
fi
