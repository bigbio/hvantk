#!/bin/bash
# Ancestry Inference CLI Examples
#
# This script demonstrates various CLI usage patterns for ancestry inference.
# Replace paths with your actual data locations.

set -e  # Exit on error

# ==============================================================================
# Configuration - Update these paths for your environment
# ==============================================================================

QUERY_MT="/path/to/your/cohort.mt"
REFERENCE_MT="/path/to/1kg_reference.mt"
OUTPUT_DIR="./ancestry_results"

# ==============================================================================
# Example 1: Basic Usage
# ==============================================================================
echo "=== Example 1: Basic Ancestry Inference ==="

hvantk ancestry-inference \
  -q "$QUERY_MT" \
  -r "$REFERENCE_MT" \
  --ancestry-col super_pop \
  -o "${OUTPUT_DIR}/basic/predictions.ht"

# ==============================================================================
# Example 2: With HTML Report and TSV Export
# ==============================================================================
echo ""
echo "=== Example 2: With Report and Export ==="

hvantk ancestry-inference \
  -q "$QUERY_MT" \
  -r "$REFERENCE_MT" \
  --ancestry-col super_pop \
  -o "${OUTPUT_DIR}/with_report/predictions.ht" \
  --output-dir "${OUTPUT_DIR}/with_report" \
  --generate-report \
  --export-tsv

# ==============================================================================
# Example 3: Conservative Assignment (Higher Confidence)
# ==============================================================================
echo ""
echo "=== Example 3: Conservative Assignment ==="

hvantk ancestry-inference \
  -q "$QUERY_MT" \
  -r "$REFERENCE_MT" \
  --ancestry-col super_pop \
  -o "${OUTPUT_DIR}/conservative/predictions.ht" \
  --min-prob 0.90 \
  --generate-report

# ==============================================================================
# Example 4: Custom Variant Filtering
# ==============================================================================
echo ""
echo "=== Example 4: Custom Variant Filtering ==="

hvantk ancestry-inference \
  -q "$QUERY_MT" \
  -r "$REFERENCE_MT" \
  --ancestry-col super_pop \
  -o "${OUTPUT_DIR}/custom_filter/predictions.ht" \
  --min-af 0.05 \
  --max-af 0.95 \
  --min-call-rate 0.99 \
  --apply-hwe-filter \
  --hwe-p 1e-6 \
  --generate-report

# ==============================================================================
# Example 5: Increased PCA Components
# ==============================================================================
echo ""
echo "=== Example 5: More PCA Components ==="

hvantk ancestry-inference \
  -q "$QUERY_MT" \
  -r "$REFERENCE_MT" \
  --ancestry-col super_pop \
  -o "${OUTPUT_DIR}/more_pcs/predictions.ht" \
  --n-pcs 30 \
  --n-pcs-classify 15 \
  --generate-report

# ==============================================================================
# Example 6: With Checkpointing (for large datasets)
# ==============================================================================
echo ""
echo "=== Example 6: With Checkpointing ==="

hvantk ancestry-inference \
  -q "$QUERY_MT" \
  -r "$REFERENCE_MT" \
  --ancestry-col super_pop \
  -o "${OUTPUT_DIR}/checkpointed/predictions.ht" \
  --checkpoint-path "${OUTPUT_DIR}/checkpointed/checkpoints" \
  --generate-report

# ==============================================================================
# Example 7: Skip Validation (faster, no CV metrics)
# ==============================================================================
echo ""
echo "=== Example 7: Skip Validation ==="

hvantk ancestry-inference \
  -q "$QUERY_MT" \
  -r "$REFERENCE_MT" \
  --ancestry-col super_pop \
  -o "${OUTPUT_DIR}/no_validation/predictions.ht" \
  --skip-validation \
  --generate-report

# ==============================================================================
# Example 8: Save Model and Loadings (for projection)
# ==============================================================================
echo ""
echo "=== Example 8: Save Model and Loadings ==="

hvantk ancestry-inference \
  -q "$QUERY_MT" \
  -r "$REFERENCE_MT" \
  --ancestry-col super_pop \
  -o "${OUTPUT_DIR}/with_model/predictions.ht" \
  --output-dir "${OUTPUT_DIR}/with_model" \
  --save-model \
  --save-loadings \
  --generate-report \
  --export-tsv

# ==============================================================================
# Example 9: Full Pipeline with All Options
# ==============================================================================
echo ""
echo "=== Example 9: Full Pipeline ==="

hvantk ancestry-inference \
  -q "$QUERY_MT" \
  -r "$REFERENCE_MT" \
  --ancestry-col super_pop \
  -o "${OUTPUT_DIR}/full_pipeline/predictions.ht" \
  --output-dir "${OUTPUT_DIR}/full_pipeline" \
  --min-af 0.01 \
  --max-af 0.99 \
  --min-call-rate 0.98 \
  --ld-r2 0.2 \
  --ld-window 500000 \
  --n-pcs 20 \
  --n-pcs-classify 10 \
  --n-estimators 100 \
  --min-prob 0.75 \
  --seed 42 \
  --n-cv-folds 5 \
  --checkpoint-path "${OUTPUT_DIR}/full_pipeline/checkpoints" \
  --generate-report \
  --export-tsv \
  --save-model \
  --save-loadings \
  --log-level INFO

# ==============================================================================
# Example 10: Debug Mode
# ==============================================================================
echo ""
echo "=== Example 10: Debug Mode ==="

hvantk ancestry-inference \
  -q "$QUERY_MT" \
  -r "$REFERENCE_MT" \
  --ancestry-col super_pop \
  -o "${OUTPUT_DIR}/debug/predictions.ht" \
  --log-level DEBUG

# ==============================================================================
# Summary
# ==============================================================================
echo ""
echo "=== All Examples Complete ==="
echo "Results saved to: $OUTPUT_DIR"
echo ""
echo "Output structure:"
echo "  ${OUTPUT_DIR}/"
echo "  ├── basic/predictions.ht"
echo "  ├── with_report/"
echo "  │   ├── predictions.ht"
echo "  │   ├── predictions.tsv"
echo "  │   └── ancestry_report.html"
echo "  ├── conservative/"
echo "  ├── custom_filter/"
echo "  ├── more_pcs/"
echo "  ├── checkpointed/"
echo "  ├── no_validation/"
echo "  ├── with_model/"
echo "  │   ├── predictions.ht"
echo "  │   ├── rf_model.pkl"
echo "  │   └── pca_loadings.ht"
echo "  └── full_pipeline/"
