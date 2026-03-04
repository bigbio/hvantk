#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Build a Hail MatrixTable from the NYGC/CCDG 1000 Genomes (2020) per-chromosome
# joint-genotyped VCF files using hvantk.
#
# Usage:
#   # Interactive (foreground)
#   bash build_1kg_nygc.sh /path/to/vcf_dir /path/to/output.mt
#
#   # Background (survives SSH disconnect)
#   nohup bash build_1kg_nygc.sh /path/to/vcf_dir /path/to/output.mt &
#   # Logs go to: /path/to/output.mt.log
#
#   # Subset of chromosomes
#   bash build_1kg_nygc.sh /path/to/vcf_dir /path/to/output.mt chr1,chr2,chr22
#
# The VCF directory should contain the CCDG files:
#   20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr*.recalibrated_variants.vcf.gz
#   20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr*.recalibrated_variants.vcf.gz.tbi
#
# Since the directory also contains annotated VCFs (*.annotated.vcf.gz) and an
# "others" contig file, this script symlinks only the recalibrated genotype VCFs
# (chr1-chr22, chrX, chrY) into a staging directory before invoking hvantk.
# ---------------------------------------------------------------------------
set -euo pipefail

# ── Arguments ──────────────────────────────────────────────────────────────
VCF_DIR="${1:?Usage: $0 <vcf_dir> <output_mt> [chromosomes]}"
OUTPUT_MT="${2:?Usage: $0 <vcf_dir> <output_mt> [chromosomes]}"
CHROMOSOMES="${3:-}"          # optional: comma-separated, e.g. "chr1,chr2,chrX"

# ── Configuration ──────────────────────────────────────────────────────────
REFERENCE_GENOME="GRCh38"
OVERWRITE=false               # set to true to replace existing output
LOG_FILE="${OUTPUT_MT}.log"
SAMPLE_ANNOTATIONS=""         # path to sample annotations file (e.g. PED)
SAMPLE_ANNOTATIONS_DELIMITER=""  # delimiter for annotations file (e.g. " " for PED)

# ── Logging ────────────────────────────────────────────────────────────────
exec > >(tee -a "$LOG_FILE") 2>&1
echo "=========================================="
echo "  1000G NYGC MatrixTable Build"
echo "  Started: $(date)"
echo "  VCF dir: ${VCF_DIR}"
echo "  Output:  ${OUTPUT_MT}"
echo "  Log:     ${LOG_FILE}"
echo "=========================================="

# ── Stage VCFs ─────────────────────────────────────────────────────────────
# Create a staging directory with symlinks to only the recalibrated genotype
# VCFs (exclude *.annotated.* and *_others.*).
STAGE_DIR=$(mktemp -d "${TMPDIR:-/tmp}/1kg_stage_XXXXXX")
trap 'rm -rf "$STAGE_DIR"' EXIT

echo ""
echo "Staging recalibrated genotype VCFs into ${STAGE_DIR} ..."

count=0
for vcf in "${VCF_DIR}"/*_chr*.recalibrated_variants.vcf.gz; do
    [[ "$vcf" == *".annotated."* ]] && continue
    [[ "$vcf" == *"_others."* ]] && continue

    ln -s "$(realpath "$vcf")" "${STAGE_DIR}/$(basename "$vcf")"
    if [[ -f "${vcf}.tbi" ]]; then
        ln -s "$(realpath "${vcf}.tbi")" "${STAGE_DIR}/$(basename "${vcf}.tbi")"
    else
        echo "WARNING: Missing tabix index for $(basename "$vcf")" >&2
    fi
    count=$((count + 1))
done

echo "Staged ${count} VCF file(s)."

if [[ "$count" -eq 0 ]]; then
    echo "ERROR: No recalibrated genotype VCFs found in ${VCF_DIR}" >&2
    exit 1
fi

# ── Build MatrixTable ──────────────────────────────────────────────────────
ARGS=(
    --input-vcfs "$STAGE_DIR"
    --output-mt "$OUTPUT_MT"
    --reference-genome "$REFERENCE_GENOME"
)

[[ -n "$CHROMOSOMES" ]] && ARGS+=(--chromosomes "$CHROMOSOMES")
[[ "$OVERWRITE" == true ]] && ARGS+=(--overwrite)
[[ -n "$SAMPLE_ANNOTATIONS" ]] && ARGS+=(--sample-annotations "$SAMPLE_ANNOTATIONS")
[[ -n "$SAMPLE_ANNOTATIONS_DELIMITER" ]] && ARGS+=(--sample-annotations-delimiter "$SAMPLE_ANNOTATIONS_DELIMITER")

echo ""
echo "Running: hvantk build-1k-genome ${ARGS[*]}"
echo ""

hvantk build-1k-genome "${ARGS[@]}"

echo ""
echo "=========================================="
echo "  Completed: $(date)"
echo "  MatrixTable: ${OUTPUT_MT}"
echo "=========================================="
