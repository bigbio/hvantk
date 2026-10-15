#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=hgc-contig
# One array task per contig: combine that contig's gVCFs into a VDS, densify, export.
# Submit from the login node, with the site values on the command line, for example:
#
#   source "$HVANTK_SITE_PROFILE" && sbatch --partition="$HVANTK_CPU_PARTITION" \
#     --array=1-24%6 --cpus-per-task=32 --mem=160G --time=0-02:00:00 \
#     --output="${HVANTK_LOG_DIR:-$HVANTK_WORK/hvantk-logs}/%x_%A_%a.out" --export=ALL \
#     --export=ALL,CONTIG_LIST=<file with one contig per line>,GVCF_ROOT=<dir>,OUT_ROOT=<dir>,N_PARTITIONS=50 \
#     "$HVANTK_REPO/skills/joint-genotyping/templates/per-contig-array.sbatch"
#
# GVCF_ROOT/<contig>/ holds that contig's gVCFs (indexed, GRCh38 contig names). Outputs go to
# OUT_ROOT/<contig>.vds, <contig>.mt, <contig>.vcf.bgz. The container environment comes
# from hpc-hail-job's template, which this script hands each command to.
set -euo pipefail
: "${HVANTK_SITE_PROFILE:?}" "${CONTIG_LIST:?}" "${GVCF_ROOT:?}" "${OUT_ROOT:?}"
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
HAIL_JOB="$HVANTK_REPO/skills/hpc-hail-job/templates/hail-job.sbatch"
CONTIG="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CONTIG_LIST")"
[[ -n "$CONTIG" ]] || { echo "ERROR: no line $SLURM_ARRAY_TASK_ID in $CONTIG_LIST" >&2; exit 2; }
[[ -d "$GVCF_ROOT/$CONTIG" ]] || { echo "ERROR: no gVCF directory for $CONTIG under $GVCF_ROOT" >&2; exit 2; }
mkdir -p "$OUT_ROOT"
for f in "$OUT_ROOT/$CONTIG.vds" "$OUT_ROOT/$CONTIG.mt" "$OUT_ROOT/$CONTIG.vcf.bgz"; do
  [[ -e "$f" ]] && { echo "ERROR: refusing to overwrite $f" >&2; exit 1; }
done
echo "== contig $CONTIG (task $SLURM_ARRAY_TASK_ID) =="
bash "$HAIL_JOB" hvantk hgc gvcf-combine -g "$GVCF_ROOT/$CONTIG" -o "$OUT_ROOT/$CONTIG.vds"
NP=()
[[ -n "${N_PARTITIONS:-}" ]] && NP=(--n-partitions "$N_PARTITIONS")
bash "$HAIL_JOB" hvantk hgc vds2mt -i "$OUT_ROOT/$CONTIG.vds" -o "$OUT_ROOT/$CONTIG.mt" ${NP[@]+"${NP[@]}"}
bash "$HAIL_JOB" hvantk hgc compute-qc -i "$OUT_ROOT/$CONTIG.mt" -o "$OUT_ROOT/qc" --prefix "$CONTIG" --no-save-mt
bash "$HAIL_JOB" hvantk hgc mt2vcf -i "$OUT_ROOT/$CONTIG.mt" -o "$OUT_ROOT/$CONTIG.vcf.bgz"
echo "== done $CONTIG =="
