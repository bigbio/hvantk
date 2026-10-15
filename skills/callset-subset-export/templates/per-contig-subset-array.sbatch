#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=subset-contig
# One array task per contig: subset that contig's VDS by the sample list, densify,
# export. Submit from the login node with the site values on the command line:
#
#   source "$HVANTK_SITE_PROFILE" && sbatch --partition="$HVANTK_CPU_PARTITION" \
#     --array=1-24%6 --cpus-per-task=32 --mem=160G --time=0-02:00:00 \
#     --output="${HVANTK_LOG_DIR:-$HVANTK_WORK/hvantk-logs}/%x_%A_%a.out" \
#     --export=ALL,CONTIG_LIST=<file>,VDS_ROOT=<dir>,SAMPLES=<list>,OUT_ROOT=<dir>,N_PARTITIONS=50 \
#     "$HVANTK_REPO/skills/callset-subset-export/templates/per-contig-subset-array.sbatch"
#
# VDS_ROOT/<contig>.vds are the per-contig callsets. Outputs: OUT_ROOT/<contig>.vds,
# <contig>.mt, <contig>.vcf.bgz. The container environment comes from hpc-hail-job's
# template, which this script hands each command to.
set -euo pipefail
: "${HVANTK_SITE_PROFILE:?}" "${CONTIG_LIST:?}" "${VDS_ROOT:?}" "${SAMPLES:?}" "${OUT_ROOT:?}"
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
HAIL_JOB="$HVANTK_REPO/skills/hpc-hail-job/templates/hail-job.sbatch"
T="$HVANTK_REPO/skills/callset-subset-export/templates"
CONTIG="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CONTIG_LIST")"
[[ -n "$CONTIG" ]] || { echo "ERROR: no line $SLURM_ARRAY_TASK_ID in $CONTIG_LIST" >&2; exit 2; }
[[ -d "$VDS_ROOT/$CONTIG.vds" ]] || { echo "ERROR: no $CONTIG.vds under $VDS_ROOT" >&2; exit 2; }
mkdir -p "$OUT_ROOT"
NP=(); [[ -n "${N_PARTITIONS:-}" ]] && NP=(--n-partitions "$N_PARTITIONS")
echo "== contig $CONTIG (task $SLURM_ARRAY_TASK_ID) =="
bash "$HAIL_JOB" python "$T/subset-vds.py" --vds "$VDS_ROOT/$CONTIG.vds" --samples "$SAMPLES" --out "$OUT_ROOT/$CONTIG.vds"
bash "$HAIL_JOB" hvantk hgc vds2mt -i "$OUT_ROOT/$CONTIG.vds" -o "$OUT_ROOT/$CONTIG.mt" "${NP[@]}"
bash "$HAIL_JOB" hvantk hgc mt2vcf -i "$OUT_ROOT/$CONTIG.mt" -o "$OUT_ROOT/$CONTIG.vcf.bgz"
echo "== done $CONTIG =="
