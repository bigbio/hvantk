#!/usr/bin/env bash
# Joint-genotype a directory of gVCFs as a chain of four single-node jobs through the
# hpc-hail-job template, or describe a finished run. Login-node orchestrator.
#
#   pipeline-job.sh --gvcf-dir DIR --out DIR --prefix NAME [--n-partitions N] [--min-ac 1]
#                   [--cpus 8] [--mem-gb 32] [--time 0-00:30:00] [--partition <P>] [--name NAME] [--plan]
#   pipeline-job.sh --describe DIR --prefix NAME [--cpus 8] [--mem-gb 32] [--time 0-00:30:00]
#
# The chain: `hvantk hgc gvcf-combine` (gVCFs to a VDS), `hvantk hgc vds2mt` (dense
# MatrixTable), `hvantk hgc compute-qc` (sample and variant QC under <out>/qc), `hvantk
# hgc mt2vcf` (biallelic cohort VCF). Each job waits for the previous one and the chain
# stops at the first failure. Prints the job ids, one per line. --plan prints the chain
# without submitting. The combiner works on GRCh38 inputs; see the skill for GRCh37.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SUBMIT="$HERE/../../hpc-hail-job/templates/submit-hail-job.sh"
[[ -f "$SUBMIT" ]] || { echo "ERROR: hpc-hail-job's submit template not found next to this skill" >&2; exit 2; }
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"; HERE_REL="${HERE#"$R"/}"; T="$HVANTK_REPO/$HERE_REL"

GVCF_DIR=""; OUT=""; PREFIX="cohort"; NPART=""; MIN_AC=1; DESCRIBE=""
CPUS=8; MEM_GB=32; TIME="0-00:30:00"; PARTITION=""; NAME=""; PLAN=0
while (( $# )); do
  case "$1" in
    --gvcf-dir) GVCF_DIR="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --n-partitions) NPART="$2"; shift 2 ;;
    --min-ac) MIN_AC="$2"; shift 2 ;;
    --describe) DESCRIBE="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem-gb) MEM_GB="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    --plan) PLAN=1; shift ;;
    *) echo "ERROR: unknown option $1" >&2; exit 2 ;;
  esac
done
COMMON=(--cpus "$CPUS" --mem-gb "$MEM_GB" --time "$TIME" --wait)
[[ -n "$PARTITION" ]] && COMMON+=(--partition "$PARTITION")

if [[ -n "$DESCRIBE" ]]; then
  [[ "$DESCRIBE" == /* ]] || { echo "ERROR: --describe needs an absolute path" >&2; exit 2; }
  exec bash "$SUBMIT" --cpus "$CPUS" --mem-gb "$MEM_GB" --time "$TIME" --wait --name "${NAME:-hgc-describe}" -- \
    python "$T/describe-callset.py" "$DESCRIBE" --prefix "$PREFIX"
fi
[[ -n "$GVCF_DIR" && -n "$OUT" ]] || { echo "ERROR: --gvcf-dir and --out are required" >&2; exit 2; }
[[ "$GVCF_DIR" == /* && "$OUT" == /* ]] || { echo "ERROR: give absolute paths; jobs run from the log directory" >&2; exit 2; }
NAME="${NAME:-hgc-$PREFIX}"
NP=(); [[ -n "$NPART" ]] && NP=(--n-partitions "$NPART")
STEPS=(
  "hvantk hgc gvcf-combine -g $GVCF_DIR -o $OUT/$PREFIX.vds"
  "hvantk hgc vds2mt -i $OUT/$PREFIX.vds -o $OUT/$PREFIX.mt ${NP[*]:-}"
  "hvantk hgc compute-qc -i $OUT/$PREFIX.mt -o $OUT/qc --prefix $PREFIX"
  "hvantk hgc mt2vcf -i $OUT/$PREFIX.mt -o $OUT/$PREFIX.vcf.bgz --min-ac $MIN_AC"
)
if (( PLAN )); then printf '%s\n' "${STEPS[@]}"; exit 0; fi
for f in "$OUT/$PREFIX.vds" "$OUT/$PREFIX.mt" "$OUT/$PREFIX.vcf.bgz"; do
  [[ -e "$f" ]] && { echo "ERROR: refusing to overwrite $f" >&2; exit 1; }
done
mkdir -p "$OUT"
printf '%s\n' "${STEPS[@]}" > "$OUT/recipe.txt"
i=0
for step in "${STEPS[@]}"; do
  i=$((i + 1))
  # shellcheck disable=SC2206
  CMD=($step)
  bash "$SUBMIT" "${COMMON[@]}" --name "$NAME-$i" -- "${CMD[@]}"
done
