#!/usr/bin/env bash
# Run the ancestry pipeline, or prepare a 1000 Genomes reference, as a single-node job
# through the hpc-hail-job template.
#
#   ancestry-job.sh --query <query.mt> --reference <reference.mt> --out <dir>
#                   [--label-col super_pop] [--extra "<more hvantk ancestry-inference flags>"]
#                   [--cpus 16] [--mem-gb 64] [--time 0-01:00:00] [--partition <P>] [--name NAME] [--wait]
#   ancestry-job.sh --prepare --vcf <chr.vcf.gz> [--vcf ...] --samples-table <samples.ht> --out <reference.mt>
#                   [--build GRCh38] [--min-maf 0.05] [--target 20000] [--label-col super_pop] [... resources]
#
# Prints the job id. Give absolute paths: the job runs from the log directory.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SUBMIT="$HERE/../../hpc-hail-job/templates/submit-hail-job.sh"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"; HERE_REL="${HERE#"$R"/}"; T="$HVANTK_REPO/$HERE_REL"

PREPARE=0; QUERY=""; REFERENCE=""; OUT=""; LABEL="super_pop"; EXTRA=""
VCFS=(); SAMPLES_TABLE=""; BUILD="GRCh38"; MIN_MAF="0.05"; TARGET="20000"
CPUS=16; MEM_GB=64; TIME="0-01:00:00"; PARTITION=""; NAME=""; WAIT=0
while (( $# )); do
  case "$1" in
    --prepare) PREPARE=1; shift ;;
    --query) QUERY="$2"; shift 2 ;;
    --reference) REFERENCE="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --label-col) LABEL="$2"; shift 2 ;;
    --extra) EXTRA="$2"; shift 2 ;;
    --vcf) VCFS+=("$2"); shift 2 ;;
    --samples-table) SAMPLES_TABLE="$2"; shift 2 ;;
    --build) BUILD="$2"; shift 2 ;;
    --min-maf) MIN_MAF="$2"; shift 2 ;;
    --target) TARGET="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem-gb) MEM_GB="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    --wait) WAIT=1; shift ;;
    *) echo "ERROR: unknown option $1" >&2; exit 2 ;;
  esac
done
COMMON=(--cpus "$CPUS" --mem-gb "$MEM_GB" --time "$TIME")
[[ -n "$PARTITION" ]] && COMMON+=(--partition "$PARTITION")
(( WAIT )) && COMMON+=(--wait)
[[ -n "$OUT" && "$OUT" == /* ]] || { echo "ERROR: --out is required and must be absolute" >&2; exit 2; }

if (( PREPARE )); then
  (( ${#VCFS[@]} > 0 )) && [[ -n "$SAMPLES_TABLE" ]] || { echo "ERROR: --prepare needs --vcf and --samples-table" >&2; exit 2; }
  CMD=(python "$T/prepare-1kg-reference.py" --samples-table "$SAMPLES_TABLE" --out "$OUT" --build "$BUILD" --min-maf "$MIN_MAF" --target "$TARGET" --label-col "$LABEL")
  for v in "${VCFS[@]}"; do CMD+=(--vcf "$v"); done
  exec bash "$SUBMIT" "${COMMON[@]}" --name "${NAME:-ancestry-prepare}" -- "${CMD[@]}"
fi
[[ -n "$QUERY" && -n "$REFERENCE" ]] || { echo "ERROR: --query and --reference are required" >&2; exit 2; }
mkdir -p "$OUT"
CMD=(hvantk ancestry-inference -q "$QUERY" -r "$REFERENCE" --ancestry-col "$LABEL" -o "$OUT/predictions.ht" --output-dir "$OUT" --export-tsv --generate-report)
# shellcheck disable=SC2206
[[ -n "$EXTRA" ]] && CMD+=($EXTRA)
exec bash "$SUBMIT" "${COMMON[@]}" --name "${NAME:-ancestry}" -- "${CMD[@]}"
