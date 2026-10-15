#!/usr/bin/env bash
# Run the HGC pipeline, or describe its outputs, as a single-node job through the
# hpc-hail-job template.
#
#   pipeline-job.sh --gvcf-dir DIR --out DIR --prefix NAME [--build GRCh38|GRCh37]
#                   [--n-partitions N] [--extra "<more hvantk hgc pipeline flags>"]
#                   [--cpus 8] [--mem-gb 32] [--time 0-00:30:00] [--partition <P>]
#                   [--name NAME] [--dry-run-pipeline] [--wait]
#   pipeline-job.sh --describe DIR --prefix NAME [--cpus 2] [--mem-gb 8] [--wait]
#
# --dry-run-pipeline submits a job that only prints the pipeline's plan. --describe
# submits a job that prints the sample set, the MatrixTable row count and the VCF record
# count of a finished run (and writes <out>/counts.json). Prints the job id.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SUBMIT="$HERE/../../hpc-hail-job/templates/submit-hail-job.sh"
[[ -f "$SUBMIT" ]] || { echo "ERROR: hpc-hail-job's submit template not found next to this skill" >&2; exit 2; }

GVCF_DIR=""; OUT=""; PREFIX="cohort"; BUILD="GRCh38"; NPART=""; EXTRA=""; DESCRIBE=""
CPUS=8; MEM_GB=32; TIME="0-00:30:00"; PARTITION=""; NAME=""; PLAN=0; WAIT=0
while (( $# )); do
  case "$1" in
    --gvcf-dir) GVCF_DIR="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --build) BUILD="$2"; shift 2 ;;
    --n-partitions) NPART="$2"; shift 2 ;;
    --extra) EXTRA="$2"; shift 2 ;;
    --describe) DESCRIBE="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem-gb) MEM_GB="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    --dry-run-pipeline) PLAN=1; shift ;;
    --wait) WAIT=1; shift ;;
    *) echo "ERROR: unknown option $1" >&2; exit 2 ;;
  esac
done
COMMON=(--cpus "$CPUS" --mem-gb "$MEM_GB" --time "$TIME")
[[ -n "$PARTITION" ]] && COMMON+=(--partition "$PARTITION")
(( WAIT )) && COMMON+=(--wait)

if [[ -n "$DESCRIBE" ]]; then
  COMMON+=(--name "${NAME:-hgc-describe}")
  exec bash "$SUBMIT" "${COMMON[@]}" -- python "$HERE/describe-callset.py" "$DESCRIBE" --prefix "$PREFIX"
fi
[[ -n "$GVCF_DIR" && -n "$OUT" ]] || { echo "ERROR: --gvcf-dir and --out are required" >&2; exit 2; }
[[ "$GVCF_DIR" == /* && "$OUT" == /* ]] || { echo "ERROR: give absolute paths; the job runs from the log directory" >&2; exit 2; }
case "$BUILD" in GRCh38|GRCh37) ;; *) echo "ERROR: --build must be GRCh38 or GRCh37" >&2; exit 2 ;; esac
CMD=(hvantk hgc pipeline -i "$GVCF_DIR" -o "$OUT" --reference-genome "$BUILD" --output-prefix "$PREFIX")
[[ -n "$NPART" ]] && CMD+=(--n-partitions "$NPART")
(( PLAN )) && CMD+=(--dry-run)
# shellcheck disable=SC2206
[[ -n "$EXTRA" ]] && CMD+=($EXTRA)
COMMON+=(--name "${NAME:-hgc-$PREFIX}")
exec bash "$SUBMIT" "${COMMON[@]}" -- "${CMD[@]}"
