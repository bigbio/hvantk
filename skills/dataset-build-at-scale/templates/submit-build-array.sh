#!/usr/bin/env bash
# Submit a Slurm array job that builds many hvantk datasets, one per task.
#
#   submit-build-array.sh --list FILE [--name NAME] [--cpus N] [--mem-gb G]
#                         [--time D-HH:MM:SS] [--partition <P>] [--max-concurrent K]
#                         [--tasks i,j] [--sif <image>] [--dry-run] [--wait]
#
# FILE is a TSV, one build per line, absolute paths:
#   <provider:dataset>\t<raw dir>\t<output path>[\t<extra reprocess options>]
# The fourth field, when present, is extra `hvantk reprocess` options
# (e.g. `--plugin-arg reference_genome=GRCh38`), word-split as written.
#
# Reads every site value from $HVANTK_SITE_PROFILE (skills/_site/SITE-PROFILE.md). One
# request sizes every task in the array: size it for the list's largest dataset, and
# give an unusually large dataset its own hpc-hail-job submission instead of a line
# here. Prints the array job id. With --dry-run, runs `sbatch --test-only`. With --wait,
# polls accounting every 30s until no task is PENDING/RUNNING/REQUEUED/COMPLETING,
# prints the per-task table, and exits 0 only if every task reached COMPLETED.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SITE_DIR="$(cd "$HERE/../../_site" && pwd -P)"
bash "$SITE_DIR/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
HAIL_JOB="$(cd "$HERE/../../hpc-hail-job/templates" && pwd -P)/hail-job.sbatch"
[[ -f "$HAIL_JOB" ]] || { echo "ERROR: hpc-hail-job's hail-job.sbatch not found next to this skill: $HAIL_JOB" >&2; exit 2; }

LIST=""; NAME="hvantk-build"; CPUS=4; MEM_GB=16; TIME="0-00:30:00"
PARTITION="$HVANTK_CPU_PARTITION"; MAX_CONCURRENT=4; TASKS=""; SIF_OVERRIDE=""
DRY=0; WAIT=0
while (( $# )); do
  case "$1" in
    --list) LIST="$2"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem-gb) MEM_GB="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --max-concurrent) MAX_CONCURRENT="$2"; shift 2 ;;
    --tasks) TASKS="$2"; shift 2 ;;
    --sif) SIF_OVERRIDE="$2"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    --wait) WAIT=1; shift ;;
    *) echo "ERROR: unknown option $1" >&2; exit 2 ;;
  esac
done
[[ -n "$LIST" ]] || { echo "ERROR: --list FILE is required" >&2; exit 2; }
[[ -f "$LIST" ]] || { echo "ERROR: task list not found: $LIST" >&2; exit 2; }
PARTITION="${PARTITION:-$HVANTK_CPU_PARTITION}"   # an empty fat-partition value falls back
if [[ "$PARTITION" == "$HVANTK_CPU_PARTITION" ]]; then
  (( CPUS <= HVANTK_NODE_CPUS )) || { echo "ERROR: --cpus $CPUS exceeds the CPU partition's $HVANTK_NODE_CPUS per node" >&2; exit 2; }
  (( MEM_GB <= HVANTK_NODE_MEM_GB )) || { echo "ERROR: --mem-gb $MEM_GB exceeds the CPU partition's $HVANTK_NODE_MEM_GB per node" >&2; exit 2; }
fi

# --- validate the list: non-empty lines, >=3 tab-separated fields, a provider:dataset
# shape, an existing absolute raw directory, an absolute output path. Checked here
# (login node) because a bad line would otherwise fail deep inside one array task. ---
N=0
while IFS=$'\t' read -r DATASET RAW OUT _EXTRA || [[ -n "${DATASET:-}" ]]; do
  [[ -n "${DATASET:-}" ]] || continue   # tolerate a trailing blank line
  N=$((N + 1))
  [[ -n "${RAW:-}" && -n "${OUT:-}" ]] || { echo "ERROR: $LIST line $N: need at least three tab-separated fields" >&2; exit 2; }
  [[ "$DATASET" =~ ^[a-z0-9_-]+:[a-z0-9_-]+$ ]] || { echo "ERROR: $LIST line $N: dataset must look like provider:dataset, got '$DATASET'" >&2; exit 2; }
  [[ "$RAW" == /* ]] || { echo "ERROR: $LIST line $N: raw dir must be an absolute path, got '$RAW'" >&2; exit 2; }
  [[ -d "$RAW" ]] || { echo "ERROR: $LIST line $N: raw dir does not exist: $RAW" >&2; exit 2; }
  [[ "$OUT" == /* ]] || { echo "ERROR: $LIST line $N: output must be an absolute path, got '$OUT'" >&2; exit 2; }
done < "$LIST"
(( N > 0 )) || { echo "ERROR: $LIST has no build lines" >&2; exit 2; }

ARRAY_SPEC="${TASKS:-1-$N}"
LOG_DIR="${HVANTK_LOG_DIR:-$HVANTK_WORK/hvantk-logs}"
mkdir -p "$LOG_DIR"

ARGS=(--parsable --job-name="$NAME" --partition="$PARTITION"
      --cpus-per-task="$CPUS" --mem="${MEM_GB}G" --time="$TIME"
      --array="${ARRAY_SPEC}%${MAX_CONCURRENT}"
      --output="$LOG_DIR/%x_%A_%a.out" --export=ALL)
if [[ -n "$SIF_OVERRIDE" ]]; then
  [[ "$SIF_OVERRIDE" == /* && -e "$SIF_OVERRIDE" ]] || { echo "ERROR: --sif needs an absolute path to an existing image: $SIF_OVERRIDE" >&2; exit 2; }
fi
export HVANTK_SIF_OVERRIDE="$SIF_OVERRIDE"   # exported even when empty, so a leftover value cannot redirect a job
[[ -n "${HVANTK_ACCOUNT:-}" ]] && ARGS+=(--account="$HVANTK_ACCOUNT")
[[ -n "${HVANTK_QOS:-}" ]] && ARGS+=(--qos="$HVANTK_QOS")
[[ -n "${SLURM_EMAIL_NOTIFICATION:-}" ]] && ARGS+=(--mail-type=END,FAIL --mail-user="$SLURM_EMAIL_NOTIFICATION")

echo "request: name=$NAME partition=$PARTITION cpus=$CPUS mem=${MEM_GB}G time=$TIME array=${ARRAY_SPEC}%${MAX_CONCURRENT} tasks=$N image=${SIF_OVERRIDE:-profile}" >&2
echo "list: $LIST" >&2
if (( DRY )); then
  sbatch --test-only "${ARGS[@]}" "$HERE/build-array.sbatch" "$LIST" "$HAIL_JOB"
  exit 0
fi
JOB_ID="$(sbatch "${ARGS[@]}" "$HERE/build-array.sbatch" "$LIST" "$HAIL_JOB")"
JOB_ID="${JOB_ID%%;*}"
echo "$JOB_ID"
(( WAIT )) || exit 0

# Wait on accounting, not on squeue (see hpc-hail-job/templates/submit-hail-job.sh for
# why). A still-pending task can appear as one collapsed row (state PENDING) before it
# is assigned its own <A>_<i> id, so the in-flight check matches on state, not shape.
IN_FLIGHT='^(PENDING|RUNNING|REQUEUED|COMPLETING)$'
STATES=""; UNKNOWN=0
while :; do
  sleep 30
  STATES="$(sacct -j "$JOB_ID" -X -n -P -o JobID,State 2>/dev/null || true)"
  if [[ -z "$STATES" ]]; then
    (( ++UNKNOWN < 10 )) || { echo "ERROR: array job $JOB_ID has no accounting rows after repeated polling; is accounting enabled?" >&2; exit 1; }
    continue
  fi
  UNKNOWN=0
  cut -d'|' -f2 <<< "$STATES" | grep -Eq "$IN_FLIGHT" || break
done
sacct -j "$JOB_ID" -X -o JobID,JobName%30,State,ExitCode,Elapsed,ReqMem,AllocCPUS >&2
echo "log: $LOG_DIR/${NAME}_${JOB_ID}_<task>.out (one file per task)" >&2
! cut -d'|' -f2 <<< "$STATES" | grep -qv '^COMPLETED$'
