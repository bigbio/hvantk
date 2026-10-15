#!/usr/bin/env bash
# Submit a Slurm array job that builds many hvantk datasets, one per task.
#
#   submit-build-array.sh --list FILE [--name NAME] [--cpus N] [--mem-gb G]
#                         [--time D-HH:MM:SS] [--partition <P>] [--max-concurrent K]
#                         [--tasks i,j] [--sif <image>] [--dry-run] [--wait]
#
# FILE is a TSV, one build per line, absolute paths, NO blank or comment lines
# (build-array.sbatch reads task i as FILE's own physical line i via `sed -n ip`, so
# every physical line must be a real build):
#   <provider:dataset>\t<raw dir>\t<output path>[\t<extra reprocess options>]
# The fourth field, when present, is extra `hvantk reprocess` options
# (e.g. `--plugin-arg reference_genome=GRCh38`), word-split as written.
#
# Reads every site value from $HVANTK_SITE_PROFILE (skills/_site/SITE-PROFILE.md). One
# request sizes every task in the array: size it for the list's largest dataset, and
# give an unusually large dataset its own hpc-hail-job submission instead of a line
# here. --tasks is a literal sbatch array spec (e.g. "2,3" or "2-3"), not a convenience
# syntax of this script's own -- sbatch's own `--array=` rules apply. Prints the array
# job id. With --dry-run, runs `sbatch --test-only`. With --wait, polls accounting
# every 30s until every task has reached a terminal state, prints the per-task table,
# and exits 0 only if every task reached COMPLETED.
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
# Resolved to an absolute path before it is handed to sbatch: build-array.sbatch
# itself runs in the submission directory (Slurm's default when no --chdir is given,
# which this script does not pass), not the log directory -- only hail-job.sbatch,
# which it calls, changes into the log directory. A relative $LIST would therefore
# resolve against wherever sbatch happens to be invoked from, not against anything
# predictable inside the job.
LIST="$(cd "$(dirname "$LIST")" && pwd -P)/$(basename "$LIST")"
PARTITION="${PARTITION:-$HVANTK_CPU_PARTITION}"   # an empty fat-partition value falls back
if [[ "$PARTITION" == "$HVANTK_CPU_PARTITION" ]]; then
  (( CPUS <= HVANTK_NODE_CPUS )) || { echo "ERROR: --cpus $CPUS exceeds the CPU partition's $HVANTK_NODE_CPUS per node" >&2; exit 2; }
  (( MEM_GB <= HVANTK_NODE_MEM_GB )) || { echo "ERROR: --mem-gb $MEM_GB exceeds the CPU partition's $HVANTK_NODE_MEM_GB per node" >&2; exit 2; }
fi

# --- validate the list: no blank or comment line (every physical line becomes one
# array task; a skipped line here would desync "line N" from "task N"), >=3
# tab-separated fields, a provider:dataset shape, an existing absolute raw directory,
# an absolute output path not repeated on another line. Checked here (login node)
# because a bad line would otherwise fail deep inside one array task, or -- for a
# repeated output -- race silently between two tasks. ---
N=0
declare -A SEEN_OUTPUTS=()
while IFS=$'\t' read -r DATASET RAW OUT _EXTRA || [[ -n "${DATASET:-}" ]]; do
  N=$((N + 1))
  # A CRLF line ending lands on whichever field has no tab after it (OUT, or _EXTRA
  # when present); read -r never strips \r. Caught here, not deep inside the field
  # checks below, where a trailing \r would just make every message confusing.
  case "$DATASET$RAW$OUT${_EXTRA:-}" in
    *$'\r'*) echo "ERROR: $LIST line $N contains a carriage return (CRLF line ending); convert the file to LF line endings" >&2; exit 2 ;;
  esac
  [[ -n "$DATASET" && "$DATASET" != '#'* ]] || { echo "ERROR: $LIST line $N is blank or a comment; every physical line must be a build (array task N reads physical line N)" >&2; exit 2; }
  [[ -n "${RAW:-}" && -n "${OUT:-}" ]] || { echo "ERROR: $LIST line $N: need at least three tab-separated fields" >&2; exit 2; }
  [[ "$DATASET" =~ ^[a-z0-9_-]+:[a-z0-9_-]+$ ]] || { echo "ERROR: $LIST line $N: dataset must look like provider:dataset, got '$DATASET'" >&2; exit 2; }
  [[ "$RAW" == /* ]] || { echo "ERROR: $LIST line $N: raw dir must be an absolute path, got '$RAW'" >&2; exit 2; }
  [[ -d "$RAW" ]] || { echo "ERROR: $LIST line $N: raw dir does not exist: $RAW" >&2; exit 2; }
  [[ "$OUT" == /* ]] || { echo "ERROR: $LIST line $N: output must be an absolute path, got '$OUT'" >&2; exit 2; }
  [[ -z "${SEEN_OUTPUTS[$OUT]:-}" ]] || { echo "ERROR: $LIST line $N: output $OUT repeats line ${SEEN_OUTPUTS[$OUT]}; two tasks would race on writing it" >&2; exit 2; }
  SEEN_OUTPUTS[$OUT]=$N
done < "$LIST"
(( N > 0 )) || { echo "ERROR: $LIST has no build lines" >&2; exit 2; }

# --- validate --tasks before submitting: only digits, commas and dashes (a literal
# sbatch array spec, never a convenience syntax of this script's own), and every
# index or range endpoint between 1 and N -- sbatch would otherwise either reject the
# whole submission or, for an out-of-range index sbatch itself tolerates, silently
# queue a task that reads a line past the end of the list. ---
if [[ -n "$TASKS" ]]; then
  [[ "$TASKS" =~ ^[0-9,-]+$ ]] || { echo "ERROR: --tasks must contain only digits, commas and dashes (a literal sbatch array spec), got '$TASKS'" >&2; exit 2; }
  IFS=',' read -r -a TASK_PARTS <<< "$TASKS"
  for PART in "${TASK_PARTS[@]}"; do
    if [[ "$PART" == *-* ]]; then
      LO="${PART%%-*}"; HI="${PART##*-}"
    else
      LO="$PART"; HI="$PART"
    fi
    [[ "$LO" =~ ^[0-9]+$ && "$HI" =~ ^[0-9]+$ ]] || { echo "ERROR: --tasks '$PART' is not a valid index or range" >&2; exit 2; }
    (( LO >= 1 && LO <= N && HI >= 1 && HI <= N && LO <= HI )) || { echo "ERROR: --tasks index/range '$PART' is out of bounds 1-$N" >&2; exit 2; }
  done
fi
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

# Wait on accounting, not on squeue, exactly as submit-hail-job.sh does for a single
# job: an empty squeue answer can also be a transient RPC failure, and accounting is
# the record the verdict is graded from. Terminal states are an allow-list (prefix
# match, as submit-hail-job.sh uses), not an "in-flight" list: an unrecognized state
# then means "keep waiting", the safe default, rather than "done". A still-pending
# task can also be reported collapsed into one range row (<A>_[i-j]) before it is
# split into its own <A>_<i> row, so that shape is treated as in-flight on sight,
# whatever its State column says.
TERMINAL='^(COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|BOOT_FAIL|DEADLINE|REVOKED)'
STATES=""; UNKNOWN=0
while :; do
  sleep 30
  STATES="$(sacct -j "$JOB_ID" -X -n -P -o JobID,State 2>/dev/null || true)"
  if [[ -z "$STATES" ]]; then
    if [[ -z "$(squeue -j "$JOB_ID" -h -o %T 2>/dev/null)" ]]; then
      (( ++UNKNOWN < 10 )) || { echo "ERROR: array job $JOB_ID is unknown to both sacct and squeue; is accounting enabled?" >&2; exit 1; }
    fi
    continue
  fi
  UNKNOWN=0
  # Count-based, never `grep -q`/`grep -Eq` on a pipe: on a large array, grep's early
  # exit on the first match can send the upstream command (cut) a SIGPIPE, which
  # set -o pipefail then reports as the whole pipeline's exit status -- independent of
  # what grep itself found, and able to flip a true "still waiting" into a false "done".
  RANGE_ROWS=$(grep -c '_\[' <<< "$STATES" || true)
  NON_TERMINAL=$(cut -d'|' -f2 <<< "$STATES" | grep -vEc "$TERMINAL" || true)
  (( RANGE_ROWS == 0 && NON_TERMINAL == 0 )) && break
done
sacct -j "$JOB_ID" -X -o JobID,JobName%30,State,ExitCode,Elapsed,ReqMem,AllocCPUS >&2
echo "log: $LOG_DIR/${NAME}_${JOB_ID}_<task>.out (one file per task)" >&2
NOT_COMPLETED=$(cut -d'|' -f2 <<< "$STATES" | grep -vc '^COMPLETED$' || true)
(( NOT_COMPLETED == 0 ))
