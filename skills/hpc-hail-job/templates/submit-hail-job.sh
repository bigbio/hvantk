#!/usr/bin/env bash
# Submit one hvantk command as a single-node Hail job.
#
#   submit-hail-job.sh [--name NAME] [--cpus N] [--mem-gb G] [--time D-HH:MM:SS]
#                      [--partition <P>] [--exclusive] [--dry-run] [--wait]
#                      -- <hvantk ...|python ...>
#
# Reads every site value from $HVANTK_SITE_PROFILE (skills/_site/SITE-PROFILE.md).
# Prints the job id. With --wait, polls accounting until the job reaches a terminal
# state and exits 0 only if that state is COMPLETED. With --dry-run, runs
# `sbatch --test-only`. The node caps (HVANTK_NODE_CPUS, HVANTK_NODE_MEM_GB) describe
# the CPU partition and are enforced only for it; a fat partition sizes itself.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SITE_DIR="$(cd "$HERE/../../_site" && pwd -P)"
bash "$SITE_DIR/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"

NAME="hvantk-hail"; CPUS=4; MEM_GB=16; TIME="0-01:00:00"; PARTITION="$HVANTK_CPU_PARTITION"
EXCLUSIVE=0; DRY=0; WAIT=0
while (( $# )); do
  case "$1" in
    --name) NAME="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem-gb) MEM_GB="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --exclusive) EXCLUSIVE=1; shift ;;
    --dry-run) DRY=1; shift ;;
    --wait) WAIT=1; shift ;;
    --) shift; break ;;
    *) echo "ERROR: unknown option $1 (put the command after --)" >&2; exit 2 ;;
  esac
done
(( $# > 0 )) || { echo "ERROR: no command after --" >&2; exit 2; }
PARTITION="${PARTITION:-$HVANTK_CPU_PARTITION}"   # an empty fat-partition value falls back
case "$1" in hvantk|python|python3) ;; *) echo "ERROR: command must start with hvantk or python, got '$1'" >&2; exit 2 ;; esac
if [[ "$PARTITION" == "$HVANTK_CPU_PARTITION" ]]; then
  (( CPUS <= HVANTK_NODE_CPUS )) || { echo "ERROR: --cpus $CPUS exceeds the CPU partition's $HVANTK_NODE_CPUS per node" >&2; exit 2; }
  (( MEM_GB <= HVANTK_NODE_MEM_GB )) || { echo "ERROR: --mem-gb $MEM_GB exceeds the CPU partition's $HVANTK_NODE_MEM_GB per node" >&2; exit 2; }
fi

LOG_DIR="${HVANTK_LOG_DIR:-$HVANTK_WORK/hvantk-logs}"
mkdir -p "$LOG_DIR"

ARGS=(--parsable --job-name="$NAME" --partition="$PARTITION"
      --cpus-per-task="$CPUS" --mem="${MEM_GB}G" --time="$TIME"
      --output="$LOG_DIR/%x_%j.out" --export=ALL)
[[ -n "${HVANTK_ACCOUNT:-}" ]] && ARGS+=(--account="$HVANTK_ACCOUNT")
[[ -n "${HVANTK_QOS:-}" ]] && ARGS+=(--qos="$HVANTK_QOS")
(( EXCLUSIVE )) && ARGS+=(--exclusive)
[[ -n "${SLURM_EMAIL_NOTIFICATION:-}" ]] && ARGS+=(--mail-type=END,FAIL --mail-user="$SLURM_EMAIL_NOTIFICATION")

echo "request: name=$NAME partition=$PARTITION cpus=$CPUS mem=${MEM_GB}G time=$TIME exclusive=$EXCLUSIVE" >&2
echo "command: $*" >&2
if (( DRY )); then
  sbatch --test-only "${ARGS[@]}" "$HERE/hail-job.sbatch" "$@"
  exit 0
fi
JOB_ID="$(sbatch "${ARGS[@]}" "$HERE/hail-job.sbatch" "$@")"
JOB_ID="${JOB_ID%%;*}"
echo "$JOB_ID"
(( WAIT )) || exit 0

# Wait on accounting, not on squeue: an empty squeue answer can also be a transient
# RPC failure, and accounting is the record the verdict is graded from.
TERMINAL='^(COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|BOOT_FAIL|DEADLINE|REVOKED)'
STATE=""; UNKNOWN=0
until [[ "$STATE" =~ $TERMINAL ]]; do
  sleep 30
  STATE="$(sacct -j "$JOB_ID" -X -n -o State 2>/dev/null | head -1 | tr -d ' ' || true)"
  if [[ -z "$STATE" && -z "$(squeue -j "$JOB_ID" -h -o %T 2>/dev/null)" ]]; then
    (( ++UNKNOWN < 10 )) || { echo "ERROR: job $JOB_ID is unknown to both sacct and squeue; is accounting enabled?" >&2; exit 1; }
  else
    UNKNOWN=0
  fi
done
sacct -j "$JOB_ID" -X -o JobID,JobName%30,State,ExitCode,Elapsed,ReqMem,AllocCPUS >&2
echo "log: $LOG_DIR/${NAME}_${JOB_ID}.out" >&2
[[ "$STATE" == "COMPLETED" ]]
