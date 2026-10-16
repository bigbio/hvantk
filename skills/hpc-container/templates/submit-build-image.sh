#!/usr/bin/env bash
# Submit an image build from the login node, where git exists, so the job carries the
# commit, branch and tree state it builds from.
#
#   submit-build-image.sh --out <image.sif> [--def <def file>] [--name NAME]
#                         [--cpus 8] [--mem-gb 32] [--time 0-01:45:00] [--partition <P>]
#                         [--dry-run] [--wait]
#
# Prints the job id. With --wait, polls accounting until the job reaches a terminal
# state and exits 0 only if it is COMPLETED.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" --allow-missing-image >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
command -v git >/dev/null || { echo "ERROR: git is needed on the submit host to record provenance" >&2; exit 2; }

OUT=""; DEF="$HVANTK_REPO/containers/hvantk.def"; NAME="hvantk-build-image"
CPUS=8; MEM_GB=32; TIME="0-01:45:00"; PARTITION="$HVANTK_CPU_PARTITION"; DRY=0; WAIT=0
while (( $# )); do
  case "$1" in
    --out) OUT="$2"; shift 2 ;;
    --def) DEF="$2"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem-gb) MEM_GB="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    --wait) WAIT=1; shift ;;
    *) echo "ERROR: unknown option $1" >&2; exit 2 ;;
  esac
done
PARTITION="${PARTITION:-$HVANTK_CPU_PARTITION}"
[[ -n "$OUT" ]] || { echo "ERROR: --out <image.sif> is required" >&2; exit 2; }
[[ "$OUT" == *.sif ]] || { echo "ERROR: --out must end in .sif" >&2; exit 2; }
[[ -f "$DEF" ]] || { echo "ERROR: def file not found: $DEF" >&2; exit 2; }
[[ "$OUT" == /* && "$DEF" == /* ]] || { echo "ERROR: give absolute paths for --out and --def; the job runs from the checkout" >&2; exit 2; }
case "$(basename "$OUT")" in hvantk.sif) echo "ERROR: never build straight onto the default name; build a versioned file and promote it" >&2; exit 2 ;; esac
if [[ -e "$HVANTK_SIF" && "$(readlink -f "$OUT" 2>/dev/null)" == "$(readlink -f "$HVANTK_SIF")" ]]; then
  echo "ERROR: $OUT is the file the profile's default image resolves to; build under another name and promote it" >&2; exit 2
fi

export BUILD_OUT="$OUT" BUILD_DEF="$DEF"
BUILD_GIT_SHA="$(git -C "$HVANTK_REPO" rev-parse HEAD)" || { echo "ERROR: cannot read the checkout's commit" >&2; exit 2; }
BUILD_GIT_BRANCH="$(git -C "$HVANTK_REPO" rev-parse --abbrev-ref HEAD)"
STATUS="$(git -C "$HVANTK_REPO" status --porcelain)" || { echo "ERROR: cannot read the checkout's status" >&2; exit 2; }
BUILD_GIT_DIRTY="$([ -z "$STATUS" ] && echo clean || echo dirty)"
BUILD_LOCK_MD5="$(md5sum "$HVANTK_REPO/poetry.lock" | cut -d' ' -f1)"
export BUILD_GIT_SHA BUILD_GIT_BRANCH BUILD_GIT_DIRTY BUILD_LOCK_MD5
LOG_DIR="${HVANTK_LOG_DIR:-$HVANTK_WORK/hvantk-logs}"
mkdir -p "$LOG_DIR"

ARGS=(--parsable --job-name="$NAME" --partition="$PARTITION"
      --cpus-per-task="$CPUS" --mem="${MEM_GB}G" --time="$TIME"
      --output="$LOG_DIR/%x_%j.out" --export=ALL)
[[ -n "${HVANTK_ACCOUNT:-}" ]] && ARGS+=(--account="$HVANTK_ACCOUNT")
[[ -n "${HVANTK_QOS:-}" ]] && ARGS+=(--qos="$HVANTK_QOS")
[[ -n "${SLURM_EMAIL_NOTIFICATION:-}" ]] && ARGS+=(--mail-type=END,FAIL --mail-user="$SLURM_EMAIL_NOTIFICATION")

echo "build:   $BUILD_GIT_SHA ($BUILD_GIT_BRANCH, $BUILD_GIT_DIRTY) -> $OUT" >&2
echo "request: name=$NAME partition=$PARTITION cpus=$CPUS mem=${MEM_GB}G time=$TIME" >&2
if (( DRY )); then
  sbatch --test-only "${ARGS[@]}" "$HERE/build-image.sbatch"
  exit 0
fi
JOB_ID="$(sbatch "${ARGS[@]}" "$HERE/build-image.sbatch")"
JOB_ID="${JOB_ID%%;*}"
echo "$JOB_ID"
(( WAIT )) || exit 0

TERMINAL='^(COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|BOOT_FAIL|DEADLINE|REVOKED)'
STATE=""; UNKNOWN=0
until [[ "$STATE" =~ $TERMINAL ]]; do
  sleep 30
  STATE="$(sacct -j "$JOB_ID" -X -n -o State 2>/dev/null | head -1 | tr -d ' ' || true)"
  if [[ -z "$STATE" && -z "$(squeue -j "$JOB_ID" -h -o %T 2>/dev/null)" ]]; then
    (( ++UNKNOWN < 10 )) || { echo "ERROR: job $JOB_ID is unknown to both sacct and squeue" >&2; exit 1; }
  else
    UNKNOWN=0
  fi
done
sacct -j "$JOB_ID" -X -o JobID,JobName%30,State,ExitCode,Elapsed,ReqMem,AllocCPUS >&2
echo "log: $LOG_DIR/${NAME}_${JOB_ID}.out" >&2
[[ "$STATE" == "COMPLETED" ]]
