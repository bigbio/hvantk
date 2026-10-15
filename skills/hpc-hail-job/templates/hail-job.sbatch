#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=hvantk-hail
# Resources, partition, account, QOS, time, output and mail are passed on the command
# line by submit-hail-job.sh. No site value is written here.
#
# Runs ONE command inside the hvantk container image with Hail in Spark local mode:
#   node-local scratch for Spark spill and Hail temp, Spark pinned to the allocation,
#   driver memory at 80 percent of the request, the work root bound into the image.
# PYTHONPATH points at the checkout, so the job runs the checkout's code inside the
# image's environment: the commit is what runs, the image is what it runs on.
# The command is the positional arguments and starts with `hvantk` or `python`. The
# working directory is the log directory, so Hail's own hail-*.log lands there; pass
# absolute paths in the command.
set -euo pipefail

: "${HVANTK_SITE_PROFILE:?HVANTK_SITE_PROFILE not set}"
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
: "${HVANTK_REPO:?}" "${HVANTK_SIF:?}" "${HVANTK_WORK:?}"
# An image other than the profile's default, for verifying a freshly built one: the
# submit script passes it as HVANTK_SIF_OVERRIDE, because sourcing the profile above
# would otherwise reset HVANTK_SIF.
HVANTK_SIF="${HVANTK_SIF_OVERRIDE:-$HVANTK_SIF}"
[[ -e "$HVANTK_SIF" ]] || { echo "ERROR: image not found: $HVANTK_SIF" >&2; exit 2; }
(( $# > 0 )) || { echo "ERROR: no command given (expected: hvantk ... | python ...)" >&2; exit 2; }

CONTAINER_BIN="${HVANTK_CONTAINER_BIN:-singularity}"
command -v "$CONTAINER_BIN" >/dev/null || { echo "ERROR: $CONTAINER_BIN not on PATH" >&2; exit 2; }
LOG_DIR="${HVANTK_LOG_DIR:-$HVANTK_WORK/hvantk-logs}"
mkdir -p "$LOG_DIR"

# --- node-local scratch (never home, never the parallel filesystem) ---
SCRATCH_ROOT="${HVANTK_SCRATCH_ROOT:-${SLURM_TMPDIR:-/tmp}}"
SCRATCH="$SCRATCH_ROOT/hvantk-$USER-${SLURM_JOB_ID:-$$}"
mkdir -p "$SCRATCH"
cleanup() { rm -rf "$SCRATCH"; }
trap cleanup EXIT   # no `exec` below, so this trap does run

# --- cap Spark at the allocation ---
CPUS="${SLURM_CPUS_PER_TASK:-1}"
MEM_MB="${SLURM_MEM_PER_NODE:-4096}"
DRIVER_MB=$(( MEM_MB * 80 / 100 ))
PYSPARK_SUBMIT_ARGS="--master local[${CPUS}] --driver-memory ${DRIVER_MB}m --conf spark.local.dir=${SCRATCH} pyspark-shell"
JAVA_TOOL_OPTIONS="-Duser.language=en -Duser.country=US"

# --- binds: scratch, the work root, and anything outside the work root the job needs ---
BINDS=(-B "$SCRATCH:$SCRATCH" -B "$HVANTK_WORK:$HVANTK_WORK")
bind_outside_work() {
  local p="$1"
  [[ -d "$p" ]] || return 0
  case "$p" in "$HVANTK_WORK"|"$HVANTK_WORK"/*) ;; *) BINDS+=(-B "$p:$p") ;; esac
}
bind_outside_work "$HVANTK_REPO"
[[ -n "${HVANTK_EVAL_DATA:-}" ]] && bind_outside_work "$HVANTK_EVAL_DATA"
[[ -n "${HVANTK_BIND:-}" ]] && BINDS+=(-B "$HVANTK_BIND")

[[ -d "$HVANTK_REPO/hvantk" ]] || { echo "ERROR: checkout $HVANTK_REPO is not visible on $(hostname); the job would silently run the image's own copy" >&2; exit 2; }
cd "$LOG_DIR"
echo "== hvantk hail job =="
echo "date:      $(date -Is)"
echo "host:      $(hostname)"
echo "job:       ${SLURM_JOB_ID:-?}  cpus=$CPUS mem_mb=$MEM_MB driver_mb=$DRIVER_MB"
echo "image:     $(readlink -f "$HVANTK_SIF")"
echo "repo:      $HVANTK_REPO"
echo "scratch:   $SCRATCH"
echo "logs:      $LOG_DIR"
echo "command:   $*"
echo "====================="

# The image's %environment resets thread counts and proxies, so every value the job
# depends on goes through --env.
"$CONTAINER_BIN" exec "${BINDS[@]}" \
  --env SPARK_LOCAL_DIRS="$SCRATCH" --env TMPDIR="$SCRATCH" --env HAIL_TMPDIR="$SCRATCH" \
  --env PYSPARK_SUBMIT_ARGS="$PYSPARK_SUBMIT_ARGS" \
  --env JAVA_TOOL_OPTIONS="$JAVA_TOOL_OPTIONS" \
  --env PYTHONPATH="$HVANTK_REPO" --env HAIL_LOG_DIR="$LOG_DIR" \
  --env OMP_NUM_THREADS="$CPUS" --env OPENBLAS_NUM_THREADS="$CPUS" \
  --env MKL_NUM_THREADS="$CPUS" --env NUMEXPR_NUM_THREADS="$CPUS" \
  "$HVANTK_SIF" "$@"
