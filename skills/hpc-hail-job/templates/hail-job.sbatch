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
# The command is the positional arguments and starts with `hvantk` or `python`.
set -euo pipefail

: "${HVANTK_SITE_PROFILE:?HVANTK_SITE_PROFILE not set}"
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
: "${HVANTK_REPO:?}" "${HVANTK_SIF:?}" "${HVANTK_WORK:?}"
(( $# > 0 )) || { echo "ERROR: no command given (expected: hvantk ... | python ...)" >&2; exit 2; }

CONTAINER_BIN="${HVANTK_CONTAINER_BIN:-singularity}"
command -v "$CONTAINER_BIN" >/dev/null || { echo "ERROR: $CONTAINER_BIN not on PATH" >&2; exit 2; }

# --- node-local scratch (never home, never the parallel filesystem) ---
SCRATCH_ROOT="${HVANTK_SCRATCH_ROOT:-${SLURM_TMPDIR:-/tmp}}"
SCRATCH="$SCRATCH_ROOT/hvantk-$USER-${SLURM_JOB_ID:-$$}"
mkdir -p "$SCRATCH"
cleanup() { rm -rf "$SCRATCH"; }
trap cleanup EXIT
export SPARK_LOCAL_DIRS="$SCRATCH" TMPDIR="$SCRATCH" HAIL_TMPDIR="$SCRATCH"

# --- cap Spark at the allocation ---
CPUS="${SLURM_CPUS_PER_TASK:-1}"
MEM_MB="${SLURM_MEM_PER_NODE:-4096}"
DRIVER_MB=$(( MEM_MB * 80 / 100 ))
export PYSPARK_SUBMIT_ARGS="--master local[${CPUS}] --driver-memory ${DRIVER_MB}m --conf spark.local.dir=${SCRATCH} pyspark-shell"
export OMP_NUM_THREADS="$CPUS" OPENBLAS_NUM_THREADS="$CPUS" MKL_NUM_THREADS="$CPUS" NUMEXPR_NUM_THREADS="$CPUS"
export JAVA_TOOL_OPTIONS="-Duser.language=en -Duser.country=US"
export NO_PROXY="localhost,127.0.0.1,0.0.0.0,::1${NO_PROXY:+,$NO_PROXY}"

# --- binds: scratch, the work root, and whatever the profile adds ---
BINDS=(-B "$SCRATCH:$SCRATCH" -B "$HVANTK_WORK:$HVANTK_WORK" -B "$HVANTK_REPO:$HVANTK_REPO")
[[ -n "${HVANTK_EVAL_DATA:-}" ]] && BINDS+=(-B "$HVANTK_EVAL_DATA:$HVANTK_EVAL_DATA")
[[ -n "${HVANTK_BIND:-}" ]] && BINDS+=(-B "$HVANTK_BIND")

cd "$HVANTK_REPO"
echo "== hvantk hail job =="
echo "date:      $(date -Is)"
echo "host:      $(hostname)"
echo "job:       ${SLURM_JOB_ID:-?}  cpus=$CPUS mem_mb=$MEM_MB driver_mb=$DRIVER_MB"
echo "image:     $(readlink -f "$HVANTK_SIF")"
echo "repo:      $HVANTK_REPO"
echo "scratch:   $SCRATCH"
echo "command:   $*"
echo "====================="

exec "$CONTAINER_BIN" exec "${BINDS[@]}" \
  --env SPARK_LOCAL_DIRS="$SCRATCH" --env TMPDIR="$SCRATCH" --env HAIL_TMPDIR="$SCRATCH" \
  --env PYSPARK_SUBMIT_ARGS="$PYSPARK_SUBMIT_ARGS" \
  --env JAVA_TOOL_OPTIONS="$JAVA_TOOL_OPTIONS" \
  "$HVANTK_SIF" "$@"
