#!/bin/bash
# Run hvantk from the Apptainer/Singularity image with Spark scratch set up correctly.
#
# Without SPARK_LOCAL_DIRS on a writable NODE-LOCAL dir, Hail init inside the container
# fails with "DiskBlockManager: Failed to create any local dir" followed by a misleading
# "[Errno 111] Connection refused" from py4j. Always go through this wrapper.
#
#   bash hvantk_run.sh utils check-install
#   bash hvantk_run.sh reprocess clinvar:variants --raw-dir ... --output ...
#   HVANTK_SIF=/path/to/hvantk.sif bash hvantk_run.sh --help
set -euo pipefail

HVANTK_SIF="${HVANTK_SIF:-${WORK:?WORK unset}/containers/hvantk.sif}"
[[ -f "$HVANTK_SIF" ]] || { echo "ERROR: image not found: $HVANTK_SIF" >&2; exit 1; }

# node-local scratch, never home/Lustre/GPFS (guide section 4.1)
SPARK_SCRATCH="${SLURM_TMPDIR:-/tmp/spark-$USER-${SLURM_JOB_ID:-$$}}"
mkdir -p "$SPARK_SCRATCH"
cleanup(){ rm -rf "$SPARK_SCRATCH"; }
trap cleanup EXIT

# Bind any extra data roots the command needs, e.g. HVANTK_BIND="$WORK:$WORK"
BINDS=(-B "$SPARK_SCRATCH:$SPARK_SCRATCH")
[[ -n "${HVANTK_BIND:-}" ]] && BINDS+=(-B "$HVANTK_BIND")

exec singularity exec "${BINDS[@]}" \
  --env SPARK_LOCAL_DIRS="$SPARK_SCRATCH" \
  --env TMPDIR="$SPARK_SCRATCH" \
  --env HAIL_TMPDIR="$SPARK_SCRATCH" \
  "$HVANTK_SIF" hvantk "$@"
