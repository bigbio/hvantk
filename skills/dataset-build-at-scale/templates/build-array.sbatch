#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=hvantk-build
# One array task builds one dataset. Resources, partition, array spec, account, QOS,
# time, output and mail are passed on the command line by submit-build-array.sh. No
# site value is written here.
#
# $1 is the task list (TSV: <provider:dataset>\t<raw dir>\t<output>[\t<extra options>]);
# $2 is the path to hpc-hail-job's hail-job.sbatch, which runs the actual build inside
# the container image. Both arrive as arguments because $0 inside an sbatch script is
# Slurm's own spooled copy, not a path this script can find siblings from.
set -euo pipefail
LIST="${1:?task list required}"
HAIL_JOB="${2:?path to hail-job.sbatch required}"
[[ -f "$LIST" ]] || { echo "ERROR: task list not found: $LIST" >&2; exit 2; }
[[ -f "$HAIL_JOB" ]] || { echo "ERROR: hail-job.sbatch not found: $HAIL_JOB" >&2; exit 2; }
: "${SLURM_ARRAY_TASK_ID:?this script runs as an array task (sbatch --array=...)}"

LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")"
[[ -n "$LINE" ]] || { echo "ERROR: no line $SLURM_ARRAY_TASK_ID in $LIST" >&2; exit 2; }
IFS=$'\t' read -r DATASET RAW OUT EXTRA <<< "$LINE"
[[ -n "$DATASET" && -n "$RAW" && -n "$OUT" ]] || { echo "ERROR: line $SLURM_ARRAY_TASK_ID of $LIST is missing a field: $LINE" >&2; exit 2; }

echo "== build task $SLURM_ARRAY_TASK_ID (array ${SLURM_ARRAY_JOB_ID:-?}, job ${SLURM_JOB_ID:-?}): $DATASET =="
echo "raw:    $RAW"
echo "output: $OUT"
echo "extra:  ${EXTRA:-(none)}"
echo "==============================================="
# EXTRA is deliberately unquoted: it is the task-list author's own extra reprocess
# options (e.g. --plugin-arg k=v --intermediate /abs/path), meant to word-split the way
# they wrote it, not a single opaque token. The split is on plain whitespace with no
# quote handling at all (there is no shell re-parse here, just bash's own word
# splitting): a value containing a space -- `--plugin-arg key="a value with spaces"` --
# cannot be expressed in this column. Keep every value in EXTRA space-free.
bash "$HAIL_JOB" hvantk reprocess "$DATASET" --raw-dir "$RAW" --output "$OUT" --skip-download --no-check-drift $EXTRA
echo "== done $DATASET =="
