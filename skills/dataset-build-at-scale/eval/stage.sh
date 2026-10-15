#!/usr/bin/env bash
# Stage the public inputs for the dataset-build-at-scale eval, once per cluster.
# Login-node work: it downloads, nothing more. Writes staged.sha256 next to the files.
#
#   HVANTK_SITE_PROFILE=... bash skills/dataset-build-at-scale/eval/stage.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/dataset-build-at-scale"
RAW="$STAGE/raw"
mkdir -p "$RAW/hgnc" "$RAW/clingen" "$RAW/gencc"
CONTAINER_BIN="${HVANTK_CONTAINER_BIN:-singularity}"
BINDS=(-B "$HVANTK_WORK:$HVANTK_WORK")
for p in "$HVANTK_REPO" "$STAGE"; do
  case "$p" in "$HVANTK_WORK"|"$HVANTK_WORK"/*) ;; *) BINDS+=(-B "$p:$p") ;; esac
done
run() { "$CONTAINER_BIN" exec "${BINDS[@]}" --env PYTHONPATH="$HVANTK_REPO" "$HVANTK_SIF" "$@"; }

# Every raw dir must hold exactly one file matching *pattern* (a bash glob, case
# `*)` syntax) and nothing else. A failed or partial download does not disappear; it
# leaves a file behind under a different name (stream_to_file's own temp-file
# convention is <name>.<random>), which this would otherwise silently accept as "the"
# raw file on the next run.
validate_raw_dir() {
  local dir="$1" pattern="$2" label="$3"
  local -a matching=() extra=()
  local f
  while IFS= read -r -d '' f; do
    case "$(basename "$f")" in
      $pattern) matching+=("$f") ;;
      *) extra+=("$f") ;;
    esac
  done < <(find "$dir" -maxdepth 1 -type f ! -name '.*' -print0)
  if (( ${#extra[@]} > 0 )); then
    echo "ERROR: $label raw dir ($dir) holds an unexpected file: ${extra[0]} (a partial download looks like <name>.<random>; remove it and re-stage)" >&2
    exit 1
  fi
  if (( ${#matching[@]} != 1 )); then
    echo "ERROR: $label raw dir ($dir) must hold exactly one file matching '$pattern'; found ${#matching[@]}: ${matching[*]:-(none)}" >&2
    exit 1
  fi
}

# Each plugin's own downloader from the checkout, so the URL has one home
# (hvantk/skills/<provider>/shared) and the checkout's code is what runs. Skip a
# provider whose directory already holds its expected file: re-staging on another day
# changes ClinGen's and GenCC's dated filename, which would move the reference out
# from under a control that has not been re-run with WRITE_REFERENCE=1.
if [[ -s "$RAW/hgnc/hgnc_complete_set.txt" ]]; then
  echo "already staged: $RAW/hgnc/hgnc_complete_set.txt"
else
  run hvantk download hgnc --output "$RAW/hgnc/hgnc_complete_set.txt"
fi
validate_raw_dir "$RAW/hgnc" "hgnc_complete_set.txt" hgnc

shopt -s nullglob
clingen_existing=("$RAW"/clingen/*.csv)
gencc_existing=("$RAW"/gencc/*.tsv)
shopt -u nullglob

if (( ${#clingen_existing[@]} > 0 )); then
  echo "already staged: $RAW/clingen"
else
  run hvantk download clingen --output-dir "$RAW/clingen"
fi
validate_raw_dir "$RAW/clingen" "*.csv" clingen

if (( ${#gencc_existing[@]} > 0 )); then
  echo "already staged: $RAW/gencc"
else
  run hvantk download gencc --output-dir "$RAW/gencc"
fi
validate_raw_dir "$RAW/gencc" "*.tsv" gencc

( cd "$STAGE" && sha256sum raw/*/* > staged.sha256 )
cat "$STAGE/staged.sha256"
echo "staged under $STAGE"
