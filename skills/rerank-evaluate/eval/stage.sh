#!/usr/bin/env bash
# Stage the rerank-evaluate eval, once per cluster. Login-node work: picks the newest
# PASSING cohort-attach-burden control (its gene table and gene index copy), and
# downloads the ClinGen and HGNC raw files through the image. Writes staged.sha256.
#
#   HVANTK_SITE_PROFILE=... bash skills/rerank-evaluate/eval/stage.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"

EVAL_ROOT="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}"
STAGE="$EVAL_ROOT/rerank-evaluate"
CAB="$EVAL_ROOT/cohort-attach-burden"
mkdir -p "$STAGE/cohort" "$STAGE/upstream" "$STAGE/raw/clingen" "$STAGE/raw/hgnc"
CONTAINER_BIN="${HVANTK_CONTAINER_BIN:-singularity}"
BINDS=(-B "$HVANTK_WORK:$HVANTK_WORK")
for p in "$HVANTK_REPO" "$STAGE"; do
  case "$p" in "$HVANTK_WORK"|"$HVANTK_WORK"/*) ;; *) BINDS+=(-B "$p:$p") ;; esac
done
run() { "$CONTAINER_BIN" exec "${BINDS[@]}" --env PYTHONPATH="$HVANTK_REPO" "$HVANTK_SIF" "$@"; }

# Every raw dir must hold exactly one file matching *pattern* and nothing else (copied
# from dataset-build-at-scale/eval/stage.sh: a partial download leaves a
# <name>.<random> temp file behind rather than disappearing).
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

# --- the newest cohort-attach-burden control directory whose positive arm passed ---
PICK="$(python3 -c '
import json, os, sys
root = sys.argv[1]
if os.path.isdir(root):
    for name in sorted(os.listdir(root), reverse=True):
        verdict = os.path.join(root, name, "pos", "verdict.json")
        if not os.path.isfile(verdict):
            continue
        try:
            v = json.load(open(verdict))
        except ValueError:
            continue
        if v.get("pass") is True:
            print(os.path.join(root, name))
            sys.exit(0)
sys.exit(1)
' "$CAB/control")" || {
  echo "ERROR: no cohort-attach-burden control directory under $CAB/control has a passing positive control; run skills/cohort-attach-burden/eval/control.sh first" >&2
  exit 1
}
echo "using cohort-attach-burden control: $PICK"

cp "$PICK/pos/gene_burden.tsv" "$STAGE/cohort/gene_burden.tsv"
# The upstream manifest's `table:` field names wherever ITS control wrote the burden
# TSV (an absolute path there too). Rewrite it to THIS stage's own absolute copy --
# hvantk resolves a manifest's `table:` against the reading process's working
# directory, never the manifest file's own location (hvantk/algorithms/cohort/
# frame.py), and a Slurm job's working directory is its log directory, so a bare
# filename or a path relative to $STAGE/cohort/ would not resolve inside a job
# at all. `|` as the sed delimiter, not `/`, since the replacement is a path. Only
# the `table:` line changes -- name, key, key_column, prior and cohort_axes (the
# architecture audit's columns) are carried verbatim.
sed -E "s|^table:.*|table: $STAGE/cohort/gene_burden.tsv|" "$PICK/pos/cohort.yaml" > "$STAGE/cohort/cohort.yaml"
cp "$PICK/pos/receipt.json" "$STAGE/upstream/receipt.json"
[[ -s "$CAB/upstream/gene-index.receipt.json" ]] || {
  echo "ERROR: $CAB/upstream/gene-index.receipt.json not found; run skills/cohort-attach-burden/eval/stage.sh first" >&2
  exit 1
}
cp "$CAB/upstream/gene-index.receipt.json" "$STAGE/upstream/gene-index.receipt.json"

[[ -e "$CAB/gene-index.ht/_SUCCESS" ]] || {
  echo "ERROR: $CAB/gene-index.ht not found (no _SUCCESS marker); run skills/cohort-attach-burden/eval/stage.sh first" >&2
  exit 1
}
# --delete, unconditionally run every time (not skip-if-exists), so a later
# cohort-attach-burden re-stage (a new upstream receipt) is picked up here too --
# this is what "refresh the copies when an upstream receipt changes" means in
# practice: there is nothing to compare against beforehand, so this just always
# re-syncs, which is idempotent and cheap when nothing changed.
rsync -a --delete "$CAB/gene-index.ht/" "$STAGE/gene-index.ht/"

# --- raw downloads, through the image, each plugin's own downloader ---
shopt -s nullglob
clingen_existing=("$STAGE"/raw/clingen/*.csv)
shopt -u nullglob
if (( ${#clingen_existing[@]} > 0 )); then
  echo "already staged: $STAGE/raw/clingen"
else
  run hvantk download clingen --output-dir "$STAGE/raw/clingen"
fi
validate_raw_dir "$STAGE/raw/clingen" "*.csv" clingen

if [[ -s "$STAGE/raw/hgnc/hgnc_complete_set.txt" ]]; then
  echo "already staged: $STAGE/raw/hgnc/hgnc_complete_set.txt"
else
  run hvantk download hgnc --output "$STAGE/raw/hgnc/hgnc_complete_set.txt"
fi
validate_raw_dir "$STAGE/raw/hgnc" "hgnc_complete_set.txt" hgnc

( cd "$STAGE" && sha256sum raw/clingen/*.csv raw/hgnc/hgnc_complete_set.txt \
    cohort/gene_burden.tsv cohort/cohort.yaml upstream/receipt.json \
    upstream/gene-index.receipt.json > staged.sha256 )
cat "$STAGE/staged.sha256"
echo "staged under $STAGE (gene index copied from $CAB/gene-index.ht)"
