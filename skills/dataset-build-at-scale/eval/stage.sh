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

# Each plugin's own downloader from the checkout, so the URL has one home
# (hvantk/skills/<provider>/shared) and the checkout's code is what runs. Skip a
# provider whose directory already holds a file: re-staging on another day changes
# ClinGen's and GenCC's dated filename, which would move the reference out from under
# a control that has not been re-run with WRITE_REFERENCE=1.
if find "$RAW/hgnc" -maxdepth 1 -type f ! -name '.*' | grep -q .; then
  echo "already staged: $RAW/hgnc"
else
  run hvantk download hgnc --output "$RAW/hgnc/hgnc_complete_set.txt"
fi
if find "$RAW/clingen" -maxdepth 1 -type f ! -name '.*' | grep -q .; then
  echo "already staged: $RAW/clingen"
else
  run hvantk download clingen --output-dir "$RAW/clingen"
fi
if find "$RAW/gencc" -maxdepth 1 -type f ! -name '.*' | grep -q .; then
  echo "already staged: $RAW/gencc"
else
  run hvantk download gencc --output-dir "$RAW/gencc"
fi

( cd "$STAGE" && sha256sum raw/*/* > staged.sha256 )
cat "$STAGE/staged.sha256"
echo "staged under $STAGE"
