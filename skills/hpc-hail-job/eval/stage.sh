#!/usr/bin/env bash
# Stage the public input for the hpc-hail-job eval, once per cluster. Login-node work:
# it downloads, nothing more. Writes staged.sha256 next to the file.
#
#   HVANTK_SITE_PROFILE=... bash skills/hpc-hail-job/eval/stage.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/hpc-hail-job"
RAW="$STAGE/raw/hgnc"
mkdir -p "$RAW"
CONTAINER_BIN="${HVANTK_CONTAINER_BIN:-singularity}"

if [[ -s "$RAW/hgnc_complete_set.txt" ]]; then
  echo "already staged: $RAW/hgnc_complete_set.txt"
else
  # The plugin's own downloader, so the URL has one home (hvantk/skills/hgnc/shared).
  "$CONTAINER_BIN" exec -B "$HVANTK_WORK:$HVANTK_WORK" "$HVANTK_SIF" \
    hvantk download hgnc --output "$RAW/hgnc_complete_set.txt"
fi
( cd "$STAGE" && sha256sum raw/hgnc/hgnc_complete_set.txt > staged.sha256 )
cat "$STAGE/staged.sha256"
echo "staged under $STAGE"
