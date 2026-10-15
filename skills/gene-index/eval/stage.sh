#!/usr/bin/env bash
# Stage the public inputs for the gene-index eval, once per cluster. Login-node work,
# through the image: it downloads the two raw inputs `annotate spine` needs and nothing
# more. Writes staged.sha256 next to them. Idempotent: skips what is already staged.
#
#   HVANTK_SITE_PROFILE=... bash skills/gene-index/eval/stage.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
RUN_IN_IMAGE="$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh"

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/gene-index"
RAW="$STAGE/raw"
mkdir -p "$RAW/ensembl-gene" "$RAW/hgnc"

# The Ensembl GTF's filename carries the pinned release number, so match it by
# extension rather than hard-coding a release that will move. The `hgnc:lookup`
# downloader always writes the one fixed filename.
shopt -s nullglob
ensembl_existing=("$RAW"/ensembl-gene/*.gtf.gz)
shopt -u nullglob

if (( ${#ensembl_existing[@]} > 0 )); then
  echo "already staged: ${ensembl_existing[0]}"
else
  bash "$RUN_IN_IMAGE" -- hvantk download ensembl-structure --raw-dir "$RAW/ensembl-gene"
fi
shopt -s nullglob
ensembl_existing=("$RAW"/ensembl-gene/*.gtf.gz)
shopt -u nullglob
(( ${#ensembl_existing[@]} == 1 )) || { echo "ERROR: $RAW/ensembl-gene must hold exactly one .gtf.gz after staging; found ${#ensembl_existing[@]}" >&2; exit 1; }

if [[ -s "$RAW/hgnc/hgnc_complete_set.txt" ]]; then
  echo "already staged: $RAW/hgnc/hgnc_complete_set.txt"
else
  bash "$RUN_IN_IMAGE" -- hvantk download hgnc --output "$RAW/hgnc/hgnc_complete_set.txt"
fi

( cd "$STAGE" && sha256sum raw/*/* > staged.sha256 )
cat "$STAGE/staged.sha256"
echo "staged under $STAGE"
