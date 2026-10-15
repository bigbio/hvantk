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

# The Ensembl GTF's filename carries the pinned release number
# (hvantk/resources/ensembl_release.py), not a calendar date. `hvantk download
# ensembl-structure` already checks for that exact pinned filename and silently
# reuses it when present, so call it unconditionally -- it is a safe no-op on a
# re-run. The point of the assertion AFTER it, rather than an existence check BEFORE
# it, is to fail loudly on a stale GTF from an older pin: a leftover
# Homo_sapiens.GRCh38.<old>.gtf.gz does not satisfy the current downloader's own
# cache check, so it downloads the current pin fresh and leaves two .gtf.gz files
# behind, which the assertion below catches instead of silently building from
# whichever one happens to glob first.
bash "$RUN_IN_IMAGE" -- hvantk download ensembl-structure --raw-dir "$RAW/ensembl-gene"
shopt -s nullglob
ensembl_existing=("$RAW"/ensembl-gene/*.gtf.gz)
shopt -u nullglob
(( ${#ensembl_existing[@]} == 1 )) || { echo "ERROR: $RAW/ensembl-gene must hold exactly one .gtf.gz after staging (a stale GTF from an older Ensembl pin alongside the current one?); found ${#ensembl_existing[@]}: ${ensembl_existing[*]:-(none)}" >&2; exit 1; }

# Unlike the Ensembl downloader, `hvantk download hgnc --output` raises
# (FileExistsError) if its target already exists, even an empty one left behind by an
# interrupted previous run -- `os.path.exists` does not care about size. Remove an
# empty file first so the download call below succeeds; a non-empty file is kept and
# skipped.
if [[ -e "$RAW/hgnc/hgnc_complete_set.txt" && ! -s "$RAW/hgnc/hgnc_complete_set.txt" ]]; then
  echo "removing empty (interrupted?) download: $RAW/hgnc/hgnc_complete_set.txt"
  rm -f "$RAW/hgnc/hgnc_complete_set.txt"
fi
if [[ -s "$RAW/hgnc/hgnc_complete_set.txt" ]]; then
  echo "already staged: $RAW/hgnc/hgnc_complete_set.txt"
else
  bash "$RUN_IN_IMAGE" -- hvantk download hgnc --output "$RAW/hgnc/hgnc_complete_set.txt"
fi

( cd "$STAGE" && sha256sum raw/*/* > staged.sha256 )
cat "$STAGE/staged.sha256"
echo "staged under $STAGE"
