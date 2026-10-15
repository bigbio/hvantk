#!/usr/bin/env bash
# Stage the qtl-cascade eval's inputs, once per cluster. Login-node work (it needs
# egress for the HGNC download and the eQTL Catalogue fetch; it does not need Hail):
# the HGNC complete-set TSV through the plugin's own downloader (skipped once
# staged -- `hvantk download hgnc` exits 1 on an existing file), then one real eQTL
# Catalogue window (GTEx liver, dataset QTD000266) turned into a GTEx v8-layout eQTL
# raw file and a SYNTHETIC Fang-layout pQTL raw file by `stage_qtl.py`. The two QTL
# files are regenerated on every run, after the HGNC step: `stage_qtl.py` writes
# them byte-for-byte deterministically, so re-running costs one quick network fetch
# and never leaves a pQTL file stale against a re-downloaded HGNC file (the pQTL
# file's gene-symbol mapping depends on it). Writes `staged.sha256` next to the raw
# files.
#
#   HVANTK_SITE_PROFILE=... bash skills/qtl-cascade/eval/stage.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"
[[ "$HERE" == "$R"/* ]] || { echo "ERROR: this eval folder ($HERE) is not inside HVANTK_REPO ($R)" >&2; exit 1; }
REL="${HERE#"$R"/}"
RUN_IN_IMAGE="$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh"
STAGE_QTL="$HVANTK_REPO/$REL/stage_qtl.py"

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/qtl-cascade"
RAW="$STAGE/raw"
mkdir -p "$RAW/hgnc" "$RAW/gtex-eqtl" "$RAW/pqtl"

HGNC_FILE="$RAW/hgnc/hgnc_complete_set.txt"
if [[ -s "$HGNC_FILE" ]]; then
  echo "already staged: $HGNC_FILE"
else
  bash "$RUN_IN_IMAGE" -- hvantk download hgnc --output "$HGNC_FILE"
fi

EQTL_FILE="$RAW/gtex-eqtl/Liver.v8.allpairs.txt.gz"
PQTL_FILE="$RAW/pqtl/Liver.allpairs_nobsGE72.txt.gz"
# pysam's remote tabix (what fetch_eqtl_region uses) downloads the remote index
# (a .tbi) into the CURRENT WORKING DIRECTORY, not a temp location -- run from a
# scratch directory under $STAGE (bound into the container the same way $STAGE
# is) so the stray index lands there, not in the checkout, and remove it after.
SCRATCH="$STAGE/tmp"
mkdir -p "$SCRATCH"
( cd "$SCRATCH" && bash "$RUN_IN_IMAGE" -- python "$STAGE_QTL" \
    --hgnc "$HGNC_FILE" --out-eqtl "$EQTL_FILE" --out-pqtl "$PQTL_FILE" )
rm -f "$SCRATCH"/*   # the stray .tbi (and anything else left in the scratch dir)

( cd "$STAGE" && sha256sum raw/*/* > staged.sha256 )
cat "$STAGE/staged.sha256"
echo "staged under $STAGE"
