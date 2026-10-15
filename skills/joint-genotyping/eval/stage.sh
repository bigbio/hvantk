#!/usr/bin/env bash
# Stage the inputs for the joint-genotyping eval: the repository's three HGC test gVCFs
# (single-sample, chr1 slices, GRCh38) copied from the checkout, plus a two-sample
# subset for the negative control. No download. Login-node work.
#
#   HVANTK_SITE_PROFILE=... bash skills/joint-genotyping/eval/stage.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
SRC="$HVANTK_REPO/hvantk/tests/testdata/hgc_data/gvcfs"
[[ -d "$SRC" ]] || { echo "ERROR: test gVCFs not found under $SRC" >&2; exit 1; }
STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/joint-genotyping"
rm -rf "$STAGE/gvcfs" "$STAGE/gvcfs_two"
mkdir -p "$STAGE/gvcfs" "$STAGE/gvcfs_two"
cp "$SRC"/*.g.vcf.gz "$SRC"/*.g.vcf.gz.tbi "$STAGE/gvcfs/"
n=0
for f in "$STAGE"/gvcfs/*.g.vcf.gz; do
  n=$((n + 1))
  (( n <= 2 )) && cp "$f" "$f.tbi" "$STAGE/gvcfs_two/"
done
( cd "$STAGE" && sha256sum gvcfs/*.g.vcf.gz > staged.sha256 )
cat "$STAGE/staged.sha256" | cut -c1-16,66-
echo "staged under $STAGE ($n gVCFs; gvcfs_two holds the first two for the negative control)"
