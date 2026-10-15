#!/usr/bin/env bash
# Stage the callset-subset-export eval: a three-sample callset built from the
# repository's HGC test gVCFs by the joint-genotyping skill's own template (one job),
# the two-sample list in the callset's spelling, and a one-sample list for the negative
# control. Login-node orchestrator.
#
#   HVANTK_SITE_PROFILE=... bash skills/callset-subset-export/eval/stage.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"; REL="${HERE#"$R"/}"
PIPE="$HVANTK_REPO/skills/joint-genotyping/templates/pipeline-job.sh"
SRC="$HVANTK_REPO/hvantk/tests/testdata/hgc_data/gvcfs"
STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/callset-subset-export"
mkdir -p "$STAGE/gvcfs"
cp -n "$SRC"/*.g.vcf.gz "$SRC"/*.g.vcf.gz.tbi "$STAGE/gvcfs/" 2>/dev/null || true

# Sample names from the gVCF headers, in sorted order; the first two are the subset.
: > "$STAGE/samples_all.txt"
for f in "$STAGE"/gvcfs/*.g.vcf.gz; do
  { zcat "$f" || true; } | grep -m1 '^#CHROM' | cut -f10- >> "$STAGE/samples_all.txt"
done
sort -o "$STAGE/samples_all.txt" "$STAGE/samples_all.txt"
head -2 "$STAGE/samples_all.txt" > "$STAGE/samples.txt"
head -1 "$STAGE/samples_all.txt" > "$STAGE/samples_one.txt"

if [[ ! -d "$STAGE/cohort/cohort.vds" ]]; then
  echo "building the callset (four jobs through the joint-genotyping template)"
  bash "$PIPE" --gvcf-dir "$STAGE/gvcfs" --out "$STAGE/cohort" --prefix cohort \
    --cpus 8 --mem-gb 32 --time 0-00:30:00 --name stage-cse-hgc >/dev/null
fi
[[ -e "$STAGE/cohort/cohort.vds/variant_data/_SUCCESS" ]] || { echo "ERROR: the staged callset is incomplete" >&2; exit 1; }
( cd "$STAGE" && sha256sum gvcfs/*.g.vcf.gz samples.txt > staged.sha256 )
echo "samples: $(tr '\n' ' ' < "$STAGE/samples.txt")"
echo "staged under $STAGE"
