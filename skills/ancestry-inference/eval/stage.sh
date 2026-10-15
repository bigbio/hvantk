#!/usr/bin/env bash
# Stage the ancestry eval: download one 1000 Genomes chromosome and its index (login
# node), build the sample table with the onek-genomes plugin (one job: it downloads the
# panel and needs Hail), then prepare the reference, query, truth and shuffled reference
# (one job). Idempotent: skips what exists.
#
#   HVANTK_SITE_PROFILE=... bash skills/ancestry-inference/eval/stage.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"; REL="${HERE#"$R"/}"
SUBMIT="$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh"
PREP="$HVANTK_REPO/$REL/prepare_eval.py"
STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/ancestry-inference"
mkdir -p "$STAGE/raw"

# The chromosome 22 file of the high-coverage phased release (GRCh38, 3202 samples),
# as recorded in docs_site/guide/data-sources.md for the onek-genomes:variants dataset.
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
FILE="CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz"
for f in "$FILE" "$FILE.tbi"; do
  if [[ ! -s "$STAGE/raw/$f" ]]; then
    echo "downloading $f"
    curl -fsSL --retry 3 -C - -o "$STAGE/raw/$f" "$BASE/$f"
  fi
done

if [[ ! -e "$STAGE/samples.ht/_SUCCESS" ]]; then
  echo "building the sample table (one job; downloads the panel)"
  bash "$SUBMIT" --name stage-anc-samples --cpus 2 --mem-gb 8 --time 0-00:20:00 --wait -- \
    hvantk reprocess onek-genomes:samples --raw-dir "$STAGE/raw/panel" --output "$STAGE/samples.ht" --no-check-drift >/dev/null
fi
if [[ ! -e "$STAGE/query.mt/_SUCCESS" || ! -e "$STAGE/reference_shuffled.mt/_SUCCESS" ]]; then
  echo "preparing reference, query, truth and shuffled reference (one job)"
  bash "$SUBMIT" --name stage-anc-prepare --cpus 16 --mem-gb 120 --time 0-01:30:00 --wait -- \
    python "$PREP" --vcf "$STAGE/raw/$FILE" --samples-table "$STAGE/samples.ht" --out "$STAGE" >/dev/null
fi
( cd "$STAGE" && sha256sum "raw/$FILE" truth.tsv > staged.sha256 )
echo "staged under $STAGE: $(wc -l < "$STAGE/truth.tsv") truth rows (incl. header)"
