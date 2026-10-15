#!/usr/bin/env bash
# Stage the ancestry eval: download one 1000 Genomes chromosome and its index and the
# sample panel (login node), build the sample table with the onek-genomes plugin (one
# job, no download inside it), then prepare the reference and query (one job). The
# held-out truth and the shuffled reference go to a sibling directory the task never
# names. Idempotent: skips what exists; downloads go to .part files first.
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
EVAL_ROOT="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}"
STAGE="$EVAL_ROOT/ancestry-inference"
CTL="$EVAL_ROOT/ancestry-inference-control"
mkdir -p "$STAGE/raw/panel" "$CTL"
rm -rf "$STAGE/truth.tsv" "$STAGE/reference_shuffled.mt"   # never leave the answer key next to the inputs

fetch() {  # fetch <url> <target>: resumable, atomic, and the size is checked against the server's
  local want have
  want="$(curl -sIL "$1" | tr -d '\r' | awk 'tolower($1)=="content-length:"{n=$2} END{print n}')"
  if [[ -s "$2" ]]; then
    have="$(stat -c %s "$2" 2>/dev/null || stat -f %z "$2")"
    [[ -z "$want" || "$have" == "$want" ]] && return 0
    echo "incomplete $(basename "$2") ($have of $want bytes); resuming"
    mv "$2" "$2.part"
  fi
  echo "downloading $(basename "$2")"
  curl -fsSL --retry 3 -C - -o "$2.part" "$1"
  mv "$2.part" "$2"
}
# Chromosome 22 of the high-coverage phased release (GRCh38, 3202 samples), as recorded
# in docs_site/guide/data-sources.md for the onek-genomes:variants dataset; and the
# 2504-sample panel the onek-genomes:samples plugin downloads (same public URL).
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
FILE="CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz"
PANEL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
fetch "$BASE/$FILE" "$STAGE/raw/$FILE"
fetch "$BASE/$FILE.tbi" "$STAGE/raw/$FILE.tbi"
if [[ ! -s "$STAGE/raw/panel/igsr_samples.tsv" ]]; then
  fetch "$PANEL" "$STAGE/raw/panel/igsr_samples.panel"
  sed -E 's/[[:space:]]+$//' "$STAGE/raw/panel/igsr_samples.panel" > "$STAGE/raw/panel/igsr_samples.tsv"
  rm -f "$STAGE/raw/panel/igsr_samples.panel"
fi

if [[ ! -e "$STAGE/samples.ht/_SUCCESS" ]]; then
  echo "building the sample table (one job, from the downloaded panel)"
  bash "$SUBMIT" --name stage-anc-samples --cpus 2 --mem-gb 8 --time 0-00:20:00 --wait -- \
    hvantk reprocess onek-genomes:samples --raw-dir "$STAGE/raw/panel" --output "$STAGE/samples.ht" --skip-download --no-check-drift >/dev/null
fi
if [[ ! -e "$STAGE/query.mt/_SUCCESS" || ! -e "$CTL/reference_shuffled.mt/_SUCCESS" ]]; then
  echo "preparing reference and query, truth and shuffled reference (one job)"
  bash "$SUBMIT" --name stage-anc-prepare --cpus 16 --mem-gb 64 --time 0-02:00:00 --wait -- \
    python "$PREP" --vcf "$STAGE/raw/$FILE" --samples-table "$STAGE/samples.ht" --out "$STAGE" --control-out "$CTL" >/dev/null
fi
( cd "$EVAL_ROOT" && sha256sum "ancestry-inference/raw/$FILE" "ancestry-inference/raw/panel/igsr_samples.tsv" "ancestry-inference-control/truth.tsv" > "$STAGE/staged.sha256" )
echo "staged under $STAGE (truth and shuffled reference under $CTL): $(($(wc -l < "$CTL/truth.tsv") - 1)) held-out samples"
