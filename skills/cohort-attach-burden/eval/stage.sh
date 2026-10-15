#!/usr/bin/env bash
# Stage the cohort-attach-burden eval: the newest PASSING gene-index positive control's
# own table, the 1000 Genomes high-coverage chr22 VCF (reusing ancestry-inference's
# staged copy when it is already complete, rather than downloading a multi-GB file
# twice on one cluster), and the synthetic case/control cohort MatrixTable built from
# both. Login-node work except for the one Hail prep job. The gene index is always
# re-synced from the newest passing upstream control (cheap); the cohort MatrixTable
# is only rebuilt when that upstream changed, was incomplete, or is itself missing.
#
#   HVANTK_SITE_PROFILE=... bash skills/cohort-attach-burden/eval/stage.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"
[[ "$HERE" == "$R"/* ]] || { echo "ERROR: this eval folder ($HERE) is not inside HVANTK_REPO ($R)" >&2; exit 1; }
REL="${HERE#"$R"/}"
SUBMIT="$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh"
PREP="$HVANTK_REPO/$REL/prepare_cohort.py"

EVAL_ROOT="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}"
STAGE="$EVAL_ROOT/cohort-attach-burden"
mkdir -p "$STAGE/raw" "$STAGE/upstream"

# --- 1. the newest gene-index positive control that actually passed ---
GI_STAGE="$EVAL_ROOT/gene-index"
GI_POS="$(python3 -c "
import json, os, sys

root = sys.argv[1]
ctl = os.path.join(root, 'control')
best = None
if os.path.isdir(ctl):
    for name in sorted(os.listdir(ctl)):
        pos = os.path.join(ctl, name, 'pos')
        verdict_path = os.path.join(pos, 'verdict.json')
        if not os.path.isfile(verdict_path):
            continue
        try:
            with open(verdict_path) as fh:
                verdict = json.load(fh)
        except ValueError:
            continue
        if verdict.get('pass') is True:
            best = pos
if best is None:
    sys.exit(1)
sys.stdout.write(best)
" "$GI_STAGE")" || {
  echo "ERROR: no passing gene-index positive control found under $GI_STAGE/control; run the gene-index control first" >&2
  exit 1
}
echo "using gene-index control: $GI_POS"
# Detect a changed (re-baselined) or previously incomplete upstream BEFORE
# overwriting the stage's own copy of the receipt, so the comparison is against what
# was actually staged last time, not against what we are about to replace it with.
GI_WAS_INCOMPLETE=1
[[ -e "$STAGE/gene-index.ht/_SUCCESS" ]] && GI_WAS_INCOMPLETE=0
OLD_RECEIPT_SHA=""
[[ -f "$STAGE/upstream/gene-index.receipt.json" ]] \
  && OLD_RECEIPT_SHA="$(sha256sum "$STAGE/upstream/gene-index.receipt.json" | cut -d' ' -f1)"
NEW_RECEIPT_SHA="$(sha256sum "$GI_POS/receipt.json" | cut -d' ' -f1)"

# Always sync, never skip: --delete removes any part file or stale leftover from an
# earlier, different upstream control so the staged copy is always exactly $GI_POS's.
rsync -a --delete "$GI_POS/gene-index.ht/" "$STAGE/gene-index.ht/"
cp "$GI_POS/receipt.json" "$STAGE/upstream/gene-index.receipt.json"

if (( GI_WAS_INCOMPLETE )) || [[ "$OLD_RECEIPT_SHA" != "$NEW_RECEIPT_SHA" ]]; then
  echo "upstream gene index changed (or the staged copy was incomplete); removing the stale cohort MatrixTable so it rebuilds"
  rm -rf "$STAGE/cohort.mt" "$STAGE/cohort.summary.json"
fi

# --- 2. the 1000 Genomes high-coverage chr22 VCF and its index (the same file
#     ancestry-inference's eval stages; reuse its copy first) ---
ANC_RAW="$EVAL_ROOT/ancestry-inference/raw"
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
FILE="CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz"

fetch() {  # fetch <url> <target>: resumable, atomic, size checked against the server's
           # (same helper as ancestry-inference/eval/stage.sh); reuses
           # ancestry-inference's own staged copy first when it is already complete.
  local want have
  want="$(curl -sIL "$1" | tr -d '\r' | awk 'tolower($1)=="content-length:"{n=$2} END{print n}')"
  if [[ -s "$2" ]]; then
    have="$(stat -c %s "$2" 2>/dev/null || stat -f %z "$2")"
    [[ -z "$want" || "$have" == "$want" ]] && return 0
    echo "incomplete $(basename "$2") ($have of $want bytes); resuming"
    mv "$2" "$2.part"
  fi
  if [[ -s "$ANC_RAW/$(basename "$2")" ]]; then
    have="$(stat -c %s "$ANC_RAW/$(basename "$2")" 2>/dev/null || stat -f %z "$ANC_RAW/$(basename "$2")")"
    if [[ -n "$want" && "$have" == "$want" ]]; then
      echo "reusing $(basename "$2") from ancestry-inference's staged copy"
      cp "$ANC_RAW/$(basename "$2")" "$2"
      return 0
    fi
  fi
  echo "downloading $(basename "$2")"
  curl -fsSL --retry 3 -C - -o "$2.part" "$1"
  mv "$2.part" "$2"
}
fetch "$BASE/$FILE" "$STAGE/raw/$FILE"
fetch "$BASE/$FILE.tbi" "$STAGE/raw/$FILE.tbi"

# --- 3. the synthetic cohort MatrixTable (one job, no download inside it) ---
if [[ ! -e "$STAGE/cohort.mt/_SUCCESS" || ! -s "$STAGE/cohort.summary.json" ]]; then
  echo "preparing the synthetic case/control cohort (one job, up to ~1h)"
  bash "$SUBMIT" --name stage-cab-cohort --cpus 16 --mem-gb 64 --time 0-01:00:00 --wait -- \
    python "$PREP" --vcf "$STAGE/raw/$FILE" --gene-index "$STAGE/gene-index.ht" \
      --out "$STAGE/cohort.mt" --summary "$STAGE/cohort.summary.json" >/dev/null
fi

( cd "$STAGE" && sha256sum "raw/$FILE" "raw/$FILE.tbi" cohort.summary.json \
    upstream/gene-index.receipt.json > staged.sha256 )
cat "$STAGE/staged.sha256"
echo "staged under $STAGE"
