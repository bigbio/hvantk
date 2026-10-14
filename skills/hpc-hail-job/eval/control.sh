#!/usr/bin/env bash
# Positive and negative control for the hpc-hail-job eval. Login-node orchestrator: it
# only submits jobs through the skill's own submit template and waits on them.
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/hpc-hail-job/eval/control.sh
#
# Positive: build hgnc:lookup from the staged file (the skill's procedure, no agent),
#   capture Slurm accounting on this node, then check.py must pass (and with
#   WRITE_REFERENCE=1, first write reference.json).
# Negative: build from the plugin's tiny committed fixture, which has the same schema
#   and a different row count; check.py --expect-fail must pass, and the verdict must
#   show M4 true and M5 false, proving the check discriminates on content.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SKILL="$(cd "$HERE/.." && pwd -P)"
SUBMIT="$SKILL/templates/submit-hail-job.sh"
bash "$SKILL/../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"
[[ "$HERE" == "$R"/* ]] || { echo "ERROR: this eval folder ($HERE) is not inside HVANTK_REPO ($R)" >&2; exit 1; }

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/hpc-hail-job"
RAW="$STAGE/raw/hgnc"
[[ -s "$RAW/hgnc_complete_set.txt" ]] || { echo "ERROR: not staged; run eval/stage.sh first" >&2; exit 1; }
INPUT_SHA="$(cut -d' ' -f1 "$STAGE/staged.sha256")"
OUT="$STAGE/control/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT/pos" "$OUT/neg"
RUN_ID="ctl-hhj-$(date -u +%H%M%S)"
SINCE="$(date -u +%Y-%m-%dT%H:%M:%S)"
COMMIT="$(git -C "$R" rev-parse HEAD)"
SIF_SHA="$(sha256sum "$(readlink -f "$HVANTK_SIF")" | cut -d' ' -f1)"   # always the resolved file
CHECK="$HERE/check.py"
REFERENCE="$HERE/reference.json"
FIXTURE_DIR="$R/hvantk/skills/hgnc/tests/testdata/raw/hgnc"
SACCT_FIELDS="JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"

capture_sacct() { sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$1"; }
ref_field() { python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['produced_by'][sys.argv[2]])" "$REFERENCE" "$1"; }
verdict_stage() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v[sys.argv[2]])" "$1" "$2"; }

echo "control run $RUN_ID -> $OUT (commit $COMMIT)"

# --- positive: the skill's procedure, no agent ---
JOB_POS="$("$SUBMIT" --name "$RUN_ID-build" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
  hvantk reprocess hgnc:lookup --raw-dir "$RAW" --output "$OUT/pos/hgnc.ht" --skip-download --no-check-drift | tail -1)"
printf 'commit: %s\ncontainer: %s\njob: %s\n' "$COMMIT" "$SIF_SHA" "$JOB_POS" > "$OUT/pos/run-notes.md"

if [[ "${WRITE_REFERENCE:-0}" == "1" ]]; then
  "$SUBMIT" --name "$RUN_ID-ref" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
    python "$CHECK" --output "$OUT/pos/hgnc.ht" --reference "$REFERENCE" --write-reference \
      --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" --input-sha256 "$INPUT_SHA" >/dev/null
  echo "reference written: $REFERENCE"
fi
[[ -s "$REFERENCE" ]] || { echo "ERROR: $REFERENCE missing; rerun with WRITE_REFERENCE=1" >&2; exit 1; }
REF_INPUT="$(ref_field input_sha256)"
if [[ "$REF_INPUT" != "$INPUT_SHA" ]]; then
  echo "ERROR: the staged input ($INPUT_SHA) differs from the one the reference was produced from ($REF_INPUT)." >&2
  echo "       HGNC changed upstream. Re-stage deliberately and rerun with WRITE_REFERENCE=1 to re-baseline." >&2
  exit 1
fi

capture_sacct "$OUT/pos/sacct.psv"
"$SUBMIT" --name "$RUN_ID-check" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
  python "$CHECK" --output "$OUT/pos/hgnc.ht" --run-notes "$OUT/pos/run-notes.md" \
    --reference "$REFERENCE" --run-id "$RUN_ID" --expect-commit "$COMMIT" \
    --expect-sif-sha256 "$SIF_SHA" --site-profile "$HVANTK_SITE_PROFILE" \
    --sacct-file "$OUT/pos/sacct.psv" >/dev/null
echo "positive control: PASS (M1..M5 = $(for s in M1 M2 M3 M4 M5; do verdict_stage "$OUT/pos/verdict.json" $s; done | tr '\n' ' '))"

# --- negative: same schema, wrong content ---
JOB_NEG="$("$SUBMIT" --name "$RUN_ID-neg" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
  hvantk reprocess hgnc:lookup --raw-dir "$FIXTURE_DIR" --output "$OUT/neg/hgnc.ht" --skip-download --no-check-drift | tail -1)"
printf 'commit: %s\ncontainer: %s\njob: %s\n' "$COMMIT" "$SIF_SHA" "$JOB_NEG" > "$OUT/neg/run-notes.md"
capture_sacct "$OUT/neg/sacct.psv"
"$SUBMIT" --name "$RUN_ID-negcheck" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
  python "$CHECK" --output "$OUT/neg/hgnc.ht" --run-notes "$OUT/neg/run-notes.md" \
    --reference "$REFERENCE" --run-id "$RUN_ID" --expect-commit "$COMMIT" \
    --expect-sif-sha256 "$SIF_SHA" --site-profile "$HVANTK_SITE_PROFILE" \
    --sacct-file "$OUT/neg/sacct.psv" --expect-fail >/dev/null
[[ "$(verdict_stage "$OUT/neg/verdict.json" M4)" == "True" && "$(verdict_stage "$OUT/neg/verdict.json" M5)" == "False" ]] \
  || { echo "ERROR: negative control did not fail on content alone (M4=$(verdict_stage "$OUT/neg/verdict.json" M4) M5=$(verdict_stage "$OUT/neg/verdict.json" M5))" >&2; exit 1; }
echo "negative control: PASS (M4 true, M5 false: the check failed the fixture-built table on content)"
echo "verdicts: $OUT/pos/verdict.json $OUT/neg/verdict.json"
