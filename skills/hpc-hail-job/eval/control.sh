#!/usr/bin/env bash
# Positive and negative control for the hpc-hail-job eval. Login-node orchestrator: it
# only submits jobs through the skill's own submit template and waits on them.
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/hpc-hail-job/eval/control.sh
#
# Positive: build hgnc:lookup from the staged file (the skill's procedure, no agent),
#   then check.py must pass (and with WRITE_REFERENCE=1, first write reference.json).
# Negative: build from the plugin's tiny committed fixture, which has the same schema
#   and a different row count; check.py --expect-fail must pass, proving the check
#   discriminates on content and not only on existence.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL="$(cd "$HERE/.." && pwd)"
SUBMIT="$SKILL/templates/submit-hail-job.sh"
bash "$SKILL/../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/hpc-hail-job"
RAW="$STAGE/raw/hgnc"
[[ -s "$RAW/hgnc_complete_set.txt" ]] || { echo "ERROR: not staged; run eval/stage.sh first" >&2; exit 1; }
INPUT_SHA="$(cut -d' ' -f1 "$STAGE/staged.sha256")"
OUT="$STAGE/control/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT/pos" "$OUT/neg"
RUN_ID="ctl-hhj-$(date -u +%H%M%S)"
COMMIT="$(git -C "$HVANTK_REPO" rev-parse HEAD)"
SIF_SHA="${HVANTK_SIF_SHA256:-$(sha256sum "$(readlink -f "$HVANTK_SIF")" | cut -d' ' -f1)}"
REPO_REL="$(python3 -c "import os,sys; print(os.path.relpath(sys.argv[1], sys.argv[2]))" "$HERE" "$HVANTK_REPO")"
CHECK="$REPO_REL/check.py"
REFERENCE="$HERE/reference.json"
FIXTURE_DIR="hvantk/skills/hgnc/tests/testdata/raw/hgnc"

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

"$SUBMIT" --name "$RUN_ID-check" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
  python "$CHECK" --output "$OUT/pos/hgnc.ht" --run-notes "$OUT/pos/run-notes.md" \
    --reference "$REFERENCE" --run-id "$RUN_ID" --expect-commit "$COMMIT" \
    --expect-sif-sha256 "$SIF_SHA" --site-profile "$HVANTK_SITE_PROFILE" >/dev/null
echo "positive control: PASS ($(cat "$OUT/pos/verdict.json" | python3 -c 'import json,sys; v=json.load(sys.stdin); print({k:v[k] for k in ("M1","M2","M3","M4","M5")})'))"

# --- negative: same schema, wrong content ---
JOB_NEG="$("$SUBMIT" --name "$RUN_ID-neg" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
  hvantk reprocess hgnc:lookup --raw-dir "$FIXTURE_DIR" --output "$OUT/neg/hgnc.ht" --skip-download --no-check-drift | tail -1)"
printf 'commit: %s\ncontainer: %s\njob: %s\n' "$COMMIT" "$SIF_SHA" "$JOB_NEG" > "$OUT/neg/run-notes.md"
"$SUBMIT" --name "$RUN_ID-negcheck" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
  python "$CHECK" --output "$OUT/neg/hgnc.ht" --run-notes "$OUT/neg/run-notes.md" \
    --reference "$REFERENCE" --run-id "$RUN_ID" --expect-commit "$COMMIT" \
    --expect-sif-sha256 "$SIF_SHA" --site-profile "$HVANTK_SITE_PROFILE" --expect-fail >/dev/null
echo "negative control: PASS (the check failed the fixture-built table, as it must)"
echo "verdicts: $OUT/pos/verdict.json $OUT/neg/verdict.json"
