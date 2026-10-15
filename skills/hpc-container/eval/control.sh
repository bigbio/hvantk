#!/usr/bin/env bash
# Positive and negative control for the hpc-container eval. Login-node orchestrator.
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/hpc-container/eval/control.sh
#
# Positive: the skill's own build job into the control directory, accounting captured
#   here, then check.py (through the image, no Hail) must pass.
# Negative: check.py --expect-fail against the image just built, with a copy of its
#   record naming another commit; M4 must be false while the image and record are
#   present, so the failure is about content, not about a missing file.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SKILL="$(cd "$HERE/.." && pwd -P)"
bash "$SKILL/../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"
[[ "$HERE" == "$R"/* ]] || { echo "ERROR: this eval folder ($HERE) is not inside HVANTK_REPO ($R)" >&2; exit 1; }
REL="${HERE#"$R"/}"
CHECK="$HVANTK_REPO/$REL/check.py"
REFERENCE="$HVANTK_REPO/$REL/reference.json"
SUBMIT="$HVANTK_REPO/${REL%/eval}/templates/submit-build-image.sh"

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/hpc-container"
[[ -f "$STAGE/staged.sha256" ]] || { echo "ERROR: run eval/stage.sh first" >&2; exit 1; }
OUT="$STAGE/control/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT/pos" "$OUT/neg"
RUN_ID="ctl-hc-$(date -u +%H%M%S)"
SINCE="$(date +%Y-%m-%dT%H:%M:%S)"
COMMIT="$(git -C "$R" rev-parse HEAD)"
SACCT_FIELDS="JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"
CONTAINER_BIN="${HVANTK_CONTAINER_BIN:-singularity}"
BINDS=(-B "$HVANTK_WORK:$HVANTK_WORK")
for p in "$HVANTK_REPO" "$STAGE"; do case "$p" in "$HVANTK_WORK"|"$HVANTK_WORK"/*) ;; *) BINDS+=(-B "$p:$p") ;; esac; done
check() { "$CONTAINER_BIN" exec "${BINDS[@]}" --env PYTHONPATH="$HVANTK_REPO" "$HVANTK_SIF" python3 "$CHECK" "$@"; }
verdict_stage() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v[sys.argv[2]])" "$1" "$2"; }

echo "control run $RUN_ID -> $OUT (commit $COMMIT)"

# --- positive: the skill's own build ---
JOB="$(bash "$SUBMIT" --out "$OUT/pos/hvantk-eval.sif" --name "$RUN_ID-build" --wait | tail -1)"
NEW_SHA="$(sha256sum "$OUT/pos/hvantk-eval.sif" | cut -d' ' -f1)"
printf 'commit: %s\ncontainer: %s\njob: %s\n' "$COMMIT" "$NEW_SHA" "$JOB" > "$OUT/pos/run-notes.md"
if [[ "${WRITE_REFERENCE:-0}" == "1" ]]; then
  check --output "$OUT/pos/hvantk-eval.sif" --reference "$REFERENCE" --write-reference --expect-commit "$COMMIT"
fi
[[ -s "$REFERENCE" ]] || { echo "ERROR: $REFERENCE missing; rerun with WRITE_REFERENCE=1" >&2; exit 1; }
sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$OUT/pos/sacct.psv"
check --output "$OUT/pos/hvantk-eval.sif" --run-notes "$OUT/pos/run-notes.md" --reference "$REFERENCE" \
  --run-id "$RUN_ID" --expect-commit "$COMMIT" --expect-sif-sha256 "$NEW_SHA" \
  --site-profile "$HVANTK_SITE_PROFILE" --sacct-file "$OUT/pos/sacct.psv" >/dev/null || true
for s in M1 M2 M3 M4 M5; do [[ "$(verdict_stage "$OUT/pos/verdict.json" $s)" == "True" ]] || { echo "ERROR: positive control $s is not true" >&2; exit 1; }; done
echo "positive control: PASS (M1..M5 true)"

# --- negative: the same image with a record that names another commit ---
ln -sfn "$OUT/pos/hvantk-eval.sif" "$OUT/neg/hvantk-eval.sif"
python3 - "$OUT/pos/hvantk-eval.sif.provenance.json" "$OUT/neg/hvantk-eval.sif.provenance.json" <<'PY'
import json, sys
rec = json.load(open(sys.argv[1]))
rec["git_commit"] = "0" * 40
json.dump(rec, open(sys.argv[2], "w"), indent=2)
PY
printf 'commit: %s\ncontainer: %s\njob: %s\n' "$COMMIT" "$NEW_SHA" "$JOB" > "$OUT/neg/run-notes.md"
check --output "$OUT/neg/hvantk-eval.sif" --run-notes "$OUT/neg/run-notes.md" --reference "$REFERENCE" \
  --run-id "$RUN_ID" --expect-commit "$COMMIT" --site-profile "$HVANTK_SITE_PROFILE" \
  --sacct-file "$OUT/pos/sacct.psv" --expect-fail >/dev/null || true
NEG_IMG="$(python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v['output']['image'], v['output']['record'])" "$OUT/neg/verdict.json")"
[[ "$(verdict_stage "$OUT/neg/verdict.json" M4)" == "False" && "$NEG_IMG" == "True True" ]] \
  || { echo "ERROR: negative control did not fail M4 on content (M4=$(verdict_stage "$OUT/neg/verdict.json" M4), image/record present: $NEG_IMG)" >&2; exit 1; }
echo "negative control: PASS (M4 false on a present image whose record names another commit)"
echo "verdicts: $OUT/pos/verdict.json $OUT/neg/verdict.json"
