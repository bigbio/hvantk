#!/usr/bin/env bash
# Positive and negative control for the hpc-container eval. Login-node orchestrator.
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/hpc-container/eval/control.sh
#
# Positive: the skill's own build job into the control directory, accounting captured
#   here, then check.py (through the image, no Hail) must pass.
# Negative: check.py --expect-fail against the site's current default image and its
#   record, which names another commit; M4 must be false (the image is real but not the
#   one asked for).
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
JOB="$(bash "$SUBMIT" --out "$OUT/pos/hvantk.sif" --name "$RUN_ID-build" --wait | tail -1)"
NEW_SHA="$(sha256sum "$OUT/pos/hvantk.sif" | cut -d' ' -f1)"
printf 'commit: %s\ncontainer: %s\njob: %s\n' "$COMMIT" "$NEW_SHA" "$JOB" > "$OUT/pos/run-notes.md"
if [[ "${WRITE_REFERENCE:-0}" == "1" ]]; then
  check --output "$OUT/pos/hvantk.sif" --reference "$REFERENCE" --write-reference --expect-commit "$COMMIT"
fi
[[ -s "$REFERENCE" ]] || { echo "ERROR: $REFERENCE missing; rerun with WRITE_REFERENCE=1" >&2; exit 1; }
sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$OUT/pos/sacct.psv"
check --output "$OUT/pos/hvantk.sif" --run-notes "$OUT/pos/run-notes.md" --reference "$REFERENCE" \
  --run-id "$RUN_ID" --expect-commit "$COMMIT" --expect-sif-sha256 "$NEW_SHA" \
  --site-profile "$HVANTK_SITE_PROFILE" --sacct-file "$OUT/pos/sacct.psv" >/dev/null
for s in M1 M2 M3 M4 M5; do [[ "$(verdict_stage "$OUT/pos/verdict.json" $s)" == "True" ]] || { echo "ERROR: positive control $s is not true" >&2; exit 1; }; done
echo "positive control: PASS (M1..M5 true)"

# --- negative: the default image, built from another commit ---
DEFAULT_REAL="$(readlink -f "$HVANTK_SIF")"
[[ -f "$DEFAULT_REAL.provenance.json" ]] || { echo "ERROR: the default image has no provenance record; cannot run the negative control" >&2; exit 1; }
ln -sfn "$DEFAULT_REAL" "$OUT/neg/hvantk.sif"
ln -sfn "$DEFAULT_REAL.provenance.json" "$OUT/neg/hvantk.sif.provenance.json"
printf 'commit: %s\ncontainer: %s\njob: %s\n' "$COMMIT" "$(sha256sum "$DEFAULT_REAL" | cut -d' ' -f1)" "$JOB" > "$OUT/neg/run-notes.md"
check --output "$OUT/neg/hvantk.sif" --run-notes "$OUT/neg/run-notes.md" --reference "$REFERENCE" \
  --run-id "$RUN_ID" --expect-commit "$COMMIT" --site-profile "$HVANTK_SITE_PROFILE" \
  --sacct-file "$OUT/pos/sacct.psv" --expect-fail >/dev/null
[[ "$(verdict_stage "$OUT/neg/verdict.json" M4)" == "False" ]] || { echo "ERROR: negative control did not fail M4 on the default image's commit" >&2; exit 1; }
echo "negative control: PASS (M4 false: the default image's record names another commit)"
echo "verdicts: $OUT/pos/verdict.json $OUT/neg/verdict.json"
