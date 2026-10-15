#!/usr/bin/env bash
# Positive and negative control for the callset-subset-export eval. Login-node
# orchestrator; the chain runs as jobs through the skill's own template, the check runs
# through the image on this node (no Hail needed).
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/callset-subset-export/eval/control.sh
#
# Positive: the two-sample delivery; check.py must pass.
# Negative: a one-sample delivery; check.py --expect-fail must pass with M4 true (a
#   complete delivery) and M5 false (the wrong sample set and count).
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
SUBSET="$HVANTK_REPO/${REL%/eval}/templates/subset-job.sh"

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/callset-subset-export"
[[ -e "$STAGE/cohort/cohort.vds/variant_data/_SUCCESS" && -s "$STAGE/samples.txt" ]] || { echo "ERROR: not staged; run eval/stage.sh first" >&2; exit 1; }
( cd "$STAGE" && sha256sum -c --quiet staged.sha256 ) || { echo "ERROR: staged inputs do not match staged.sha256; run eval/stage.sh again" >&2; exit 1; }
INPUT_SHA="$(sha256sum "$STAGE/staged.sha256" | cut -d' ' -f1)"
OUT="$STAGE/control/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT/pos" "$OUT/neg"
RUN_ID="ctl-cse-$(date -u +%H%M%S)"
SINCE="$(date +%Y-%m-%dT%H:%M:%S)"
COMMIT="$(git -C "$R" rev-parse HEAD)"
SIF_SHA="$(sha256sum "$(readlink -f "$HVANTK_SIF")" | cut -d' ' -f1)"
SACCT_FIELDS="JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"
CONTAINER_BIN="${HVANTK_CONTAINER_BIN:-singularity}"
BINDS=(-B "$HVANTK_WORK:$HVANTK_WORK")
for p in "$HVANTK_REPO" "$STAGE"; do case "$p" in "$HVANTK_WORK"|"$HVANTK_WORK"/*) ;; *) BINDS+=(-B "$p:$p") ;; esac; done
check() { "$CONTAINER_BIN" exec "${BINDS[@]}" --env PYTHONPATH="$HVANTK_REPO" "$HVANTK_SIF" python3 "$CHECK" "$@"; }
ref_field() { python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['produced_by'][sys.argv[2]])" "$REFERENCE" "$1"; }
verdict_stage() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v[sys.argv[2]])" "$1" "$2"; }

echo "control run $RUN_ID -> $OUT (commit $COMMIT)"
if [[ "${WRITE_REFERENCE:-0}" != "1" && -s "$REFERENCE" && "$(ref_field input_sha256)" != "$INPUT_SHA" ]]; then
  echo "ERROR: the staged inputs differ from the ones the reference was produced from; re-stage and rerun with WRITE_REFERENCE=1" >&2; exit 1
fi

# --- positive: two samples ---
mapfile -t JOBS_POS < <(bash "$SUBSET" --vds "$STAGE/cohort/cohort.vds" --samples "$STAGE/samples.txt" --out "$OUT/pos/subset" \
  --prefix subset --cpus 8 --mem-gb 32 --time 0-00:30:00 --name "$RUN_ID-pos")
{ printf 'commit: %s\ncontainer: %s\n' "$COMMIT" "$SIF_SHA"; printf 'job: %s\n' "${JOBS_POS[@]}"; } > "$OUT/pos/run-notes.md"
if [[ "${WRITE_REFERENCE:-0}" == "1" ]]; then
  check --output "$OUT/pos/subset" --reference "$REFERENCE" --write-reference \
    --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" --input-sha256 "$INPUT_SHA"
fi
[[ -s "$REFERENCE" ]] || { echo "ERROR: $REFERENCE missing; rerun with WRITE_REFERENCE=1" >&2; exit 1; }
sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$OUT/pos/sacct.psv"
check --output "$OUT/pos/subset" --run-notes "$OUT/pos/run-notes.md" --reference "$REFERENCE" \
  --run-id "$RUN_ID" --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" \
  --site-profile "$HVANTK_SITE_PROFILE" --sacct-file "$OUT/pos/sacct.psv" >/dev/null || true
for s in M1 M2 M3 M4 M5; do [[ "$(verdict_stage "$OUT/pos/verdict.json" $s)" == "True" ]] || { echo "ERROR: positive control $s is not true" >&2; exit 1; }; done
echo "positive control: PASS (M1..M5 true)"

# --- negative: one sample ---
mapfile -t JOBS_NEG < <(bash "$SUBSET" --vds "$STAGE/cohort/cohort.vds" --samples "$STAGE/samples_one.txt" --out "$OUT/neg/subset" \
  --prefix subset --cpus 8 --mem-gb 32 --time 0-00:30:00 --name "$RUN_ID-neg")
{ printf 'commit: %s\ncontainer: %s\n' "$COMMIT" "$SIF_SHA"; printf 'job: %s\n' "${JOBS_NEG[@]}"; } > "$OUT/neg/run-notes.md"
sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$OUT/neg/sacct.psv"
check --output "$OUT/neg/subset" --run-notes "$OUT/neg/run-notes.md" --reference "$REFERENCE" \
  --run-id "$RUN_ID" --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" \
  --site-profile "$HVANTK_SITE_PROFILE" --sacct-file "$OUT/neg/sacct.psv" --expect-fail >/dev/null || true
[[ "$(verdict_stage "$OUT/neg/verdict.json" M4)" == "True" && "$(verdict_stage "$OUT/neg/verdict.json" M5)" == "False" ]] \
  || { echo "ERROR: negative control did not fail on content alone (M4=$(verdict_stage "$OUT/neg/verdict.json" M4) M5=$(verdict_stage "$OUT/neg/verdict.json" M5))" >&2; exit 1; }
echo "negative control: PASS (M4 true, M5 false: a complete delivery of the wrong samples)"
echo "verdicts: $OUT/pos/verdict.json $OUT/neg/verdict.json"
