#!/usr/bin/env bash
# Positive and negative control for the ancestry-inference eval. Login-node
# orchestrator; the pipeline runs as a job through the skill's template; the check runs
# through the image on this node (no Hail needed).
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/ancestry-inference/eval/control.sh
#
# Positive: the pipeline on query.mt against reference.mt; check.py must pass.
# Negative: the same against reference_shuffled.mt (labels permuted); check.py
#   --expect-fail must pass with M4 true (predictions exist) and M5 false (chance-level
#   concordance).
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
JOB="$HVANTK_REPO/${REL%/eval}/templates/ancestry-job.sh"

EVAL_ROOT="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}"
STAGE="$EVAL_ROOT/ancestry-inference"
CTL="$EVAL_ROOT/ancestry-inference-control"
[[ -e "$STAGE/query.mt/_SUCCESS" && -s "$CTL/truth.tsv" && -e "$CTL/reference_shuffled.mt/_SUCCESS" ]] || { echo "ERROR: not staged; run eval/stage.sh first" >&2; exit 1; }
( cd "$EVAL_ROOT" && sha256sum -c --quiet "$STAGE/staged.sha256" ) || { echo "ERROR: staged inputs do not match staged.sha256; run eval/stage.sh again" >&2; exit 1; }
INPUT_SHA="$(sha256sum "$STAGE/staged.sha256" | cut -d' ' -f1)"
OUT="$STAGE/control/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT/pos" "$OUT/neg"
RUN_ID="ctl-anc-$(date -u +%H%M%S)"
SINCE="$(date +%Y-%m-%dT%H:%M:%S)"
COMMIT="$(git -C "$R" rev-parse HEAD)"
SIF_SHA="$(sha256sum "$(readlink -f "$HVANTK_SIF")" | cut -d' ' -f1)"
SACCT_FIELDS="JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"
CONTAINER_BIN="${HVANTK_CONTAINER_BIN:-singularity}"
BINDS=(-B "$HVANTK_WORK:$HVANTK_WORK")
for p in "$HVANTK_REPO" "$STAGE" "$CTL"; do case "$p" in "$HVANTK_WORK"|"$HVANTK_WORK"/*) ;; *) BINDS+=(-B "$p:$p") ;; esac; done
check() { "$CONTAINER_BIN" exec "${BINDS[@]}" --env PYTHONPATH="$HVANTK_REPO" "$HVANTK_SIF" python3 "$CHECK" "$@"; }
ref_field() { python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['produced_by'][sys.argv[2]])" "$REFERENCE" "$1"; }
verdict_stage() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v[sys.argv[2]])" "$1" "$2"; }

echo "control run $RUN_ID -> $OUT (commit $COMMIT)"
if [[ "${WRITE_REFERENCE:-0}" != "1" && -s "$REFERENCE" && "$(ref_field input_sha256)" != "$INPUT_SHA" ]]; then
  echo "ERROR: the staged inputs differ from the ones the reference was produced from; re-stage and rerun with WRITE_REFERENCE=1" >&2; exit 1
fi

# --- positive ---
JOB_POS="$(bash "$JOB" --query "$STAGE/query.mt" --reference "$STAGE/reference.mt" --out "$OUT/pos" \
  --cpus 16 --mem-gb 64 --time 0-01:00:00 --name "$RUN_ID-ancestry" --wait | tail -1)"
printf 'commit: %s\ncontainer: %s\njob: %s\n' "$COMMIT" "$SIF_SHA" "$JOB_POS" > "$OUT/pos/run-notes.md"
if [[ "${WRITE_REFERENCE:-0}" == "1" ]]; then
  check --output "$OUT/pos" --truth "$CTL/truth.tsv" --reference "$REFERENCE" --write-reference \
    --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" --input-sha256 "$INPUT_SHA"
fi
[[ -s "$REFERENCE" ]] || { echo "ERROR: $REFERENCE missing; rerun with WRITE_REFERENCE=1" >&2; exit 1; }
sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$OUT/pos/sacct.psv"
check --output "$OUT/pos" --truth "$CTL/truth.tsv" --run-notes "$OUT/pos/run-notes.md" --reference "$REFERENCE" \
  --run-id "$RUN_ID" --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" \
  --site-profile "$HVANTK_SITE_PROFILE" --sacct-file "$OUT/pos/sacct.psv" >/dev/null || true
for s in M1 M2 M3 M4 M5; do [[ "$(verdict_stage "$OUT/pos/verdict.json" $s)" == "True" ]] || { echo "ERROR: positive control $s is not true" >&2; exit 1; }; done
echo "positive control: PASS (M1..M5 true; concordance $(python3 -c "import json; print(round(json.load(open('$OUT/pos/verdict.json'))['output']['concordance'],4))"))"

# --- negative: permuted labels ---
JOB_NEG="$(bash "$JOB" --query "$STAGE/query.mt" --reference "$CTL/reference_shuffled.mt" --out "$OUT/neg" \
  --cpus 16 --mem-gb 64 --time 0-01:00:00 --name "$RUN_ID-neg" --wait | tail -1)"
printf 'commit: %s\ncontainer: %s\njob: %s\n' "$COMMIT" "$SIF_SHA" "$JOB_NEG" > "$OUT/neg/run-notes.md"
sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$OUT/neg/sacct.psv"
check --output "$OUT/neg" --truth "$CTL/truth.tsv" --run-notes "$OUT/neg/run-notes.md" --reference "$REFERENCE" \
  --run-id "$RUN_ID" --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" \
  --site-profile "$HVANTK_SITE_PROFILE" --sacct-file "$OUT/neg/sacct.psv" --expect-fail >/dev/null || true
[[ "$(verdict_stage "$OUT/neg/verdict.json" M4)" == "True" && "$(verdict_stage "$OUT/neg/verdict.json" M5)" == "False" ]] \
  || { echo "ERROR: negative control did not fail on content alone (M4=$(verdict_stage "$OUT/neg/verdict.json" M4) M5=$(verdict_stage "$OUT/neg/verdict.json" M5))" >&2; exit 1; }
echo "negative control: PASS (M4 true, M5 false: permuted labels give chance-level concordance)"
echo "verdicts: $OUT/pos/verdict.json $OUT/neg/verdict.json"
