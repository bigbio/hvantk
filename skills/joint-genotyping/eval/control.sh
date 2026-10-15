#!/usr/bin/env bash
# Positive and negative control for the joint-genotyping eval. Login-node orchestrator;
# every heavy step is a job through the skill's own pipeline template.
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/joint-genotyping/eval/control.sh
#
# Positive: the pipeline on the three staged gVCFs, accounting captured here, check.py
#   (in the image, on a compute node) must pass.
# Negative: the pipeline on the two-sample subset; check.py --expect-fail must pass with
#   M4 true (a valid callset) and M5 false (the wrong sample set and counts).
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
PIPE="$HVANTK_REPO/${REL%/eval}/templates/pipeline-job.sh"
SUBMIT="$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh"

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/joint-genotyping"
[[ -d "$STAGE/gvcfs" && -s "$STAGE/staged.sha256" ]] || { echo "ERROR: not staged; run eval/stage.sh first" >&2; exit 1; }
INPUT_SHA="$(sha256sum "$STAGE/staged.sha256" | cut -d' ' -f1)"
OUT="$STAGE/control/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT/pos" "$OUT/neg"
RUN_ID="ctl-jg-$(date -u +%H%M%S)"
SINCE="$(date +%Y-%m-%dT%H:%M:%S)"
COMMIT="$(git -C "$R" rev-parse HEAD)"
SIF_SHA="$(sha256sum "$(readlink -f "$HVANTK_SIF")" | cut -d' ' -f1)"
SACCT_FIELDS="JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"
ref_field() { python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['produced_by'][sys.argv[2]])" "$REFERENCE" "$1"; }
verdict_stage() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v[sys.argv[2]])" "$1" "$2"; }
check_job() { bash "$SUBMIT" --name "$1" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- python "$CHECK" "${@:2}" >/dev/null || true; }

echo "control run $RUN_ID -> $OUT (commit $COMMIT)"
( cd "$STAGE" && sha256sum -c --quiet staged.sha256 ) || { echo "ERROR: staged gVCFs do not match staged.sha256; run eval/stage.sh again" >&2; exit 1; }
if [[ "${WRITE_REFERENCE:-0}" != "1" && -s "$REFERENCE" && "$(ref_field input_sha256)" != "$INPUT_SHA" ]]; then
  echo "ERROR: the staged inputs differ from the ones the reference was produced from; re-stage and rerun with WRITE_REFERENCE=1" >&2; exit 1
fi

# --- positive ---
mapfile -t JOBS_POS < <(bash "$PIPE" --gvcf-dir "$STAGE/gvcfs" --out "$OUT/pos/cohort" --prefix cohort \
  --cpus 8 --mem-gb 32 --time 0-00:30:00 --name "$RUN_ID-hgc")
{ printf 'commit: %s\ncontainer: %s\n' "$COMMIT" "$SIF_SHA"; printf 'job: %s\n' "${JOBS_POS[@]}"; } > "$OUT/pos/run-notes.md"
if [[ "${WRITE_REFERENCE:-0}" == "1" ]]; then
  check_job "$RUN_ID-ref" --output "$OUT/pos/cohort" --reference "$REFERENCE" --write-reference \
    --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" --input-sha256 "$INPUT_SHA"
  echo "reference written: $REFERENCE"
fi
[[ -s "$REFERENCE" ]] || { echo "ERROR: $REFERENCE missing; rerun with WRITE_REFERENCE=1" >&2; exit 1; }
sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$OUT/pos/sacct.psv"
check_job "$RUN_ID-check" --output "$OUT/pos/cohort" --run-notes "$OUT/pos/run-notes.md" --reference "$REFERENCE" \
  --run-id "$RUN_ID" --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" \
  --site-profile "$HVANTK_SITE_PROFILE" --sacct-file "$OUT/pos/sacct.psv"
for s in M1 M2 M3 M4 M5; do [[ "$(verdict_stage "$OUT/pos/verdict.json" $s)" == "True" ]] || { echo "ERROR: positive control $s is not true" >&2; exit 1; }; done
echo "positive control: PASS (M1..M5 true)"

# --- negative: two of the three samples ---
mapfile -t JOBS_NEG < <(bash "$PIPE" --gvcf-dir "$STAGE/gvcfs_two" --out "$OUT/neg/cohort" --prefix cohort \
  --cpus 8 --mem-gb 32 --time 0-00:30:00 --name "$RUN_ID-neg")
{ printf 'commit: %s\ncontainer: %s\n' "$COMMIT" "$SIF_SHA"; printf 'job: %s\n' "${JOBS_NEG[@]}"; } > "$OUT/neg/run-notes.md"
sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$OUT/neg/sacct.psv"
check_job "$RUN_ID-negcheck" --output "$OUT/neg/cohort" --run-notes "$OUT/neg/run-notes.md" --reference "$REFERENCE" \
  --run-id "$RUN_ID" --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" \
  --site-profile "$HVANTK_SITE_PROFILE" --sacct-file "$OUT/neg/sacct.psv" --expect-fail
[[ "$(verdict_stage "$OUT/neg/verdict.json" M4)" == "True" && "$(verdict_stage "$OUT/neg/verdict.json" M5)" == "False" ]] \
  || { echo "ERROR: negative control did not fail on content alone (M4=$(verdict_stage "$OUT/neg/verdict.json" M4) M5=$(verdict_stage "$OUT/neg/verdict.json" M5))" >&2; exit 1; }
echo "negative control: PASS (M4 true, M5 false: a valid callset with the wrong samples)"
echo "verdicts: $OUT/pos/verdict.json $OUT/neg/verdict.json"
