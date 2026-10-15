#!/usr/bin/env bash
# Positive and negative control for the dataset-build-at-scale eval. Login-node
# orchestrator: it only submits jobs through the skill's own templates and waits on
# them.
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/dataset-build-at-scale/eval/control.sh
#
# Positive: the three datasets built as ONE array job through submit-build-array.sh
#   (the skill's own procedure, no agent); check.py must pass, proving M1..M5 true and
#   array_used true (with WRITE_REFERENCE=1, first write reference.json).
# Negative: the same three builds as three SEPARATE submit-hail-job.sh jobs, no array;
#   check.py --expect-fail must pass, with M4 true and M5 false because content_ok is
#   true (same raw files, same deterministic build) but array_used is false -- the
#   check must fail on "not an array", not on wrong content.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SKILL="$(cd "$HERE/.." && pwd -P)"
SUBMIT_ARRAY="$SKILL/templates/submit-build-array.sh"
SUBMIT="$SKILL/../hpc-hail-job/templates/submit-hail-job.sh"
[[ -f "$SUBMIT_ARRAY" ]] || { echo "ERROR: submit-build-array.sh not found: $SUBMIT_ARRAY" >&2; exit 2; }
[[ -f "$SUBMIT" ]] || { echo "ERROR: hpc-hail-job's submit-hail-job.sh not found: $SUBMIT" >&2; exit 2; }
bash "$SKILL/../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"
[[ "$HERE" == "$R"/* ]] || { echo "ERROR: this eval folder ($HERE) is not inside HVANTK_REPO ($R)" >&2; exit 1; }
# Paths handed to a job (anything run with `python ...` inside the container) use the
# profile's spelling of the checkout, which is what the job binds; the resolved $R
# above is only for the containment test.
REL="${HERE#"$R"/}"
CHECK="$HVANTK_REPO/$REL/check.py"
REFERENCE="$HVANTK_REPO/$REL/reference.json"

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/dataset-build-at-scale"
RAW="$STAGE/raw"
[[ -s "$STAGE/staged.sha256" ]] || { echo "ERROR: not staged; run eval/stage.sh first" >&2; exit 1; }
HGNC_RAW="$RAW/hgnc"
CLINGEN_RAW="$RAW/clingen"
GENCC_RAW="$RAW/gencc"
INPUT_SHA="$(sha256sum "$STAGE/staged.sha256" | cut -d' ' -f1)"
OUT="$STAGE/control/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT/pos" "$OUT/neg"
BASE="ctl-dbs-$(date -u +%H%M%S)"
RUN_POS="$BASE-p"
RUN_NEG="$BASE-n"
SINCE="$(date +%Y-%m-%dT%H:%M:%S)"   # local time: sacct -S reads the node's clock, not UTC
COMMIT="$(git -C "$R" rev-parse HEAD)"
SIF_SHA="$(sha256sum "$(readlink -f "$HVANTK_SIF")" | cut -d' ' -f1)"   # always the resolved file
SACCT_FIELDS="JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"

capture_sacct() { sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$1"; }
# `produced_by` is null in the placeholder reference.json ({"produced_by": null})
# the implementer commits before any control has run; `or {}` reads that the same as
# a reference with no matching field, rather than crashing on None has no .get.
ref_field() { python3 -c "import json,sys
d=json.load(open(sys.argv[1]))
p=d.get('produced_by') or {}
print(p.get(sys.argv[2], ''))" "$REFERENCE" "$1"; }
verdict_stage() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v.get(sys.argv[2]))" "$1" "$2"; }
print_verdict() { python3 -c "import json,sys
v=json.load(open(sys.argv[1]))
print('verdict:', {k: v.get(k) for k in ('M1','M2','M3','M4','M5','pass','content_ok','array_used')})" "$1"; }

# Login-node orchestrator wrapper for check.py: `|| true` so a job that legitimately
# exits non-zero (the negative control under --expect-fail is SUPPOSED to; a genuine
# check bug would too) never aborts this script under set -e before the verdict has
# been printed and the assertions below have had a chance to say exactly what failed.
run_check() {
  local name="$1" verdict_dir="$2"; shift 2
  bash "$SUBMIT" --name "$name" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
    python "$CHECK" "$@" >/dev/null || true
  print_verdict "$verdict_dir/verdict.json"
}

echo "control run $BASE -> $OUT (commit $COMMIT)"
( cd "$STAGE" && sha256sum -c --quiet staged.sha256 ) || { echo "ERROR: staged inputs do not match staged.sha256; delete them and run eval/stage.sh again" >&2; exit 1; }
# sha256sum -c only checks that the files IT lists still match; it says nothing about
# an extra file dropped into raw/ since staging, or a provider's raw dir now missing
# its file. Compare the set of paths, not just their hashes.
EXPECTED_PATHS="$(cut -d' ' -f3- "$STAGE/staged.sha256" | sort)"
ACTUAL_PATHS="$(cd "$STAGE" && printf '%s\n' raw/*/* | sort)"
[[ "$EXPECTED_PATHS" == "$ACTUAL_PATHS" ]] || { echo "ERROR: raw/ no longer holds exactly the files staged.sha256 lists; delete the stale raw dirs and run eval/stage.sh again" >&2; echo "  expected: $EXPECTED_PATHS" >&2; echo "  actual:   $ACTUAL_PATHS" >&2; exit 1; }
if [[ "${WRITE_REFERENCE:-0}" != "1" && -s "$REFERENCE" ]]; then
  REF_INPUT="$(ref_field input_sha256)"
  if [[ -n "$REF_INPUT" && "$REF_INPUT" != "$INPUT_SHA" ]]; then
    echo "ERROR: the staged inputs ($INPUT_SHA) differ from the ones the reference was produced from ($REF_INPUT)." >&2
    echo "       One of HGNC/ClinGen/GenCC changed upstream, or re-staging moved ClinGen/GenCC's dated filename." >&2
    echo "       Re-stage (eval/stage.sh) and rerun with WRITE_REFERENCE=1 to re-baseline." >&2
    exit 1
  fi
fi

# --- positive: the three datasets as one array job (the skill's own procedure) ---
TASKS_POS="$OUT/pos/tasks.tsv"
{
  printf 'hgnc:lookup\t%s\t%s\n' "$HGNC_RAW" "$OUT/pos/hgnc-lookup.ht"
  printf 'clingen:gene-disease\t%s\t%s\n' "$CLINGEN_RAW" "$OUT/pos/clingen-gene-disease.ht"
  printf 'gencc:submissions\t%s\t%s\n' "$GENCC_RAW" "$OUT/pos/gencc-submissions.ht"
} > "$TASKS_POS"

JOB_POS="$(bash "$SUBMIT_ARRAY" --list "$TASKS_POS" --name "$RUN_POS-build" --cpus 4 --mem-gb 16 \
  --time 0-00:30:00 --max-concurrent 3 --wait | tail -1)"
# The array's own job id is the one and only `job:` line: every task is accounted
# under it (<id>_<task>), which is exactly what check.py's array-aware M3 expects.
printf 'commit: %s\ncontainer: %s\njob: %s\n' "$COMMIT" "$SIF_SHA" "$JOB_POS" > "$OUT/pos/run-notes.md"

if [[ "${WRITE_REFERENCE:-0}" == "1" ]]; then
  bash "$SUBMIT" --name "$RUN_POS-ref" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
    python "$CHECK" --out-dir "$OUT/pos" --reference "$REFERENCE" --write-reference \
      --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" --input-sha256 "$INPUT_SHA"
  echo "reference written: $REFERENCE"
fi
[[ -s "$REFERENCE" ]] || { echo "ERROR: $REFERENCE missing; rerun with WRITE_REFERENCE=1" >&2; exit 1; }

capture_sacct "$OUT/pos/sacct.psv"
run_check "$RUN_POS-check" "$OUT/pos" --out-dir "$OUT/pos" --run-notes "$OUT/pos/run-notes.md" \
  --reference "$REFERENCE" --run-id "$RUN_POS" --expect-commit "$COMMIT" \
  --expect-sif-sha256 "$SIF_SHA" --site-profile "$HVANTK_SITE_PROFILE" \
  --sacct-file "$OUT/pos/sacct.psv"
for s in M1 M2 M3 M4 M5; do [[ "$(verdict_stage "$OUT/pos/verdict.json" $s)" == "True" ]] || { echo "ERROR: positive control $s is not true" >&2; exit 1; }; done
[[ "$(verdict_stage "$OUT/pos/verdict.json" array_used)" == "True" ]] || { echo "ERROR: positive control did not record array_used=true" >&2; exit 1; }
echo "positive control: PASS (M1..M5 true, array_used true)"

# --- negative: the same three builds as three SEPARATE single-job submissions ---
JOB_NEG_HGNC="$(bash "$SUBMIT" --name "$RUN_NEG-hgnc" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
  hvantk reprocess hgnc:lookup --raw-dir "$HGNC_RAW" --output "$OUT/neg/hgnc-lookup.ht" \
    --skip-download --no-check-drift | tail -1)"
JOB_NEG_CLINGEN="$(bash "$SUBMIT" --name "$RUN_NEG-clingen" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
  hvantk reprocess clingen:gene-disease --raw-dir "$CLINGEN_RAW" --output "$OUT/neg/clingen-gene-disease.ht" \
    --skip-download --no-check-drift | tail -1)"
JOB_NEG_GENCC="$(bash "$SUBMIT" --name "$RUN_NEG-gencc" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
  hvantk reprocess gencc:submissions --raw-dir "$GENCC_RAW" --output "$OUT/neg/gencc-submissions.ht" \
    --skip-download --no-check-drift | tail -1)"
{
  printf 'commit: %s\ncontainer: %s\n' "$COMMIT" "$SIF_SHA"
  printf 'job: %s\njob: %s\njob: %s\n' "$JOB_NEG_HGNC" "$JOB_NEG_CLINGEN" "$JOB_NEG_GENCC"
} > "$OUT/neg/run-notes.md"

capture_sacct "$OUT/neg/sacct.psv"
run_check "$RUN_NEG-check" "$OUT/neg" --out-dir "$OUT/neg" --run-notes "$OUT/neg/run-notes.md" \
  --reference "$REFERENCE" --run-id "$RUN_NEG" --expect-commit "$COMMIT" \
  --expect-sif-sha256 "$SIF_SHA" --site-profile "$HVANTK_SITE_PROFILE" \
  --sacct-file "$OUT/neg/sacct.psv" --expect-fail
[[ "$(verdict_stage "$OUT/neg/verdict.json" M4)" == "True" && "$(verdict_stage "$OUT/neg/verdict.json" M5)" == "False" ]] \
  || { echo "ERROR: negative control did not fail the expected way (M4=$(verdict_stage "$OUT/neg/verdict.json" M4) M5=$(verdict_stage "$OUT/neg/verdict.json" M5))" >&2; exit 1; }
[[ "$(verdict_stage "$OUT/neg/verdict.json" content_ok)" == "True" && "$(verdict_stage "$OUT/neg/verdict.json" array_used)" == "False" ]] \
  || { echo "ERROR: negative control did not fail specifically on array_used (content_ok=$(verdict_stage "$OUT/neg/verdict.json" content_ok) array_used=$(verdict_stage "$OUT/neg/verdict.json" array_used))" >&2; exit 1; }
echo "negative control: PASS (M4 true, M5 false: content_ok true but array_used false -- three separate jobs, not an array)"
echo "verdicts: $OUT/pos/verdict.json $OUT/neg/verdict.json"
