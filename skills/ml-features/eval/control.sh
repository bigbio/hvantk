#!/usr/bin/env bash
# Positive and negative control for the ml-features eval. Login-node orchestrator: it
# only submits jobs through hpc-hail-job's submit-hail-job.sh and waits on them.
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/ml-features/eval/control.sh
#
# Positive: the skill's own procedure, no agent -- the seven input datasets each built
#   as a job from the staged raw files, one `annotate prepare` job per layer1 axis, one
#   `annotate compose` job, the `training-set build` job, then the skill's receipt
#   template as a job. check.py must pass, proving M1..M5 true (with WRITE_REFERENCE=1,
#   first write reference.json from this build).
# Negative: the training set rebuilt from spec_negative.yaml (min_review_stars: 0) on the
#   SAME built inputs and gene matrix -- a real, self-consistent build with the wrong
#   label filter. Graded with --run-id "$RUN_ID-n" so M2/M3 see only its own job.
#   check.py --expect-fail must pass with M4 true and M5 false because content_ok is
#   false while reported_ok is true: the check must fail on "wrong content", not on
#   "notes don't match the dataset".
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SKILL="$(cd "$HERE/.." && pwd -P)"
SUBMIT="$SKILL/../hpc-hail-job/templates/submit-hail-job.sh"
[[ -f "$SUBMIT" ]] || { echo "ERROR: hpc-hail-job's submit-hail-job.sh not found: $SUBMIT" >&2; exit 2; }
bash "$SKILL/../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"
[[ "$HERE" == "$R"/* ]] || { echo "ERROR: this eval folder ($HERE) is not inside HVANTK_REPO ($R)" >&2; exit 1; }
# Paths handed to a job use the profile's spelling of the checkout, which is what the
# job binds; the resolved $R above is only for the containment test.
REL="${HERE#"$R"/}"
CHECK="$HVANTK_REPO/$REL/check.py"
REFERENCE="$HVANTK_REPO/$REL/reference.json"
DESCRIBE="$HVANTK_REPO/${REL%/eval}/templates/describe-training-set.py"
NEG_SPEC="$HVANTK_REPO/$REL/spec_negative.yaml"

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/ml-features"
RAW="$STAGE/raw"
[[ -s "$STAGE/staged.sha256" ]] || { echo "ERROR: not staged; run eval/stage.sh first" >&2; exit 1; }
( cd "$STAGE" && sha256sum -c --quiet staged.sha256 ) || { echo "ERROR: staged inputs do not match staged.sha256; delete them and run eval/stage.sh again" >&2; exit 1; }
EXPECTED_PATHS="$(cut -d' ' -f3- "$STAGE/staged.sha256" | grep '^raw/' | sort)"
ACTUAL_PATHS="$(cd "$STAGE" && printf '%s\n' raw/*/* | sort)"
[[ "$EXPECTED_PATHS" == "$ACTUAL_PATHS" ]] || { echo "ERROR: raw/ no longer holds exactly the files staged.sha256 lists; re-stage" >&2; exit 1; }
INPUT_SHA="$(sha256sum "$STAGE/staged.sha256" | cut -d' ' -f1)"

OUT="$STAGE/control/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT/pos/inputs" "$OUT/neg/inputs"
RUN_ID="ctl-mlf-$(date -u +%H%M%S)"
RUN_POS="$RUN_ID-p"
RUN_NEG="$RUN_ID-n"
SINCE="$(date +%Y-%m-%dT%H:%M:%S)"   # local time: sacct -S reads the node's clock, not UTC
COMMIT="$(git -C "$R" rev-parse HEAD)"
SIF_SHA="$(sha256sum "$(readlink -f "$HVANTK_SIF")" | cut -d' ' -f1)"   # always the resolved file
SACCT_FIELDS="JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"

capture_sacct() { sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$1"; }
ref_field() { python3 -c "import json,sys
d=json.load(open(sys.argv[1]))
p=d.get('produced_by') or {}
print(p.get(sys.argv[2], ''))" "$REFERENCE" "$1"; }
verdict_stage() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v.get(sys.argv[2]))" "$1" "$2"; }
print_verdict() { python3 -c "import json,sys
v=json.load(open(sys.argv[1]))
print('verdict:', {k: v.get(k) for k in ('M1','M2','M3','M4','M5','pass','content_ok','reported_ok','error')})" "$1"; }
# rows and positives from a build's manifest, as "ROWS POSITIVES"
manifest_counts() { python3 -c "import json,sys
c=json.load(open(sys.argv[1]))['counts']
print(c['rows'], c['label'].get('1', 0))" "$1"; }

# One job sized to the task budget. Under `set -e` a failing submission would abort
# on the command substitution with a bare non-zero exit; print the captured stderr
# (the submit script's request, accounting and log path) first so it is diagnosable.
submit_build() {
  local name="$1" stderr_file="$2" job_id rc
  shift 2
  set +e
  job_id="$(bash "$SUBMIT" --name "$name" --cpus 4 --mem-gb 16 --time 0-00:45:00 --wait -- "$@" 2> "$stderr_file" | tail -1)"
  rc=$?
  set -e
  if (( rc != 0 )); then
    echo "ERROR: job $name did not complete; captured stderr ($stderr_file):" >&2
    cat "$stderr_file" >&2
    exit 1
  fi
  echo "$job_id"
}

# add_job NAME STDERR CMD...: submit_build, then append the job id to JOBS. Two steps,
# because a failing command substitution aborts under `set -e` only in a plain
# assignment, not reliably inside an array append.
JOBS=()
add_job() {
  local j
  j="$(submit_build "$@")"
  JOBS+=("$j")
}

echo "control run $RUN_ID -> $OUT (commit $COMMIT)"
if [[ "${WRITE_REFERENCE:-0}" != "1" && -s "$REFERENCE" ]]; then
  REF_INPUT="$(ref_field input_sha256)"
  if [[ -n "$REF_INPUT" && "$REF_INPUT" != "$INPUT_SHA" ]]; then
    echo "ERROR: the staged inputs ($INPUT_SHA) differ from the ones the reference was produced from ($REF_INPUT)." >&2
    echo "       UniProt PTM is a live snapshot with no archived releases, so re-staging on a later day changes it." >&2
    echo "       Re-stage (eval/stage.sh) and rerun with WRITE_REFERENCE=1 to re-baseline." >&2
    exit 1
  fi
fi

# --- positive: the seven input datasets, one job each ---
I="$OUT/pos/inputs"
build_input() {  # build_input SHORT DATASET RAW_SUBDIR OUTPUT_NAME
  add_job "$RUN_POS-$1" "$OUT/pos/$1.stderr" \
    hvantk reprocess "$2" --raw-dir "$RAW/$3" --output "$I/$4" --skip-download --no-check-drift
}
build_input clinvar clinvar:variants clinvar clinvar.ht
build_input dbnsfp dbnsfp:variants dbnsfp dbnsfp.ht
build_input insider insider:variants insider_variants insider_variants.ht
build_input ptm uniprot-ptm:sites uniprot_ptm ptm_sites.ht
build_input interfaces insider:interfaces insider_interfaces insider_interfaces.ht
build_input gnomad gnomad-metrics:metrics gnomad_metrics gnomad_metrics.ht
build_input gevir gevir:metrics gevir gevir.ht

# --- the gene matrix: one prepare job per layer1 axis, then compose ---
SPEC="$STAGE/spec.yaml"
GENE_INDEX="$STAGE/upstream/gene-index.ht"
HGNC="$STAGE/upstream/hgnc-lookup.ht"
mkdir -p "$OUT/pos/prepared"
prepare_axis() {  # prepare_axis AXIS INPUT [--hgnc PATH]
  local axis="$1" input="$2"; shift 2
  add_job "$RUN_POS-prep-$axis" "$OUT/pos/prep-$axis.stderr" \
    hvantk annotate prepare --spec "$SPEC" --axis "$axis" --input "$input" \
      --spine "$GENE_INDEX" --output "$OUT/pos/prepared/$axis.ht" "$@"
}
prepare_axis gnomad_constraint "$I/gnomad_metrics.ht"
prepare_axis gevir "$I/gevir.ht"
prepare_axis insider_interfaces "$I/insider_interfaces.ht" --hgnc "$HGNC"
add_job "$RUN_POS-compose" "$OUT/pos/compose.stderr" \
  hvantk annotate compose --spec "$SPEC" --spine "$GENE_INDEX" \
    --prepared gnomad_constraint="$OUT/pos/prepared/gnomad_constraint.ht" \
    --prepared gevir="$OUT/pos/prepared/gevir.ht" \
    --prepared insider_interfaces="$OUT/pos/prepared/insider_interfaces.ht" \
    --output "$OUT/pos/gene_matrix.ht"

# --- the training set, then the skill's receipt ---
ts_command() {  # ts_command SPEC OUTDIR: sets TS to the training-set build command
  TS=(hvantk training-set build --spec "$1"
    --input clinvar="$I/clinvar.ht" --input dbnsfp="$I/dbnsfp.ht"
    --input insider_variants="$I/insider_variants.ht" --input ptm="$I/ptm_sites.ht"
    --input hgnc="$HGNC" --gene-matrix "$OUT/pos/gene_matrix.ht"
    --gene-list is_eval_gene="$STAGE/genes.txt" --output "$2/clinvar_chr22")
}
ts_command "$SPEC" "$OUT/pos"
add_job "$RUN_POS-build" "$OUT/pos/build.stderr" "${TS[@]}"
add_job "$RUN_POS-describe" "$OUT/pos/describe.stderr" \
  python "$DESCRIBE" --prefix "$OUT/pos/clinvar_chr22" \
    --receipt "$OUT/pos/clinvar_chr22.receipt.json"

read -r POS_ROWS POS_POSITIVES <<< "$(manifest_counts "$OUT/pos/clinvar_chr22.manifest.json")"
{
  printf 'commit: %s\ncontainer: %s\n' "$COMMIT" "$SIF_SHA"
  printf 'job: %s\n' "${JOBS[@]}"
  printf 'rows: %s\npositives: %s\n' "$POS_ROWS" "$POS_POSITIVES"
} > "$OUT/pos/run-notes.md"

if [[ "${WRITE_REFERENCE:-0}" == "1" ]]; then
  bash "$SUBMIT" --name "$RUN_POS-ref" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
    python "$CHECK" --out-dir "$OUT/pos" --reference "$REFERENCE" --write-reference \
      --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" --input-sha256 "$INPUT_SHA"
  echo "reference written: $REFERENCE"
fi
[[ -s "$REFERENCE" ]] || { echo "ERROR: $REFERENCE missing; rerun with WRITE_REFERENCE=1" >&2; exit 1; }

capture_sacct "$OUT/pos/sacct.psv"
bash "$SUBMIT" --name "$RUN_POS-check" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
  python "$CHECK" --out-dir "$OUT/pos" --run-notes "$OUT/pos/run-notes.md" \
    --reference "$REFERENCE" --run-id "$RUN_POS" --expect-commit "$COMMIT" \
    --expect-sif-sha256 "$SIF_SHA" --site-profile "$HVANTK_SITE_PROFILE" \
    --sacct-file "$OUT/pos/sacct.psv" >/dev/null || true
print_verdict "$OUT/pos/verdict.json"
for s in M1 M2 M3 M4 M5; do [[ "$(verdict_stage "$OUT/pos/verdict.json" $s)" == "True" ]] || { echo "ERROR: positive control $s is not true" >&2; exit 1; }; done
echo "positive control: PASS (M1..M5 true)"

# --- negative: min_review_stars 0, on the SAME built inputs and gene matrix ---
# check.py's M4 requires --out-dir to hold the input builds' sidecars (so an agent
# cannot pass by copying a finished dataset); the negative reads the positive's inputs,
# so copy their sidecars into neg/inputs/ to keep M4 true and fail on content_ok alone.
cp "$I"/*.provenance.json "$OUT/neg/inputs/"
ts_command "$NEG_SPEC" "$OUT/neg"
JOB_NEG="$(submit_build "$RUN_NEG-build" "$OUT/neg/build.stderr" "${TS[@]}")"
read -r NEG_ROWS NEG_POSITIVES <<< "$(manifest_counts "$OUT/neg/clinvar_chr22.manifest.json")"
printf 'commit: %s\ncontainer: %s\njob: %s\nrows: %s\npositives: %s\n' \
  "$COMMIT" "$SIF_SHA" "$JOB_NEG" "$NEG_ROWS" "$NEG_POSITIVES" > "$OUT/neg/run-notes.md"

capture_sacct "$OUT/neg/sacct.psv"
bash "$SUBMIT" --name "$RUN_NEG-check" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
  python "$CHECK" --out-dir "$OUT/neg" --run-notes "$OUT/neg/run-notes.md" \
    --reference "$REFERENCE" --run-id "$RUN_NEG" --expect-commit "$COMMIT" \
    --expect-sif-sha256 "$SIF_SHA" --site-profile "$HVANTK_SITE_PROFILE" \
    --sacct-file "$OUT/neg/sacct.psv" --expect-fail >/dev/null || true
print_verdict "$OUT/neg/verdict.json"
[[ "$(verdict_stage "$OUT/neg/verdict.json" M4)" == "True" && "$(verdict_stage "$OUT/neg/verdict.json" M5)" == "False" ]] \
  || { echo "ERROR: negative control did not fail the expected way (M4=$(verdict_stage "$OUT/neg/verdict.json" M4) M5=$(verdict_stage "$OUT/neg/verdict.json" M5))" >&2; exit 1; }
[[ "$(verdict_stage "$OUT/neg/verdict.json" reported_ok)" == "True" && "$(verdict_stage "$OUT/neg/verdict.json" content_ok)" == "False" ]] \
  || { echo "ERROR: negative control did not fail specifically on content_ok" >&2; exit 1; }
echo "negative control: PASS (M4 true, M5 false: reported_ok true but content_ok false)"
echo "verdicts: $OUT/pos/verdict.json $OUT/neg/verdict.json"

# --- very last action: the marker a downstream stage.sh selects on ---
TSV_SHA="$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['tsv_sha256'])" "$OUT/pos/receipt.json")"
python3 -c "import json,sys; json.dump({'training_set_tsv_sha256': sys.argv[2]}, open(sys.argv[1], 'w'), indent=2, sort_keys=True)" \
  "$OUT/control.passed" "$TSV_SHA"
echo "control.passed: $OUT/control.passed"
