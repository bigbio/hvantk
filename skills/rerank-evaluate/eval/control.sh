#!/usr/bin/env bash
# Positive and negative control for the rerank-evaluate eval. Login-node orchestrator:
# it only submits jobs through hpc-hail-job's own submit-hail-job.sh and runs check.py
# as a direct container exec (no Hail, no Slurm job needed for the check itself).
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/rerank-evaluate/eval/control.sh
#
# Positive: the skill's own procedure end to end -- clingen build, geneset extraction,
#   cohort and feature-axis symbol mapping (features RESTRICTED to the cohort's tested
#   genes), labels flattened and restricted to the same set, rerank --blocks --seed 42,
#   evaluate; check.py must pass (and with WRITE_REFERENCE=1, first write
#   reference.json).
# Negative: identical in every step EXCEPT the feature axis is built from the WHOLE
#   gene index, not restricted to the cohort's tested genes (labels stay restricted
#   the same way as the positive); check.py --expect-fail must pass, with M4 true and
#   M5 false because universe_ok alone is false -- the check must fail on "wrong
#   universe", not on a broken run.
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
# Paths handed to a job (anything run with `python ...`/`hvantk ...` inside the
# container) use the profile's spelling of the checkout, which is what the job binds;
# the resolved $R above is only for the containment test.
REL="${HERE#"$R"/}"
CHECK="$HVANTK_REPO/$REL/check.py"
REFERENCE="$HVANTK_REPO/$REL/reference.json"
PREP_ABS="$HVANTK_REPO/${REL%/eval}/templates/prepare-rerank-inputs.py"
EVAL_ABS="$HVANTK_REPO/${REL%/eval}/templates/evaluate-ranking.py"

EVAL_ROOT="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}"
STAGE="$EVAL_ROOT/rerank-evaluate"
[[ -s "$STAGE/staged.sha256" ]] || { echo "ERROR: not staged; run eval/stage.sh first" >&2; exit 1; }
( cd "$STAGE" && sha256sum -c --quiet staged.sha256 ) || { echo "ERROR: staged inputs do not match staged.sha256; delete them and run eval/stage.sh again" >&2; exit 1; }
INPUT_SHA="$(sha256sum "$STAGE/staged.sha256" | cut -d' ' -f1)"

OUT="$STAGE/control/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT/pos" "$OUT/neg"
BASE="ctl-rev-$(date -u +%H%M%S)"
RUN_POS="$BASE-p"
RUN_NEG="$BASE-n"
SINCE="$(date +%Y-%m-%dT%H:%M:%S)"   # local time: sacct -S reads the node's clock, not UTC
COMMIT="$(git -C "$R" rev-parse HEAD)"
SIF_SHA="$(sha256sum "$(readlink -f "$HVANTK_SIF")" | cut -d' ' -f1)"   # always the resolved file
CONTAINER_BIN="${HVANTK_CONTAINER_BIN:-singularity}"
BINDS=(-B "$HVANTK_WORK:$HVANTK_WORK")
for p in "$HVANTK_REPO" "$STAGE"; do case "$p" in "$HVANTK_WORK"|"$HVANTK_WORK"/*) ;; *) BINDS+=(-B "$p:$p") ;; esac; done
SACCT_FIELDS="JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"

capture_sacct() { sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$1"; }
check() { "$CONTAINER_BIN" exec "${BINDS[@]}" --env PYTHONPATH="$HVANTK_REPO" "$HVANTK_SIF" python3 "$CHECK" "$@"; }
# `produced_by` is null in the placeholder reference.json ({"produced_by": null}) the
# implementer commits before any control has run; `or {}` reads that the same as a
# reference with no matching field, rather than crashing on None has no .get.
ref_field() { python3 -c "import json,sys
d=json.load(open(sys.argv[1]))
p=d.get('produced_by') or {}
print(p.get(sys.argv[2], ''))" "$REFERENCE" "$1"; }
verdict_stage() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v.get(sys.argv[2]))" "$1" "$2"; }
content_check() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v.get('content_checks',{}).get(sys.argv[2]))" "$1" "$2"; }
print_verdict() { python3 -c "import json,sys
v=json.load(open(sys.argv[1]))
print('verdict:', {k: v.get(k) for k in ('M1','M2','M3','M4','M5','pass')})
print('content_checks:', v.get('content_checks'))" "$1"; }

# Login-node orchestrator wrapper for check.py: `|| true` so a check that legitimately
# exits non-zero (the negative control under --expect-fail is SUPPOSED to; a genuine
# check bug would too) never aborts this script under set -e before the assertions
# below have had a chance to say exactly what failed.
run_check() { check "$@" >/dev/null || true; }

echo "control run $BASE -> $OUT (commit $COMMIT)"
if [[ "${WRITE_REFERENCE:-0}" != "1" && -s "$REFERENCE" ]]; then
  REF_INPUT="$(ref_field input_sha256)"
  if [[ -n "$REF_INPUT" && "$REF_INPUT" != "$INPUT_SHA" ]]; then
    echo "ERROR: the staged inputs ($INPUT_SHA) differ from the ones the reference was produced from ($REF_INPUT)." >&2
    echo "       Re-stage (eval/stage.sh) and rerun with WRITE_REFERENCE=1 to re-baseline." >&2
    exit 1
  fi
fi

# --- one arm (pos or neg): clingen build, geneset extraction, cohort/feature/label
# --- symbol mapping, rerank, evaluate. $3 switches whether the feature axis is
# --- restricted to the cohort's tested genes (1, the positive) or not (0, the
# --- negative) -- the one difference the two arms exist to exercise.
run_arm() {
  local arm="$1" run_prefix="$2" restrict_features="$3"
  local dir="$OUT/$arm"
  local jobs=() j

  j="$(bash "$SUBMIT" --name "$run_prefix-clingen" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
    hvantk reprocess clingen:gene-disease --raw-dir "$STAGE/raw/clingen" \
      --output "$dir/clingen.ht" --skip-download --no-check-drift | tail -1)"
  jobs+=("$j")

  j="$(bash "$SUBMIT" --name "$run_prefix-geneset" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
    hvantk genesets clingen --ht "$dir/clingen.ht" --min-classification Moderate \
      -o "$dir/labels.json" | tail -1)"
  jobs+=("$j")

  j="$(bash "$SUBMIT" --name "$run_prefix-cohort" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
    python "$PREP_ABS" cohort --table "$STAGE/cohort/gene_burden.tsv" --key-column gene \
      --gene-index "$STAGE/gene-index.ht" --hgnc-tsv "$STAGE/raw/hgnc/hgnc_complete_set.txt" \
      --out "$dir/cohort_symbols.tsv" --manifest-in "$STAGE/cohort/cohort.yaml" \
      --manifest-out "$dir/cohort_symbols.yaml" | tail -1)"
  jobs+=("$j")

  j="$(bash "$SUBMIT" --name "$run_prefix-labels" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
    python "$PREP_ABS" labels --geneset "$dir/labels.json" \
      --restrict-to "$dir/cohort_symbols.tsv" --restrict-column symbol \
      --out "$dir/labels.txt" | tail -1)"
  jobs+=("$j")

  if [[ "$restrict_features" == "1" ]]; then
    j="$(bash "$SUBMIT" --name "$run_prefix-features" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
      python "$PREP_ABS" features --ht "$STAGE/gene-index.ht" \
        --fields cds_length,n_coding_exons,n_transcripts \
        --hgnc-tsv "$STAGE/raw/hgnc/hgnc_complete_set.txt" \
        --restrict-to "$dir/cohort_symbols.tsv" --restrict-column symbol \
        --out "$dir/structure.tsv" | tail -1)"
  else
    j="$(bash "$SUBMIT" --name "$run_prefix-features" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
      python "$PREP_ABS" features --ht "$STAGE/gene-index.ht" \
        --fields cds_length,n_coding_exons,n_transcripts \
        --hgnc-tsv "$STAGE/raw/hgnc/hgnc_complete_set.txt" \
        --out "$dir/structure.tsv" | tail -1)"
  fi
  jobs+=("$j")

  cat > "$dir/config.yaml" <<YAML
name: rerank-evaluate-$arm
cohort: $dir/cohort_symbols.yaml
features:
  - name: structure
    path: $dir/structure.tsv
labels:
  path: $dir/labels.txt
YAML

  j="$(bash "$SUBMIT" --name "$run_prefix-rerank" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
    hvantk rerank -c "$dir/config.yaml" -o "$dir/ranked.tsv" \
      --blocks "$STAGE/raw/hgnc/hgnc_complete_set.txt" --seed 42 | tail -1)"
  jobs+=("$j")

  j="$(bash "$SUBMIT" --name "$run_prefix-eval" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
    python "$EVAL_ABS" --ranked "$dir/ranked.tsv" --receipt "$dir/eval.json" | tail -1)"
  jobs+=("$j")

  local auc
  auc="$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['auc'])" "$dir/eval.json")"
  {
    printf 'commit: %s\n' "$COMMIT"
    printf 'container: %s\n' "$SIF_SHA"
    for j in "${jobs[@]}"; do printf 'job: %s\n' "$j"; done
    printf 'auc: %s\n' "$auc"
  } > "$dir/run-notes.md"
}

# --- positive: features restricted to the cohort's tested genes ---
echo "=== positive arm ($RUN_POS) ==="
run_arm pos "$RUN_POS" 1

if [[ "${WRITE_REFERENCE:-0}" == "1" ]]; then
  check --config "$OUT/pos/config.yaml" --ranked "$OUT/pos/ranked.tsv" --reference "$REFERENCE" \
    --write-reference --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" --input-sha256 "$INPUT_SHA"
  echo "reference written: $REFERENCE"
fi
[[ -s "$REFERENCE" ]] || { echo "ERROR: $REFERENCE missing; rerun with WRITE_REFERENCE=1" >&2; exit 1; }

capture_sacct "$OUT/pos/sacct.psv"
run_check --config "$OUT/pos/config.yaml" --ranked "$OUT/pos/ranked.tsv" --run-notes "$OUT/pos/run-notes.md" \
  --reference "$REFERENCE" --run-id "$RUN_POS" --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" \
  --site-profile "$HVANTK_SITE_PROFILE" --sacct-file "$OUT/pos/sacct.psv"
for s in M1 M2 M3 M4 M5; do
  [[ "$(verdict_stage "$OUT/pos/verdict.json" $s)" == "True" ]] \
    || { echo "ERROR: positive control $s is not true" >&2; print_verdict "$OUT/pos/verdict.json"; exit 1; }
done
echo "positive control: PASS"
print_verdict "$OUT/pos/verdict.json"

# --- negative: the feature axis built from the WHOLE gene index, not the cohort's
# --- tested genes; labels restricted the same way as the positive ---
echo "=== negative arm ($RUN_NEG): features not restricted to the cohort's genes ==="
run_arm neg "$RUN_NEG" 0

capture_sacct "$OUT/neg/sacct.psv"
run_check --config "$OUT/neg/config.yaml" --ranked "$OUT/neg/ranked.tsv" --run-notes "$OUT/neg/run-notes.md" \
  --reference "$REFERENCE" --run-id "$RUN_NEG" --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" \
  --site-profile "$HVANTK_SITE_PROFILE" --sacct-file "$OUT/neg/sacct.psv" --expect-fail
[[ "$(verdict_stage "$OUT/neg/verdict.json" M4)" == "True" && "$(verdict_stage "$OUT/neg/verdict.json" M5)" == "False" ]] \
  || { echo "ERROR: negative control did not fail the expected way (M4=$(verdict_stage "$OUT/neg/verdict.json" M4) M5=$(verdict_stage "$OUT/neg/verdict.json" M5))" >&2; exit 1; }
[[ "$(content_check "$OUT/neg/verdict.json" universe_ok)" == "False" ]] \
  || { echo "ERROR: negative control did not fail specifically on universe_ok ($(content_check "$OUT/neg/verdict.json" universe_ok))" >&2; exit 1; }
echo "negative control: PASS (M4 true, M5 false: universe_ok false -- the whole gene index was ranked, not the cohort's tested genes)"
print_verdict "$OUT/neg/verdict.json"
echo "verdicts: $OUT/pos/verdict.json $OUT/neg/verdict.json"
