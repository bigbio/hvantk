#!/usr/bin/env bash
# Positive and negative control for the gene-index eval. Login-node orchestrator: it
# only submits jobs through the skill's own templates and waits on them.
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/gene-index/eval/control.sh
#
# Positive: the skill's own procedure, no agent -- ensembl-gene:structure and
#   hgnc:lookup each built as their own job, then the protein-coding gene index (the
#   command's default --biotype) as a third job; the control reads the gene-index job's
#   own two printed lines (`spine: N genes`, `HGNC mapping rate: X%`) out of its Slurm
#   log, whose path submit-hail-job.sh --wait prints on its own stderr as `log: ...`.
#   check.py must pass, proving M1..M5 true (with WRITE_REFERENCE=1, first write
#   reference.json from this build).
# Negative: the gene index rebuilt with --biotype all, from the SAME two already-built
#   inputs (no rebuild of ensembl-gene:structure or hgnc:lookup) -- a real, internally
#   self-consistent build, just the wrong biotype filter. Graded with --run-id
#   "$RUN_ID-n" so M2/M3 see only this one job. check.py --expect-fail must pass, with
#   M4 true (the schema contract still holds) and M5 false because content_ok is false
#   (it does not reproduce the reference row count or digest) while reported_ok is true
#   (its own run-notes genes:/mapping_rate: lines agree with what it actually built) --
#   the check must fail on "wrong content", not on "notes don't match the table".
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
# Paths handed to a job (anything run with `hvantk`/`python ...` inside the container)
# use the profile's spelling of the checkout, which is what the job binds; the resolved
# $R above is only for the containment test.
REL="${HERE#"$R"/}"
CHECK="$HVANTK_REPO/$REL/check.py"
REFERENCE="$HVANTK_REPO/$REL/reference.json"

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/gene-index"
RAW="$STAGE/raw"
[[ -s "$STAGE/staged.sha256" ]] || { echo "ERROR: not staged; run eval/stage.sh first" >&2; exit 1; }
( cd "$STAGE" && sha256sum -c --quiet staged.sha256 ) || { echo "ERROR: staged inputs do not match staged.sha256; delete them and run eval/stage.sh again" >&2; exit 1; }
# sha256sum -c only checks that the files IT lists still match; it says nothing about an
# extra file dropped into raw/ since staging. Compare the set of paths too.
EXPECTED_PATHS="$(cut -d' ' -f3- "$STAGE/staged.sha256" | sort)"
ACTUAL_PATHS="$(cd "$STAGE" && printf '%s\n' raw/*/* | sort)"
[[ "$EXPECTED_PATHS" == "$ACTUAL_PATHS" ]] || { echo "ERROR: raw/ no longer holds exactly the files staged.sha256 lists; delete the stale raw dirs and run eval/stage.sh again" >&2; echo "  expected: $EXPECTED_PATHS" >&2; echo "  actual:   $ACTUAL_PATHS" >&2; exit 1; }
INPUT_SHA="$(sha256sum "$STAGE/staged.sha256" | cut -d' ' -f1)"

OUT="$STAGE/control/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT/pos" "$OUT/neg"
RUN_ID="ctl-gi-$(date -u +%H%M%S)"
RUN_POS="$RUN_ID-p"
RUN_NEG="$RUN_ID-n"
SINCE="$(date +%Y-%m-%dT%H:%M:%S)"   # local time: sacct -S reads the node's clock, not UTC
COMMIT="$(git -C "$R" rev-parse HEAD)"
SIF_SHA="$(sha256sum "$(readlink -f "$HVANTK_SIF")" | cut -d' ' -f1)"   # always the resolved file
SACCT_FIELDS="JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"

capture_sacct() { sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$1"; }
# `produced_by` is null in the placeholder reference.json ({"produced_by": null}) the
# implementer commits before any control has run; `or {}` reads that the same as a
# reference with no matching field, rather than crashing on None has no .get.
ref_field() { python3 -c "import json,sys
d=json.load(open(sys.argv[1]))
p=d.get('produced_by') or {}
print(p.get(sys.argv[2], ''))" "$REFERENCE" "$1"; }
verdict_stage() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v.get(sys.argv[2]))" "$1" "$2"; }
print_verdict() { python3 -c "import json,sys
v=json.load(open(sys.argv[1]))
print('verdict:', {k: v.get(k) for k in ('M1','M2','M3','M4','M5','pass','content_ok','reported_ok')})" "$1"; }

# Parses the gene-index job's own two printed lines (`spine: N genes`, `HGNC mapping
# rate: X%`) out of its Slurm log. The log's PATH is on submit-hail-job.sh's own
# stderr (captured to $1 by the caller), never its stdout (which carries only the job
# id); the two printed lines themselves are inside that log file, because Slurm's
# --output captures the submitted command's own stdout. Echoes "N_GENES MAPPING_PCT".
parse_gi_log() {
  local stderr_file="$1" log_path genes_line rate_line n_genes mapping_pct
  log_path="$(grep -m1 '^log: ' "$stderr_file" | cut -d' ' -f2-)"
  [[ -n "$log_path" && -f "$log_path" ]] || { echo "ERROR: could not find the gene-index job's log path in $stderr_file" >&2; return 1; }
  genes_line="$(grep -Eo 'spine: [0-9]+ genes' "$log_path" | tail -1)"
  rate_line="$(grep -Eo 'HGNC mapping rate: [0-9.]+%' "$log_path" | tail -1)"
  n_genes="$(echo "$genes_line" | grep -Eo '[0-9]+' | head -1)"
  mapping_pct="$(echo "$rate_line" | grep -Eo '[0-9.]+' | head -1)"
  [[ -n "$n_genes" && -n "$mapping_pct" ]] || { echo "ERROR: could not parse 'spine: N genes' / 'HGNC mapping rate: X%' out of $log_path (stderr: $stderr_file)" >&2; return 1; }
  echo "$n_genes $mapping_pct"
}

echo "control run $RUN_ID -> $OUT (commit $COMMIT)"
if [[ "${WRITE_REFERENCE:-0}" != "1" && -s "$REFERENCE" ]]; then
  REF_INPUT="$(ref_field input_sha256)"
  if [[ -n "$REF_INPUT" && "$REF_INPUT" != "$INPUT_SHA" ]]; then
    echo "ERROR: the staged inputs ($INPUT_SHA) differ from the ones the reference was produced from ($REF_INPUT)." >&2
    echo "       Ensembl or HGNC changed upstream, or a re-stage moved the GTF's dated release filename." >&2
    echo "       Re-stage (eval/stage.sh) and rerun with WRITE_REFERENCE=1 to re-baseline." >&2
    exit 1
  fi
fi

# --- positive: the two inputs, each its own job, then the gene index as a third job,
#     at the command's own default --biotype (protein_coding) ---
JOB_POS_ENSEMBL="$(bash "$SUBMIT" --name "$RUN_POS-ensembl" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
  hvantk reprocess ensembl-gene:structure --raw-dir "$RAW/ensembl-gene" \
    --output "$OUT/pos/ensembl-gene-structure.ht" --skip-download --no-check-drift | tail -1)"
JOB_POS_HGNC="$(bash "$SUBMIT" --name "$RUN_POS-hgnc" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
  hvantk reprocess hgnc:lookup --raw-dir "$RAW/hgnc" \
    --output "$OUT/pos/hgnc-lookup.ht" --skip-download --no-check-drift | tail -1)"

POS_STDERR="$OUT/pos/gi-submit.stderr"
JOB_POS_GI="$(bash "$SUBMIT" --name "$RUN_POS-gi" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
  hvantk annotate spine --gene-table "$OUT/pos/ensembl-gene-structure.ht" \
    --hgnc "$OUT/pos/hgnc-lookup.ht" --output "$OUT/pos/gene-index.ht" \
  2> "$POS_STDERR" | tail -1)"
read -r POS_GENES POS_RATE < <(parse_gi_log "$POS_STDERR")

printf 'commit: %s\ncontainer: %s\njob: %s\njob: %s\njob: %s\ngenes: %s\nmapping_rate: %s%%\n' \
  "$COMMIT" "$SIF_SHA" "$JOB_POS_ENSEMBL" "$JOB_POS_HGNC" "$JOB_POS_GI" "$POS_GENES" "$POS_RATE" \
  > "$OUT/pos/run-notes.md"

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
[[ "$(verdict_stage "$OUT/pos/verdict.json" content_ok)" == "True" && "$(verdict_stage "$OUT/pos/verdict.json" reported_ok)" == "True" ]] \
  || { echo "ERROR: positive control did not record content_ok and reported_ok both true" >&2; exit 1; }
echo "positive control: PASS (M1..M5 true, content_ok true, reported_ok true)"

# --- negative: --biotype all, reusing the SAME two built inputs (no rebuild) ---
NEG_STDERR="$OUT/neg/gi-submit.stderr"
JOB_NEG_GI="$(bash "$SUBMIT" --name "$RUN_NEG-gi" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
  hvantk annotate spine --gene-table "$OUT/pos/ensembl-gene-structure.ht" \
    --hgnc "$OUT/pos/hgnc-lookup.ht" --output "$OUT/neg/gene-index.ht" --biotype all \
  2> "$NEG_STDERR" | tail -1)"
read -r NEG_GENES NEG_RATE < <(parse_gi_log "$NEG_STDERR")

printf 'commit: %s\ncontainer: %s\njob: %s\ngenes: %s\nmapping_rate: %s%%\n' \
  "$COMMIT" "$SIF_SHA" "$JOB_NEG_GI" "$NEG_GENES" "$NEG_RATE" > "$OUT/neg/run-notes.md"

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
  || { echo "ERROR: negative control did not fail specifically on content_ok (content_ok=$(verdict_stage "$OUT/neg/verdict.json" content_ok) reported_ok=$(verdict_stage "$OUT/neg/verdict.json" reported_ok))" >&2; exit 1; }
echo "negative control: PASS (M4 true, M5 false: reported_ok true but content_ok false -- a real build, wrong biotype filter)"
echo "verdicts: $OUT/pos/verdict.json $OUT/neg/verdict.json"
