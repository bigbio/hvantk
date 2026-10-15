#!/usr/bin/env bash
# Positive and negative control for the cohort-attach-burden eval. Login-node
# orchestrator: it runs the skill's own procedure through its templates (no agent),
# then grades with check.py.
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/cohort-attach-burden/eval/control.sh
#
# Positive: burden (--carrier-mode het --mtc bh) -> edit the emitted manifest's name
#   -> validate -> attach onto the staged gene index; check.py must pass.
# Negative: the identical procedure on the SAME cohort MatrixTable, but
#   --carrier-mode hom instead of het. The burden op groups by (gene, route)
#   independently of carrier mode, so the row count and the gene universe are
#   unchanged; only the per-gene statistics (minp, odds_ratio, the architecture
#   reductions) differ. check.py --expect-fail must pass with M4 true, M5 false,
#   reported_ok true and burden_ok false.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SKILL="$(cd "$HERE/.." && pwd -P)"
SUBMIT="$SKILL/../hpc-hail-job/templates/submit-hail-job.sh"
RUN_IN_IMAGE="$SKILL/../hpc-hail-job/templates/run-in-image.sh"
[[ -f "$SUBMIT" ]] || { echo "ERROR: submit-hail-job.sh not found: $SUBMIT" >&2; exit 2; }
bash "$SKILL/../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"
[[ "$HERE" == "$R"/* ]] || { echo "ERROR: this eval folder ($HERE) is not inside HVANTK_REPO ($R)" >&2; exit 1; }
# Paths handed to a job (anything run with `hvantk`/`python ...` inside the container)
# use the profile's spelling of the checkout, which is what the job binds; the
# resolved $R above is only for the containment test.
REL="${HERE#"$R"/}"
CHECK="$HVANTK_REPO/$REL/check.py"
REFERENCE="$HVANTK_REPO/$REL/reference.json"

EVAL_ROOT="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}"
STAGE="$EVAL_ROOT/cohort-attach-burden"
FILE="CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz"
[[ -s "$STAGE/staged.sha256" ]] || { echo "ERROR: not staged; run eval/stage.sh first" >&2; exit 1; }
( cd "$STAGE" && sha256sum -c --quiet staged.sha256 ) \
  || { echo "ERROR: staged inputs do not match staged.sha256; delete the stage directory and run eval/stage.sh again" >&2; exit 1; }
# sha256sum -c only checks that the files IT lists still match; it says nothing about
# an extra file dropped into the stage directory since staging.
EXPECTED_PATHS="$(cut -d' ' -f3- "$STAGE/staged.sha256" | sort)"
ACTUAL_PATHS="$(cd "$STAGE" && printf '%s\n' "raw/$FILE" "raw/$FILE.tbi" cohort.summary.json upstream/gene-index.receipt.json | sort)"
[[ "$EXPECTED_PATHS" == "$ACTUAL_PATHS" ]] \
  || { echo "ERROR: the stage directory no longer holds exactly the files staged.sha256 lists; delete the stale files and run eval/stage.sh again" >&2; echo "  expected: $EXPECTED_PATHS" >&2; echo "  actual:   $ACTUAL_PATHS" >&2; exit 1; }
[[ -e "$STAGE/gene-index.ht/_SUCCESS" && -e "$STAGE/cohort.mt/_SUCCESS" ]] \
  || { echo "ERROR: gene-index.ht or cohort.mt missing under $STAGE; run eval/stage.sh again" >&2; exit 1; }
INPUT_SHA="$(sha256sum "$STAGE/staged.sha256" | cut -d' ' -f1)"
OUT="$STAGE/control/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT/pos" "$OUT/neg"
BASE="ctl-cab-$(date -u +%H%M%S)"
RUN_POS="$BASE-p"
RUN_NEG="$BASE-n"
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
print('verdict:', {k: v.get(k) for k in ('M1','M2','M3','M4','M5','pass','burden_ok','attach_ok','reported_ok')})" "$1"; }
# Renames the manifest's `name:` line in place (3.6-safe); --emit-manifest's default is
# the --output path's stem ("gene_burden"), the task wants "eval-cohort". Everything
# else is kept exactly as rendered by render_cohort_manifest (verified against
# hvantk/algorithms/burden/manifest.py before this script was written: key_column
# defaults to "gene", prior.column "minp"/"lower_is_better", one `architecture` axis
# holding n_case_var/conc/driver_af -- already what cohort.yaml needs).
rename_manifest() {
  python3 - "$1" <<'PY'
import re
import sys

path = sys.argv[1]
with open(path) as fh:
    text = fh.read()
new_text, n = re.subn(r"^name: .*$", "name: eval-cohort", text, count=1, flags=re.M)
if n != 1:
    sys.exit("expected exactly one 'name: ...' line in %s, found %d" % (path, n))
with open(path, "w") as fh:
    fh.write(new_text)
print("renamed: " + path)
PY
}
# The attach job's own click.echo'd "attach NAME: n/m genes tested -> ..." line lands
# in the Slurm job's log file, not in submit-hail-job.sh's own stdout; the log path is
# submit-hail-job.sh's own "log: <path>" line on STDERR (see --wait's polling block).
attached_fraction() {  # attached_fraction <stderr-file>: prints "n/m" or exits 1
  local log
  log="$(grep -oE '^log: .*' "$1" | sed 's/^log: //')"
  [[ -n "$log" && -f "$log" ]] || { echo "ERROR: could not find the job's log via submit-hail-job.sh's stderr ($1)" >&2; return 1; }
  local line
  line="$(grep -E '^attach .*: [0-9]+/[0-9]+ genes tested ->' "$log" | tail -1)"
  [[ -n "$line" ]] || { echo "ERROR: 'attach NAME: n/m genes tested' line not found in $log" >&2; cat "$log" >&2; return 1; }
  echo "$line" | grep -oE '[0-9]+/[0-9]+' | head -1
}

echo "control run $BASE -> $OUT (commit $COMMIT)"
if [[ "${WRITE_REFERENCE:-0}" != "1" && -s "$REFERENCE" ]]; then
  REF_INPUT="$(ref_field input_sha256)"
  if [[ -n "$REF_INPUT" && "$REF_INPUT" != "$INPUT_SHA" ]]; then
    echo "ERROR: the staged inputs ($INPUT_SHA) differ from the ones the reference was produced from ($REF_INPUT); re-stage (eval/stage.sh) and rerun with WRITE_REFERENCE=1 to re-baseline" >&2
    exit 1
  fi
fi

# --- positive: burden (het, bh) -> manifest -> validate -> attach ---
JOB_POS_BURDEN="$(bash "$SUBMIT" --name "$RUN_POS-burden" --cpus 8 --mem-gb 32 --time 0-00:30:00 --wait -- \
  hvantk cohort burden --mt "$STAGE/cohort.mt" --gene-col gene_id --route-col route --arm-col is_case \
    --key gene_id --carrier-mode het --mtc bh --output "$OUT/pos/gene_burden.tsv" \
    --emit-manifest "$OUT/pos/cohort.yaml" | tail -1)"
rename_manifest "$OUT/pos/cohort.yaml"
bash "$RUN_IN_IMAGE" -- hvantk cohort validate --cohort "$OUT/pos/cohort.yaml"

JOB_POS_ATTACH="$(bash "$SUBMIT" --name "$RUN_POS-attach" --cpus 8 --mem-gb 32 --time 0-00:30:00 --wait -- \
  hvantk cohort attach --cohort "$OUT/pos/cohort.yaml" --layer1 "$STAGE/gene-index.ht" \
    --output "$OUT/pos/attached.ht" 2>"$OUT/pos/attach.submit.stderr" | tail -1)"
ATTACHED_POS="$(attached_fraction "$OUT/pos/attach.submit.stderr")"
printf 'commit: %s\ncontainer: %s\njob: %s\njob: %s\nattached: %s\n' \
  "$COMMIT" "$SIF_SHA" "$JOB_POS_BURDEN" "$JOB_POS_ATTACH" "$ATTACHED_POS" > "$OUT/pos/run-notes.md"

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

# --- negative: the same procedure, --carrier-mode hom instead of het ---
JOB_NEG_BURDEN="$(bash "$SUBMIT" --name "$RUN_NEG-burden" --cpus 8 --mem-gb 32 --time 0-00:30:00 --wait -- \
  hvantk cohort burden --mt "$STAGE/cohort.mt" --gene-col gene_id --route-col route --arm-col is_case \
    --key gene_id --carrier-mode hom --mtc bh --output "$OUT/neg/gene_burden.tsv" \
    --emit-manifest "$OUT/neg/cohort.yaml" | tail -1)"
rename_manifest "$OUT/neg/cohort.yaml"
bash "$RUN_IN_IMAGE" -- hvantk cohort validate --cohort "$OUT/neg/cohort.yaml"

JOB_NEG_ATTACH="$(bash "$SUBMIT" --name "$RUN_NEG-attach" --cpus 8 --mem-gb 32 --time 0-00:30:00 --wait -- \
  hvantk cohort attach --cohort "$OUT/neg/cohort.yaml" --layer1 "$STAGE/gene-index.ht" \
    --output "$OUT/neg/attached.ht" 2>"$OUT/neg/attach.submit.stderr" | tail -1)"
ATTACHED_NEG="$(attached_fraction "$OUT/neg/attach.submit.stderr")"
printf 'commit: %s\ncontainer: %s\njob: %s\njob: %s\nattached: %s\n' \
  "$COMMIT" "$SIF_SHA" "$JOB_NEG_BURDEN" "$JOB_NEG_ATTACH" "$ATTACHED_NEG" > "$OUT/neg/run-notes.md"

capture_sacct "$OUT/neg/sacct.psv"
bash "$SUBMIT" --name "$RUN_NEG-check" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
  python "$CHECK" --out-dir "$OUT/neg" --run-notes "$OUT/neg/run-notes.md" \
    --reference "$REFERENCE" --run-id "$RUN_NEG" --expect-commit "$COMMIT" \
    --expect-sif-sha256 "$SIF_SHA" --site-profile "$HVANTK_SITE_PROFILE" \
    --sacct-file "$OUT/neg/sacct.psv" --expect-fail >/dev/null || true
print_verdict "$OUT/neg/verdict.json"
[[ "$(verdict_stage "$OUT/neg/verdict.json" M4)" == "True" && "$(verdict_stage "$OUT/neg/verdict.json" M5)" == "False" ]] \
  || { echo "ERROR: negative control did not fail the expected way (M4=$(verdict_stage "$OUT/neg/verdict.json" M4) M5=$(verdict_stage "$OUT/neg/verdict.json" M5))" >&2; exit 1; }
[[ "$(verdict_stage "$OUT/neg/verdict.json" reported_ok)" == "True" && "$(verdict_stage "$OUT/neg/verdict.json" burden_ok)" == "False" ]] \
  || { echo "ERROR: negative control did not fail specifically on burden content (reported_ok=$(verdict_stage "$OUT/neg/verdict.json" reported_ok) burden_ok=$(verdict_stage "$OUT/neg/verdict.json" burden_ok))" >&2; exit 1; }
echo "negative control: PASS (M4 true, M5 false: reported_ok true, burden_ok false -- same gene/route universe, different carrier-mode statistics)"
echo "verdicts: $OUT/pos/verdict.json $OUT/neg/verdict.json"
echo "pos/ holds gene_burden.tsv and cohort.yaml for rerank-evaluate's stage to copy"
