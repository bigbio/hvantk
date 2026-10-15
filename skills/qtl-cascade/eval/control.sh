#!/usr/bin/env bash
# Positive and negative control for the qtl-cascade eval. Login-node orchestrator: it
# only submits jobs through the skill's own templates and waits on them.
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/qtl-cascade/eval/control.sh
#
# Positive: hgnc, eqtl (all pairs) and pqtl (all pairs) built as three jobs, then the
#   cascade run (with coloc) as a fourth job; check.py must pass, proving M1..M5 true
#   (with WRITE_REFERENCE=1, first write reference.json).
# Negative: the SAME eqtl table and hgnc table, but the pQTL table rebuilt at
#   p_threshold=5e-8 (significant pairs only) and passed as BOTH --pqtl-ht and
#   --pqtl-allpairs -- the significant-pairs trap (skills/qtl-cascade/SKILL.md,
#   Traps). Cascade classification re-applies that same threshold internally (a
#   no-op on an already-filtered table), so the cascade output is unaffected; coloc
#   needs every pair in the window, so it is not. check.py --expect-fail must pass,
#   with M4 true and M5 false because cascade_ok is true but coloc_ok is false.
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
# Paths handed to a job (anything run with `hvantk`/`python ...` inside the
# container) use the profile's spelling of the checkout, which is what the job
# binds; the resolved $R above is only for the containment test.
REL="${HERE#"$R"/}"
CHECK="$HVANTK_REPO/$REL/check.py"
REFERENCE="$HVANTK_REPO/$REL/reference.json"

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/qtl-cascade"
RAW="$STAGE/raw"
[[ -s "$STAGE/staged.sha256" ]] || { echo "ERROR: not staged; run eval/stage.sh first" >&2; exit 1; }
HGNC_RAW="$RAW/hgnc"
EQTL_RAW="$RAW/gtex-eqtl"
PQTL_RAW="$RAW/pqtl"
INPUT_SHA="$(sha256sum "$STAGE/staged.sha256" | cut -d' ' -f1)"
OUT="$STAGE/control/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT/pos" "$OUT/neg"
BASE="ctl-qc-$(date -u +%H%M%S)"
RUN_POS="$BASE-p"
RUN_NEG="$BASE-n"
SINCE="$(date +%Y-%m-%dT%H:%M:%S)"   # local time: sacct -S reads the node's clock, not UTC
COMMIT="$(git -C "$R" rev-parse HEAD)"
SIF_SHA="$(sha256sum "$(readlink -f "$HVANTK_SIF")" | cut -d' ' -f1)"   # always the resolved file
SACCT_FIELDS="JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"
LOG_DIR="${HVANTK_LOG_DIR:-$HVANTK_WORK/hvantk-logs}"

capture_sacct() { sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$1"; }
# `produced_by` is null in the placeholder reference.json ({"produced_by": null})
# the implementer commits before any control has run; `or {}` reads that the same as
# a reference with no matching field, rather than crashing on None has no .get.
ref_field() { python3 -c "import json,sys
d=json.load(open(sys.argv[1]))
p=d.get('produced_by') or {}
print(p.get(sys.argv[2], ''))" "$REFERENCE" "$1"; }
verdict_stage() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v.get(sys.argv[2]))" "$1" "$2"; }
output_stage() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v.get('output', {}).get(sys.argv[2]))" "$1" "$2"; }
print_verdict() { python3 -c "import json,sys
v=json.load(open(sys.argv[1]))
print('verdict:', {k: v.get(k) for k in ('M1','M2','M3','M4','M5','pass')},
      {k: v.get('output', {}).get(k) for k in ('cascade_ok','coloc_ok','reported_ok')})" "$1"; }
# The cascade job's own log carries `qtlcascade run`'s printed "  Liver: N genes, M
# colocalised" line (its click.echo output lands in the job's redirected log file,
# never in this script's own captured stdout). Parsed on the login node, Python
# 3.6-safe.
report_counts() { python3 -c "import re,sys
text = open(sys.argv[1]).read()
m = re.search(r'^  ' + re.escape(sys.argv[2]) + r': (\d+) genes, (\d+) colocalised\$', text, re.M)
if not m:
    sys.exit('no \"  ' + sys.argv[2] + ': N genes, M colocalised\" line in ' + sys.argv[1])
print(m.group(1))
print(m.group(2))" "$1" "$2"; }

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
    echo "       The eQTL Catalogue window may have changed upstream, or the HGNC file was re-staged." >&2
    echo "       Re-stage (eval/stage.sh) and rerun with WRITE_REFERENCE=1 to re-baseline." >&2
    exit 1
  fi
fi

# --- positive: hgnc + eqtl + pqtl (all-pairs) builds, then the cascade run ---
HGNC_HT="$OUT/pos/hgnc-lookup.ht"
EQTL_HT="$OUT/pos/eqtl-liver.ht"
PQTL_HT="$OUT/pos/pqtl-liver.ht"
CASCADE_OUT="$OUT/pos/cascade"

JOB_HGNC="$(bash "$SUBMIT" --name "$RUN_POS-hgnc" --cpus 2 --mem-gb 8 --time 0-00:20:00 --wait -- \
  hvantk reprocess hgnc:lookup --raw-dir "$HGNC_RAW" --output "$HGNC_HT" \
    --skip-download --no-check-drift | tail -1)"
JOB_EQTL="$(bash "$SUBMIT" --name "$RUN_POS-eqtl" --cpus 4 --mem-gb 16 --time 0-00:45:00 --wait -- \
  hvantk reprocess gtex-eqtl:eqtls --raw-dir "$EQTL_RAW" --output "$EQTL_HT" \
    --skip-download --no-check-drift \
    --plugin-arg source=gtex_v8 --plugin-arg tissue=Liver --plugin-arg p_threshold=0 | tail -1)"
JOB_PQTL="$(bash "$SUBMIT" --name "$RUN_POS-pqtl" --cpus 4 --mem-gb 16 --time 0-00:45:00 --wait -- \
  hvantk reprocess pqtl:metrics --raw-dir "$PQTL_RAW" --output "$PQTL_HT" \
    --skip-download --no-check-drift \
    --plugin-arg tissue=Liver --plugin-arg hgnc_ht="$HGNC_HT" | tail -1)"
JOB_CASCADE="$(bash "$SUBMIT" --name "$RUN_POS-cascade" --cpus 4 --mem-gb 16 --time 0-00:45:00 --wait -- \
  hvantk qtlcascade run --eqtl-ht "$EQTL_HT" --pqtl-ht "$PQTL_HT" \
    --eqtl-allpairs "$EQTL_HT" --pqtl-allpairs "$PQTL_HT" \
    --tissues Liver -o "$CASCADE_OUT" | tail -1)"
mapfile -t POS_COUNTS < <(report_counts "$LOG_DIR/$RUN_POS-cascade_$JOB_CASCADE.out" Liver)
printf 'commit: %s\ncontainer: %s\njob: %s\njob: %s\njob: %s\njob: %s\ncascade_genes: %s\ncolocalised: %s\n' \
  "$COMMIT" "$SIF_SHA" "$JOB_HGNC" "$JOB_EQTL" "$JOB_PQTL" "$JOB_CASCADE" \
  "${POS_COUNTS[0]}" "${POS_COUNTS[1]}" > "$OUT/pos/run-notes.md"

if [[ "${WRITE_REFERENCE:-0}" == "1" ]]; then
  bash "$SUBMIT" --name "$RUN_POS-ref" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
    python "$CHECK" --out-dir "$CASCADE_OUT" --tissue Liver --reference "$REFERENCE" \
      --write-reference --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" \
      --input-sha256 "$INPUT_SHA"
  echo "reference written: $REFERENCE"
fi
[[ -s "$REFERENCE" ]] || { echo "ERROR: $REFERENCE missing; rerun with WRITE_REFERENCE=1" >&2; exit 1; }

capture_sacct "$OUT/pos/sacct.psv"
run_check "$RUN_POS-check" "$CASCADE_OUT" --out-dir "$CASCADE_OUT" --tissue Liver \
  --run-notes "$OUT/pos/run-notes.md" --reference "$REFERENCE" --run-id "$RUN_POS" \
  --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" \
  --site-profile "$HVANTK_SITE_PROFILE" --sacct-file "$OUT/pos/sacct.psv"
for s in M1 M2 M3 M4 M5; do [[ "$(verdict_stage "$CASCADE_OUT/verdict.json" $s)" == "True" ]] || { echo "ERROR: positive control $s is not true" >&2; exit 1; }; done
echo "positive control: PASS (M1..M5 true)"

# --- negative: the same eqtl/hgnc tables, pQTL rebuilt significant-pairs-only and
# --- passed as BOTH --pqtl-ht and --pqtl-allpairs ---
PQTL_SIGNIF_HT="$OUT/neg/pqtl-liver-signif.ht"
CASCADE_OUT_NEG="$OUT/neg/cascade"

JOB_PQTL_SIGNIF="$(bash "$SUBMIT" --name "$RUN_NEG-pqtl" --cpus 4 --mem-gb 16 --time 0-00:45:00 --wait -- \
  hvantk reprocess pqtl:metrics --raw-dir "$PQTL_RAW" --output "$PQTL_SIGNIF_HT" \
    --skip-download --no-check-drift \
    --plugin-arg tissue=Liver --plugin-arg hgnc_ht="$HGNC_HT" --plugin-arg p_threshold=5e-8 | tail -1)"
JOB_CASCADE_NEG="$(bash "$SUBMIT" --name "$RUN_NEG-cascade" --cpus 4 --mem-gb 16 --time 0-00:45:00 --wait -- \
  hvantk qtlcascade run --eqtl-ht "$EQTL_HT" --pqtl-ht "$PQTL_SIGNIF_HT" \
    --eqtl-allpairs "$EQTL_HT" --pqtl-allpairs "$PQTL_SIGNIF_HT" \
    --tissues Liver -o "$CASCADE_OUT_NEG" | tail -1)"
mapfile -t NEG_COUNTS < <(report_counts "$LOG_DIR/$RUN_NEG-cascade_$JOB_CASCADE_NEG.out" Liver)
printf 'commit: %s\ncontainer: %s\njob: %s\njob: %s\ncascade_genes: %s\ncolocalised: %s\n' \
  "$COMMIT" "$SIF_SHA" "$JOB_PQTL_SIGNIF" "$JOB_CASCADE_NEG" \
  "${NEG_COUNTS[0]}" "${NEG_COUNTS[1]}" > "$OUT/neg/run-notes.md"

capture_sacct "$OUT/neg/sacct.psv"
run_check "$RUN_NEG-check" "$CASCADE_OUT_NEG" --out-dir "$CASCADE_OUT_NEG" --tissue Liver \
  --run-notes "$OUT/neg/run-notes.md" --reference "$REFERENCE" --run-id "$RUN_NEG" \
  --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" \
  --site-profile "$HVANTK_SITE_PROFILE" --sacct-file "$OUT/neg/sacct.psv" --expect-fail
[[ "$(verdict_stage "$CASCADE_OUT_NEG/verdict.json" M4)" == "True" && "$(verdict_stage "$CASCADE_OUT_NEG/verdict.json" M5)" == "False" ]] \
  || { echo "ERROR: negative control did not fail the expected way (M4=$(verdict_stage "$CASCADE_OUT_NEG/verdict.json" M4) M5=$(verdict_stage "$CASCADE_OUT_NEG/verdict.json" M5))" >&2; exit 1; }
[[ "$(output_stage "$CASCADE_OUT_NEG/verdict.json" cascade_ok)" == "True" && "$(output_stage "$CASCADE_OUT_NEG/verdict.json" coloc_ok)" == "False" ]] \
  || { echo "ERROR: negative control did not fail specifically on coloc_ok (cascade_ok=$(output_stage "$CASCADE_OUT_NEG/verdict.json" cascade_ok) coloc_ok=$(output_stage "$CASCADE_OUT_NEG/verdict.json" coloc_ok))" >&2; exit 1; }
echo "negative control: PASS (M4 true, M5 false: cascade_ok true but coloc_ok false -- significant-pairs-only coloc input)"
echo "verdicts: $CASCADE_OUT/verdict.json $CASCADE_OUT_NEG/verdict.json"
