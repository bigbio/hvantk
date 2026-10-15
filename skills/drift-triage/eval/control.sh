#!/usr/bin/env bash
# Positive and negative control for the drift-triage eval. Login-node orchestrator: it
# copies the checkout, applies the plugin repair itself (no agent), submits the builds
# and runs the bookkeeping through the image, then grades both copies with check.py.
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/drift-triage/eval/control.sh
#
# Positive: a copy of the checkout, repaired and built from the mutated raw file, with
#   the fingerprint regenerated and the ledger row recorded and marked rebuilt; check.py
#   (in the image, on a compute node) must pass.
# Negative: a second copy, repaired the same way but with NO bookkeeping, graded against
#   the POSITIVE's own built table; check.py --expect-fail must pass, showing the
#   content is fine and the code genuinely repaired (M4 true, repair_effective true)
#   while the records are not (M5 false: fingerprint and ledger row both missing).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SKILL="$(cd "$HERE/.." && pwd -P)"
SUBMIT="$SKILL/../hpc-hail-job/templates/submit-hail-job.sh"
bash "$SKILL/../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"
[[ "$HERE" == "$R"/* ]] || { echo "ERROR: this eval folder ($HERE) is not inside HVANTK_REPO ($R)" >&2; exit 1; }
[[ -z "$(git -C "$R" status --porcelain)" ]] \
  || { echo "ERROR: $R has uncommitted changes; the repair this control applies must be the only diff in the copies it makes from it" >&2; exit 1; }
# Paths handed to jobs use the profile's spelling of the checkout, which is what the
# job binds; the resolved path above is only for the containment and clean-tree checks.
REL="${HERE#"$R"/}"
CHECK="$HVANTK_REPO/$REL/check.py"
REFERENCE="$HVANTK_REPO/$REL/reference.json"

STAGE="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/drift-triage"
[[ -s "$STAGE/staged.sha256" ]] || { echo "ERROR: not staged; run eval/stage.sh first" >&2; exit 1; }
( cd "$STAGE" && sha256sum -c --quiet staged.sha256 ) \
  || { echo "ERROR: staged inputs do not match staged.sha256; delete the stage directory and run eval/stage.sh again" >&2; exit 1; }
ORIGINAL="$(find "$STAGE/raw-original/clingen" -maxdepth 1 -name '*.csv' -print -quit)"
MUTATED="$STAGE/raw/clingen/$(basename "$ORIGINAL")"
# One digest over the whole (multi-line) staged.sha256, so it changes if ANY staged
# file does -- the original download, the mutated copy, or the staged drift report --
# not just the first line of it.
INPUT_SHA="$(sha256sum "$STAGE/staged.sha256" | cut -d' ' -f1)"
OUT="$STAGE/control/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT/pos" "$OUT/neg"
RUN_ID="ctl-dt-$(date -u +%H%M%S)"
SINCE="$(date +%Y-%m-%dT%H:%M:%S)"   # local time: sacct -S reads the node's clock, not UTC
COMMIT="$(git -C "$R" rev-parse HEAD)"
SIF_SHA="$(sha256sum "$(readlink -f "$HVANTK_SIF")" | cut -d' ' -f1)"   # always the resolved file
SACCT_FIELDS="JobID,JobName,State,ReqMem,AllocCPUS,Timelimit,Partition"
CONTAINER_BIN="${HVANTK_CONTAINER_BIN:-singularity}"
BINDS=(-B "$HVANTK_WORK:$HVANTK_WORK")
for p in "$HVANTK_REPO" "$STAGE"; do
  case "$p" in "$HVANTK_WORK"|"$HVANTK_WORK"/*) ;; *) BINDS+=(-B "$p:$p") ;; esac
done
RSYNC_EXCLUDES=(--exclude .git --exclude __pycache__ --exclude local --exclude manuscript
  --exclude .pytest_cache --exclude .ruff_cache --exclude .mypy_cache --exclude .hypothesis
  --exclude '*.egg-info' --exclude .venv --exclude site)

capture_sacct() { sacct -X -n -P -S "$SINCE" -o "$SACCT_FIELDS" > "$1"; }
# A non-zero exit here is the GRADED run failing, not this script's own plumbing --
# `|| true` lets control continue to the verdict it already wrote, so the printed
# M1..M5 breakdown (next line, every call site) explains the failure before the
# explicit assertions below it decide whether that failure was expected.
check_job() { bash "$SUBMIT" --name "$1" --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- python "$CHECK" "${@:2}" || true; }
print_verdict() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print('verdict %s: %s' % (sys.argv[1], {k: v.get(k) for k in ('M1','M2','M3','M4','M5','pass')}))" "$1"; }
verdict_stage() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v[sys.argv[2]])" "$1" "$2"; }
verdict_output() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v['output'][sys.argv[2]])" "$1" "$2"; }
# Tolerates the placeholder `{"produced_by": null}`: "or {}" makes a null
# produced_by read as empty rather than raising, so the caller sees "" (never a
# match) instead of a traceback.
ref_field() { python3 -c "import json,sys; v=(json.load(open(sys.argv[1])).get('produced_by') or {}); print(v.get(sys.argv[2], ''))" "$1" "$2"; }

# Repairs `constants.py` at $1 the same way, in place, for either copy. Python-3.6-safe
# (runs on the login node, not through the image): asserts the anchor line exists
# exactly once, inserts the new rename right after it with the same indentation, and
# verifies the result before returning.
repair_constants() {
  python3 - "$1" <<'PY'
import sys

path = sys.argv[1]
anchor = '    "GENE ID (HGNC)": "hgnc_id",'
insert = '    "HGNC ID": "hgnc_id",'

with open(path) as fh:
    lines = fh.readlines()

matches = [i for i, line in enumerate(lines) if line.rstrip("\n") == anchor]
if len(matches) != 1:
    sys.exit("expected exactly one anchor line in %s, found %d" % (path, len(matches)))

i = matches[0]
ending = lines[i][len(lines[i].rstrip("\n")):] or "\n"
lines.insert(i + 1, insert + ending)

with open(path, "w") as fh:
    fh.writelines(lines)

with open(path) as fh:
    written = fh.read()
if written.count(insert) != 1:
    sys.exit("repair verification failed: %r not found exactly once in %s" % (insert, path))
print("repaired: " + path)
PY
}

echo "control run $RUN_ID -> $OUT (commit $COMMIT)"

if [[ "${WRITE_REFERENCE:-0}" != "1" ]]; then
  REF_INPUT=""
  [[ -s "$REFERENCE" ]] && REF_INPUT="$(ref_field "$REFERENCE" input_sha256)"
  if [[ -z "$REF_INPUT" ]]; then
    echo "ERROR: $REFERENCE is missing or still the {\"produced_by\": null} placeholder; rerun with WRITE_REFERENCE=1" >&2
    exit 1
  fi
  if [[ "$REF_INPUT" != "$INPUT_SHA" ]]; then
    echo "ERROR: the staged input ($INPUT_SHA) differs from the one the reference was produced from ($REF_INPUT)." >&2
    echo "       Upstream or the mutation changed. Delete the stage directory, run eval/stage.sh again, then rerun with WRITE_REFERENCE=1 to re-baseline." >&2
    exit 1
  fi
fi

# --- optional: (re)write the reference from the real checkout's own (unrepaired)
#     builder, against the ORIGINAL (unmutated) raw file ---
if [[ "${WRITE_REFERENCE:-0}" == "1" ]]; then
  JOB_REF="$(bash "$SUBMIT" --name "$RUN_ID-ref-build" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
    hvantk reprocess clingen:gene-disease --raw-dir "$(dirname "$ORIGINAL")" \
      --output "$OUT/pos/clingen-gene-disease.ref.ht" --skip-download --no-check-drift | tail -1)"
  echo "reference build: job $JOB_REF"
  check_job "$RUN_ID-ref-check" --output "$OUT/pos/clingen-gene-disease.ref.ht" --repo "$HVANTK_REPO" \
    --reference "$REFERENCE" --write-reference --staged-file "$MUTATED" \
    --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" --input-sha256 "$INPUT_SHA" >/dev/null
  [[ "$(ref_field "$REFERENCE" input_sha256)" == "$INPUT_SHA" ]] \
    || { echo "ERROR: the reference was not (re)written: the ref-check job failed; read its log" >&2; exit 1; }
  echo "reference written: $REFERENCE"
fi
[[ -s "$REFERENCE" ]] || { echo "ERROR: $REFERENCE missing; rerun with WRITE_REFERENCE=1" >&2; exit 1; }

# --- positive: a repaired copy, built from the mutated file, fingerprint regenerated
#     and ledger row recorded + marked rebuilt ---
POS_REPO="$OUT/pos/repo"
echo "preparing the positive copy ($POS_REPO) ..."
rsync -a "${RSYNC_EXCLUDES[@]}" "$HVANTK_REPO"/ "$POS_REPO/"
repair_constants "$POS_REPO/hvantk/skills/clingen/shared/constants.py"

TMP_PROFILE="$OUT/pos/site-profile.repo-copy.sh"
cp "$HVANTK_SITE_PROFILE" "$TMP_PROFILE"
printf '\n# drift-triage control: build from the repaired copy, not the real checkout\nHVANTK_REPO="%s"\n' "$POS_REPO" >> "$TMP_PROFILE"

JOB_BUILD="$(HVANTK_SITE_PROFILE="$TMP_PROFILE" bash "$SUBMIT" --name "$RUN_ID-build" --cpus 4 --mem-gb 16 --time 0-00:30:00 --wait -- \
  hvantk reprocess clingen:gene-disease --raw-dir "$(dirname "$MUTATED")" \
    --output "$OUT/pos/clingen-gene-disease.ht" --skip-download --no-check-drift | tail -1)"

# The helper (run-in-image.sh) always sets PYTHONPATH from the profile's HVANTK_REPO,
# which here still names the real checkout -- pointing it at the copy instead would
# need yet another temporary profile swap. An explicit exec is simpler for this one
# case: it is the copy's OWN code that must regenerate the copy's OWN fingerprint and
# ledger, not the real checkout's.
echo "regenerating the fingerprint and recording the ledger row in $POS_REPO ..."
REGEN_OUT="$("$CONTAINER_BIN" exec "${BINDS[@]}" --env PYTHONPATH="$POS_REPO" "$HVANTK_SIF" \
  python "$POS_REPO/skills/drift-triage/templates/probe-offline.py" drift --dataset clingen:gene-disease \
    --file "$MUTATED" -- --regenerate --accepted-in "control $RUN_ID" 2>&1)"
echo "$REGEN_OUT"
echo "$REGEN_OUT" | grep -qF "regenerated: clingen:gene-disease" \
  || { echo "ERROR: regenerate did not report 'regenerated: clingen:gene-disease'" >&2; exit 1; }
echo "$REGEN_OUT" | grep -qF "ledger row recorded" \
  || { echo "ERROR: regenerate did not record a ledger row" >&2; exit 1; }
echo "$REGEN_OUT" | grep -qF "signal: schema" \
  || { echo "ERROR: the ledger row was not recorded with signal: schema" >&2; exit 1; }

MARK_OUT="$("$CONTAINER_BIN" exec "${BINDS[@]}" --env PYTHONPATH="$POS_REPO" "$HVANTK_SIF" \
  hvantk drift --mark-rebuilt clingen:gene-disease 2>&1)"
echo "$MARK_OUT"
echo "$MARK_OUT" | grep -qF "marked rebuilt:" \
  || { echo "ERROR: mark-rebuilt did not report 'marked rebuilt:'" >&2; exit 1; }

printf 'commit: %s\ncontainer: %s\njob: %s\n' "$COMMIT" "$SIF_SHA" "$JOB_BUILD" > "$OUT/pos/run-notes.md"
capture_sacct "$OUT/pos/sacct.psv"
check_job "$RUN_ID-check" --output "$OUT/pos/clingen-gene-disease.ht" --repo "$POS_REPO" \
  --verdict-dir "$OUT/pos" --run-notes "$OUT/pos/run-notes.md" --reference "$REFERENCE" \
  --run-id "$RUN_ID" --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" \
  --site-profile "$HVANTK_SITE_PROFILE" --sacct-file "$OUT/pos/sacct.psv" >/dev/null
print_verdict "$OUT/pos/verdict.json"
for s in M1 M2 M3 M4 M5; do [[ "$(verdict_stage "$OUT/pos/verdict.json" $s)" == "True" ]] || { echo "ERROR: positive control $s is not true" >&2; exit 1; }; done
echo "positive control: PASS (M1..M5 true)"

# --- negative: the same repair, but no bookkeeping -- graded against the POSITIVE's
#     own built table, so the only variable is the checkout's recorded state ---
NEG_REPO="$OUT/neg/repo"
echo "preparing the negative copy ($NEG_REPO, same repair, no bookkeeping) ..."
rsync -a "${RSYNC_EXCLUDES[@]}" "$HVANTK_REPO"/ "$NEG_REPO/"
repair_constants "$NEG_REPO/hvantk/skills/clingen/shared/constants.py"

printf 'commit: %s\ncontainer: %s\njob: %s\n' "$COMMIT" "$SIF_SHA" "$JOB_BUILD" > "$OUT/neg/run-notes.md"
capture_sacct "$OUT/neg/sacct.psv"
check_job "$RUN_ID-negcheck" --output "$OUT/pos/clingen-gene-disease.ht" --repo "$NEG_REPO" \
  --verdict-dir "$OUT/neg" --run-notes "$OUT/neg/run-notes.md" --reference "$REFERENCE" \
  --run-id "$RUN_ID" --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" \
  --site-profile "$HVANTK_SITE_PROFILE" --sacct-file "$OUT/neg/sacct.psv" --expect-fail >/dev/null
print_verdict "$OUT/neg/verdict.json"
[[ "$(verdict_stage "$OUT/neg/verdict.json" M4)" == "True" ]] \
  || { echo "ERROR: negative control did not keep M4 true (the content is still a valid build)" >&2; exit 1; }
[[ "$(verdict_stage "$OUT/neg/verdict.json" M5)" == "False" ]] \
  || { echo "ERROR: negative control did not fail M5 (the records should be missing)" >&2; exit 1; }
[[ "$(verdict_output "$OUT/neg/verdict.json" content_ok)" == "True" \
   && "$(verdict_output "$OUT/neg/verdict.json" repair_effective)" == "True" \
   && "$(verdict_output "$OUT/neg/verdict.json" fingerprint_regenerated)" == "False" \
   && "$(verdict_output "$OUT/neg/verdict.json" ledger_row_ok)" == "False" ]] \
  || { echo "ERROR: negative control did not fail for the expected reason (content_ok/repair_effective/fingerprint_regenerated/ledger_row_ok)" >&2; exit 1; }
echo "negative control: PASS (M4 true, M5 false: content is fine and the code is genuinely repaired, but the records are not)"
echo "verdicts: $OUT/pos/verdict.json $OUT/neg/verdict.json"
