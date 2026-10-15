#!/usr/bin/env bash
# Positive and negative control for the hpc-sync eval. Git only: runs on a login node
# or a laptop, submits no job.
#
#   HVANTK_SITE_PROFILE=... [WRITE_REFERENCE=1] bash skills/hpc-sync/eval/control.sh
#
# Positive: a fresh scenario, then the skill's own procedure (fetch both remotes, prove
#   the fast-forward, push by explicit refspec from origin/<branch>, fast-forward the
#   clone), then check.py must pass.
# Negative: a fresh scenario, then the trap the skill warns about (`git push cluster
#   main dev`, which pushes the stale local branches); check.py --expect-fail must pass
#   and the verdict must show M4 false (wrong refs) with M5 true (nothing destroyed).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"
[[ "$HERE" == "$R"/* ]] || { echo "ERROR: this eval folder ($HERE) is not inside HVANTK_REPO ($R)" >&2; exit 1; }

EVAL="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/hpc-sync"
OUT="$EVAL/control/$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT/pos" "$OUT/neg"
COMMIT="$(git -C "$R" rev-parse HEAD)"
SIF_SHA="$(sha256sum "$(readlink -f "$HVANTK_SIF")" | cut -d' ' -f1)"
REFERENCE="$HERE/reference.json"
CHECK="$HERE/check.py"
verdict_stage() { python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v[sys.argv[2]])" "$1" "$2"; }
notes() { printf 'commit: %s\ncontainer: %s\n' "$COMMIT" "$SIF_SHA" > "$1"; }

echo "control run -> $OUT (commit $COMMIT)"

# --- positive: the skill's procedure, no agent ---
bash "$HERE/stage.sh" "$OUT/pos/scenario" >/dev/null
if [[ "${WRITE_REFERENCE:-0}" == "1" ]]; then
  python3 "$CHECK" --scenario "$OUT/pos/scenario" --reference "$REFERENCE" --write-reference \
    --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA"
fi
[[ -s "$REFERENCE" ]] || { echo "ERROR: $REFERENCE missing; rerun with WRITE_REFERENCE=1" >&2; exit 1; }
(
  cd "$OUT/pos/scenario/laptop"
  git fetch -q origin && git fetch -q cluster
  for b in main dev; do
    git merge-base --is-ancestor "cluster/$b" "origin/$b" || { echo "ERROR: cluster/$b is not an ancestor of origin/$b" >&2; exit 1; }
  done
  git push -q cluster origin/main:refs/heads/main origin/dev:refs/heads/dev
  git fetch -q cluster
  for b in main dev; do
    [[ "$(git rev-parse "origin/$b")" == "$(git rev-parse "cluster/$b")" ]] || { echo "ERROR: $b differs after push" >&2; exit 1; }
  done
)
notes "$OUT/pos/run-notes.md"
python3 "$CHECK" --scenario "$OUT/pos/scenario" --run-notes "$OUT/pos/run-notes.md" --reference "$REFERENCE" \
  --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" >/dev/null
for s in M1 M4 M5; do [[ "$(verdict_stage "$OUT/pos/verdict.json" $s)" == "True" ]] || { echo "ERROR: positive control $s is not true" >&2; exit 1; }; done
echo "positive control: PASS (M1 M4 M5 true; M2 M3 null: no job)"

# --- negative: the trap, stale local branches pushed by name ---
bash "$HERE/stage.sh" "$OUT/neg/scenario" >/dev/null
( cd "$OUT/neg/scenario/laptop" && git push -q cluster main dev )
notes "$OUT/neg/run-notes.md"
python3 "$CHECK" --scenario "$OUT/neg/scenario" --run-notes "$OUT/neg/run-notes.md" --reference "$REFERENCE" \
  --expect-commit "$COMMIT" --expect-sif-sha256 "$SIF_SHA" --expect-fail >/dev/null
[[ "$(verdict_stage "$OUT/neg/verdict.json" M4)" == "False" && "$(verdict_stage "$OUT/neg/verdict.json" M5)" == "True" ]] \
  || { echo "ERROR: negative control did not fail the way the trap predicts (M4=$(verdict_stage "$OUT/neg/verdict.json" M4) M5=$(verdict_stage "$OUT/neg/verdict.json" M5))" >&2; exit 1; }
echo "negative control: PASS (M4 false: the cluster got the stale local dev; M5 true: nothing destroyed)"
echo "verdicts: $OUT/pos/verdict.json $OUT/neg/verdict.json"
