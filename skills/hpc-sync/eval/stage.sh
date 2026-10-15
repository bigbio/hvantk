#!/usr/bin/env bash
# Build the hpc-sync eval scenario: two bare remotes and a stale working clone, with
# deterministic commit hashes (fixed identity, dates and content), so the expected refs
# can be recorded once in reference.json and checked anywhere.
#
#   HVANTK_SITE_PROFILE=... bash skills/hpc-sync/eval/stage.sh [TARGET_DIR]
#
# TARGET_DIR defaults to $HVANTK_EVAL_DATA/hpc-sync/scenario and is replaced if present.
# Login-node or laptop work: git only, no network, no jobs.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"

TARGET="${1:-${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}/hpc-sync/scenario}"
rm -rf "$TARGET"
mkdir -p "$TARGET"
cd "$TARGET"

export GIT_AUTHOR_NAME=eval GIT_AUTHOR_EMAIL=eval@example.org
export GIT_COMMITTER_NAME=eval GIT_COMMITTER_EMAIL=eval@example.org
export GIT_CONFIG_GLOBAL=/dev/null GIT_CONFIG_SYSTEM=/dev/null
n=0
commit() {  # commit <file> <message>: deterministic content and timestamps
  n=$((n + 1))
  printf '%s\n' "$1 at step $n" > "$1"
  git add "$1"
  GIT_AUTHOR_DATE="2026-01-01T00:00:0${n}+0000" GIT_COMMITTER_DATE="2026-01-01T00:00:0${n}+0000" \
    git commit -q -m "$2"
}

git init -q --bare -b main github.git
git init -q --bare -b main cluster.git

git init -q -b main build
( cd build
  commit README.md "c1 readme"
  commit a.txt "c2 add a"
  commit b.txt "c3 add b"
  C3="$(git rev-parse HEAD)"
  git checkout -q -b dev
  commit c.txt "c4 add c"
  C4="$(git rev-parse HEAD)"
  commit d.txt "c5 add d"
  # GitHub: main at c3, dev at c5.
  git push -q ../github.git "refs/heads/main:refs/heads/main" "refs/heads/dev:refs/heads/dev"
  # Cluster: main at c2, dev at c3, plus its own branch on top of c3.
  git push -q ../cluster.git "${C3}~1:refs/heads/main" "${C3}:refs/heads/dev"
  git checkout -q -b cluster-work "$C3"
  commit x.txt "x cluster-side work"
  git push -q ../cluster.git "refs/heads/cluster-work:refs/heads/cluster-work"
  printf '%s\n' "$C4" > ../.stale_dev
)
rm -rf build

# The working clone: tracks github.git, has the cluster as a second remote, and its
# local dev is one commit behind GitHub, which is the trap the task is about.
git clone -q github.git laptop
( cd laptop
  git remote add cluster "$TARGET/cluster.git"
  git fetch -q cluster
  git checkout -q -b dev origin/dev
  git reset -q --hard "$(cat ../.stale_dev)"
  git checkout -q main
)
rm -f .stale_dev
echo "scenario staged under $TARGET"
echo "github:  main=$(git -C github.git rev-parse --short main) dev=$(git -C github.git rev-parse --short dev)"
echo "cluster: main=$(git -C cluster.git rev-parse --short main) dev=$(git -C cluster.git rev-parse --short dev) cluster-work=$(git -C cluster.git rev-parse --short cluster-work)"
