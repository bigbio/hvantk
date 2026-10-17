---
name: hpc-sync
description: 'Use when a cluster clone of hvantk must be brought to a target commit without losing cluster-side commits, when branches must be pushed to a second (cluster) remote, or when checking whether the GitHub remote and the cluster remote agree. Also when a job ran old code although the branch looked current, or when the container image was built from a stale lock file.'
user-invocable: true
argument-hint: "<branches to sync, e.g. main dev> [--remote <cluster remote name>]"
metadata:
  stage: substrate
  status: provisional
  scale: laptop
---

# Keep a cluster clone in step with the GitHub remote

Read the `hvantk:workflow-conventions` skill (`skills/workflow-conventions/SKILL.md`)
first; this skill assumes its site-profile rule, its section contract and its
measurement contract.

## 1. Scope and outcome

Two remotes hold the same repository: the GitHub remote, where pull requests and CI
live, and a bare repository on the cluster, which the cluster's working clone pulls
from. Local branches track GitHub, so a bare `git push` never touches the cluster and
the cluster clone goes stale without any error. The outcome of this skill is that the
named branches on the cluster remote point at exactly the commits they point at on
GitHub, every commit that existed only on the cluster is still reachable, and the
cluster's working clone is fast-forwarded to the same commit, with the proof printed.

Use it after anything lands on a shared branch, before building a container image from
the cluster clone, and whenever a job's log names a commit that is not the one you
reviewed. Not this skill: resolving merge conflicts (that is a human decision), and
pushing to GitHub itself (that is the ordinary PR workflow).

Where each step runs: steps 1 to 4 and 6 of the procedure run in the clone you push
from (a laptop, usually), and need only the cluster remote's name; preflight step 1,
procedure step 5 and verification step 2 run on the cluster, with the cluster's site
profile. A laptop has no site profile of its own, so the steps that run there never
source one: the user gives you the cluster remote's name and it replaces
`<cluster-remote>` literally in every command below.

## 2. Preflight

1. Site profile: `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`;
   stop and ask for any missing value. Every later command that reads a profile value
   starts with the same `source "$HVANTK_SITE_PROFILE" &&`.
2. Checkout and image, for the run notes: the commit the cluster checkout is at
   (`git -C "$HVANTK_REPO" rev-parse HEAD`, on the cluster) and the image digest
   (`sha256sum "$(readlink -f "$HVANTK_SIF")"`), as every workflow records them.
3. Remotes: `git remote -v` in the clone you are pushing from. The cluster remote's
   name is a site value (`HVANTK_CLUSTER_REMOTE` in the cluster's profile when set); if
   the user has not given it and the listing does not make it obvious, ask. Never guess which remote is the
   cluster.
4. Fetch both sides before reasoning about either:
   `git fetch origin && git fetch <cluster-remote>`.
   Everything below compares remote-tracking refs, never local branches.
5. Record in your run notes the commit each branch is at on GitHub
   (`git rev-parse origin/<branch>`) and every ref the cluster remote holds
   (`git ls-remote <cluster-remote>`): the first is the target, the second is what
   must survive.

## 3. Procedure

1. Compare, branch by branch, before pushing:
   ```bash
   for b in <branches>; do
     printf '%-6s origin=%s cluster=%s\n' "$b" \
       "$(git rev-parse --short "origin/$b")" "$(git rev-parse --short "<cluster-remote>/$b" 2>/dev/null || echo none)"
   done
   ```
2. Prove each push is a fast-forward, so no cluster-side commit is clobbered:
   `git merge-base --is-ancestor <cluster-remote>/<branch> origin/<branch> && echo fast-forward`.
   If it is not, the cluster has commits GitHub does not have. Stop: list them with
   `git log --oneline origin/<branch>..<cluster-remote>/<branch>`, show them to
   the user, and let them decide whether those commits go to GitHub first (a branch and
   a PR) or are discarded. Never force-push to resolve it.
3. Push the GitHub refs by explicit refspec, never local branch names:
   ```bash
   for b in <branches>; do git push <cluster-remote> "origin/$b:refs/heads/$b"; done
   ```
   `git push <cluster-remote> <branches>` pushes your *local* branches, which can be
   weeks behind while git reports success, because sending an already-present commit
   is a legitimate no-op. The explicit form is right whether or not your local
   branches are current.
4. Do not trust the push's output. Fetch and diff the refs again (step 1); every line
   must show the same hash on both sides.
5. On the cluster, fast-forward the working clone (`$HVANTK_REPO`) from the bare remote
   it tracks; `--ff-only` makes a clone with local changes fail loudly instead of
   merging silently:
   ```bash
   source "$HVANTK_SITE_PROFILE" && git -C "$HVANTK_REPO" fetch --all --prune && git -C "$HVANTK_REPO" merge --ff-only "@{upstream}" && git -C "$HVANTK_REPO" log --oneline -1
   ```
6. Write the proof into your run notes: the branch table from step 4 and the clone's
   commit from step 5. A later job's preflight compares against exactly this.

## 4. Scale-up

- Several working clones on the cluster (one per project, or worktrees for evaluation
  runs) are each fast-forwarded separately; the bare remote is the single source.
- Tags travel the same way, by explicit refspec:
  `git push <cluster-remote> refs/tags/<tag>`.
  A release that is tagged on GitHub but not on the cluster builds an image whose
  provenance names a tag the cluster cannot resolve.
- The container image is built from `poetry.lock`. After a sync that changed the lock
  file, the image must be rebuilt before the next job, or the job runs new code on
  old dependencies. Check with `git diff --stat <old>..<new> -- poetry.lock`.
- Automating the sync (a post-merge hook or a CI step with a deploy key) is possible,
  but the fast-forward check and the ref diff stay in; silent success is the failure
  mode this skill exists to prevent.

## 5. Traps

- `git push <cluster> main dev` succeeded but the cluster did not move. Cause: it pushed
  stale local branches; the commits were already there, so git said nothing. Fix:
  explicit refspecs from `origin/<branch>` (step 3), then the ref diff (step 4).
- A job's log shows a commit you never reviewed. Cause: the cluster clone was behind,
  or ahead on a local commit. Fix: step 5, and read the clone's `git status` before
  trusting it.
- `merge --ff-only` refuses. Cause: the clone has local commits, or edits to files the
  update touches (other uncommitted edits are carried along silently; verification
  step 2 catches them). Fix: show them to the user; commit them to a branch and push
  that branch to GitHub, or discard them deliberately with `git reset --hard @{upstream}`
  only after the user says so.
- The image was built, check-install passed, and the dependency set is still the old
  one. Cause: the lock file changed upstream but the cluster clone was stale at build
  time. Fix: sync first, then build; compare the lock file's hash in the image's
  provenance record with `md5sum poetry.lock` in the synced clone.
- The push is rejected with "non-fast-forward". Cause: cluster-side commits GitHub does
  not have. Fix: step 2; never `--force`.
- `git fetch` over SSH hangs or asks for a passphrase. Cause: the key agent is empty or
  the control socket is stale. Fix: load the key into the agent and retry; a sync that
  waits on a prompt in a background job never finishes.
- The two clones carry different `origin` meanings. Cause: in the laptop clone `origin`
  is GitHub, in the cluster clone `origin` is the bare cluster repository. Fix: always
  say which clone a command runs in, and read `git remote -v` first.

## 6. Verification

1. The branch table (step 1) shows identical hashes on both sides for every branch
   named in the task.
2. `git -C "$HVANTK_REPO" rev-parse HEAD` on the cluster equals that hash for the branch
   the clone tracks, and `git -C "$HVANTK_REPO" status --porcelain` is empty.
3. Every ref that existed only on the cluster before the sync still resolves:
   `git ls-remote <cluster-remote>` lists it with the hash recorded in preflight
   step 5; and GitHub holds no ref it did not hold before.
4. `references/ref-report.md` shows the one-screen report to paste into the run notes,
   with the commit each branch is at, the fast-forward verdict, and the clone's state.
