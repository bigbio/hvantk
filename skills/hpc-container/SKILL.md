---
name: hpc-container
description: 'Use when the hvantk container image must be built from the repository lock file on a cluster without root (fakeroot or proot), validated with the install check, recorded with provenance, or promoted as the default image. Also when a job runs old dependencies although the code is current, when a build fails with "no valid mapping entry found", when check-install fails inside a fresh image, or when the lock file changed.'
user-invocable: true
argument-hint: "[--out <image path>] [--promote] (builds from the checkout named in the site profile)"
metadata:
  stage: substrate
  status: provisional
  scale: node
---

# Build, validate and promote the container image

Read `skills/_conventions/SKILL.md` first; this skill assumes its site-profile rule,
its section contract and its measurement contract.

## 1. Scope and outcome

The image is the environment hvantk runs in on the cluster: Python 3.10, Java 11, Hail
and every locked dependency, built from `containers/hvantk.def` and the checkout's
`poetry.lock`. The outcome is a versioned image file, a provenance record next to it
(`<image>.provenance.json`) naming the commit, the lock file hash and the install check's
result, and, when asked, the default image pointer moved to it with the old default kept.

Build when the lock file changed, when the package's entry points changed, when the
base image or the def file changed, and before a release's first cluster run. Do not
build for a code change alone: jobs run the checkout's code inside the image
(`hpc-hail-job` exports `PYTHONPATH`), so the image only has to carry dependencies.

Not this skill: running jobs (that is `hpc-hail-job`), editing the def file (a code
change with its own review), and clusters where images are built by administrators
only (then the "build elsewhere and copy" note in section 4 applies).

## 2. Preflight

1. Site profile: `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`;
   stop and ask for any missing value. Every later command that reads a profile value
   starts with the same `source "$HVANTK_SITE_PROFILE" &&`.
2. Checkout: `git -C "$HVANTK_REPO" rev-parse HEAD` and `git -C "$HVANTK_REPO" status --porcelain`.
   Build from a clean tree at the commit the user expects; a dirty or stale tree builds
   an image whose provenance names a commit that is not what went in. If the clone is
   behind, run the `hpc-sync` skill first.
3. Unprivileged build method, once per cluster: `grep -c "^$USER:" /etc/subuid` greater
   than zero means `--fakeroot` works; zero means the build needs a static `proot` on
   `PATH` (`references/build-troubleshooting.md` says how to get one). Record which.
4. Egress: the build pulls the base image and wheels, so the node that builds needs
   outbound access. `HVANTK_EGRESS` says whether compute nodes have it; if they do not,
   build on the login node only if the site allows it, otherwise build elsewhere and copy
   the image file (section 4).
5. Space: cache and temporary files go under the work root (the template puts them in
   `$HVANTK_WORK/containers/{cache,tmp}`), never under the home directory; a build needs
   several GB there and the image is about 2 GB.
6. Container binary: `HVANTK_CONTAINER_BIN` (`singularity` or `apptainer`); the template
   uses it for the build and for the check.

## 3. Procedure

1. Choose the image name from the commit or tag: `hvantk-<tag or short sha>.sif` under
   `$HVANTK_WORK/containers/`. Never build straight onto the default name.
2. Submit the build from the login node, which captures the git facts the compute node
   cannot (compute nodes often have no git):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-container/templates/submit-build-image.sh" \
     --out "$HVANTK_WORK/containers/hvantk-<id>.sif" --name <run>-build --wait
   ```
   The job builds from the repository root (the def file's `%files` are relative),
   forces an overwrite of the named file, gates the Java and Hail versions inside the
   image, runs `hvantk utils check-install` through the scratch pattern, and writes the
   provenance record. Allow about twenty minutes on a shared partition; the request is
   8 cores, 32 GB and 1h45 so a slow pull does not kill it.
3. Read the job log and the record: `cat "<image>.provenance.json"`. The fields are in
   `references/provenance.md`; `check_install_rc` must be 0 and `hail_version` must be
   the locked one.
4. Promote only when asked, and only after step 3:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-container/templates/promote-image.sh" "$HVANTK_WORK/containers/hvantk-<id>.sif"
   ```
   It refuses while a job of yours is running, keeps the old default as
   `hvantk-<date>.sif` with its record, points `hvantk.sif` (and its record) at the new
   file, and submits a verification job that runs the install check through the
   promoted name. Say in the run notes which image is now the default.
5. Write the run notes: commit, image digest (`sha256sum` of the new file), job ids, and
   whether the default moved.

## 4. Scale-up

- Several images coexist under `$HVANTK_WORK/containers/`: one per release or per lock
  file change, each with its provenance record. Old defaults keep their dated names
  until the user removes them; a job that still names one keeps working.
- Rebuild cadence: on every lock file change, and at least at every release so the
  image's `git_commit` is a tag that the cluster can resolve (sync tags with `hpc-sync`).
- Building elsewhere: on a laptop or a build host with the same container tool, run
  `singularity build hvantk-<id>.sif containers/hvantk.def` from the repository root
  (with `--fakeroot` or `proot` as the host allows), copy the file and write the
  provenance record by hand from the same fields; then run the verification job on the
  cluster, because an image that was never checked on a compute node is not validated.
- A site that proxies the registry or pins its own base images: the proxy and registry
  variables are shell environment for the build job. Export them in the site profile,
  which the template sources, never in a committed file.
- Cache hygiene: `$HVANTK_WORK/containers/cache` grows with every pull; clear it
  between releases, never during a build.

## 5. Traps

- `could not use fakeroot: no valid mapping entry found for <user>` or `--remote,
  --fakeroot, or the proot command are required`. Cause: no `/etc/subuid` entry.
  Fix: a static `proot` on `PATH` (section 2 step 3); the template detects which to use.
- The build fails at `%files` with files not found. Cause: it was not started from the
  repository root; the def file copies relative paths. Fix: the template `cd`s to
  `$HVANTK_REPO`; if you build by hand, do the same.
- `apt-get` finds the index but every package 404s. Cause: a base image whose
  distribution left support (the old Debian bullseye base). Fix: the def file pins a
  supported base; do not "fix" it by pointing apt at a snapshot.
- The build wants `openjdk-11` and the distribution offers only 17. Cause: a base
  without Java 11 (Debian bookworm). Fix: keep the def file's base; Hail supports Java
  8 and 11 only.
- `check-install` fails with `Connection refused` right after the build. Cause: it was
  run on the login node or without node-local scratch. Fix: the template runs it in the
  job through the scratch pattern; by hand, use `containers/hvantk_run.sh` on a compute
  node.
- The build dies after pulling nothing, or hangs on the registry. Cause: no egress from
  the node. Fix: section 2 step 4; the template probes the registry first and says so.
- The home quota fills during a build. Cause: cache or temporary files under the home
  directory. Fix: the template sets the cache and tmp directories under the work root;
  by hand, export `SINGULARITY_CACHEDIR` and `SINGULARITY_TMPDIR` (or the `APPTAINER_`
  pair) before building.
- The image passes every check and the dependency set is still old. Cause: the clone
  was stale when the lock file was copied in. Fix: `hpc-sync` first; compare
  `poetry_lock_md5` in the record with `md5sum poetry.lock` in the synced clone.
- Promotion moved the symlink while a job was running. Cause: no check for running
  jobs. Fix: the promote script refuses while `squeue -u "$USER"` lists anything; wait.
- `singularity` is not on `PATH` but `apptainer` is (or the reverse). Cause: site
  naming. Fix: `HVANTK_CONTAINER_BIN` in the profile.
- The record's `host` or `slurm_job` field is empty. Cause: the build ran outside a job.
  Fix: fine for a build elsewhere; say so in the run notes.

## 6. Verification

1. `<image>.provenance.json` exists, `check_install_rc` is 0, `git_commit` is the
   commit you built from, `worktree` is `clean`, `poetry_lock_md5` equals
   `md5sum "$HVANTK_REPO/poetry.lock"`.
2. `hail_version` in the record equals the locked Hail version, and `java_version`
   starts with `11`.
3. A verification job through the final name (`hvantk.sif` after promotion, the
   versioned name otherwise):
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" --sif "$HVANTK_WORK/containers/hvantk-<id>.sif" --name <run>-verify --cpus 4 --mem-gb 16 --time 0-00:20:00 --wait -- hvantk utils check-install
   ```
   prints the Hail version and `All checks passed.`
4. The run notes carry the commit, the image digest and the job ids; `references/provenance.md`
   shows the record and how each field is read.
