---
name: hpc-hail-job
description: 'Use when running any Hail-backed hvantk command (reprocess, hgc, annotate, cohort, qtlcascade, utils check-install) on a Slurm cluster from the container image: sizing one node, node-local scratch, submitting, waiting, triaging a pending or failed job, and verifying the output. Also when a Hail job dies with "Connection refused", "Failed to create any local dir", a cgroup out-of-memory kill, or "Executor heartbeat timed out".'
stage: substrate
status: provisional
scale: node
user-invocable: true
argument-hint: "<hvantk subcommand and its arguments> [--cpus N --mem-gb G --time D-HH:MM:SS]"
---

# Run a Hail-backed hvantk command as a Slurm job

Read `skills/_conventions/SKILL.md` first; this skill assumes its site-profile rule,
its section contract and its measurement contract.

## 1. Scope and outcome

One hvantk command, one Slurm job, one node, Hail in Spark local mode inside the
container image. The outcome is the command's output (for a dataset build, a Hail table
or matrix table directory plus its `.provenance.json` sidecar), the job log under the
log directory, and a line in your run notes with the commit, the image digest and the
job id.

Use it for every Hail-backed command: `reprocess`, `hgc`, `annotate`, `cohort`,
`qtlcascade`, `utils check-install`. Use it for small commands too. Nothing runs on the
login node except `git`, file inspection, downloads and job submission.

Not this skill: building the image (a separate skill), multi-node Spark (section 4
says when, and that it is rarely worth it), commands with no Hail in them (they still
run fine through the same template, with a small request).

## 2. Preflight

1. Site profile: `source "$HVANTK_SITE_PROFILE"` and run
   `bash skills/_site/check-site-profile.sh`. If it reports a missing value, stop and
   ask the user for it. Never guess a partition, a path or an image.
2. Checkout: `git -C "$HVANTK_REPO" rev-parse HEAD` and `git -C "$HVANTK_REPO" status --porcelain`.
   Compare with the commit the user expects. A clone that is behind runs old code with
   no error of any kind; if it is behind, say so before submitting anything. The job
   runs the checkout's code inside the image's environment (the template exports
   `PYTHONPATH`), so the commit decides what runs; the image only has to carry the
   dependencies, and needs rebuilding only when the lock file changes.
3. Image: `readlink -f "$HVANTK_SIF"` (the configured path is often a symlink to a
   versioned image) and its digest: `$HVANTK_SIF_SHA256` if the profile has it, else
   `sha256sum "$(readlink -f "$HVANTK_SIF")"` (about ten seconds). Record both.
4. Inputs and outputs: inputs must be readable from compute nodes and the output must
   be written under `$HVANTK_WORK` or another path that is writable from compute nodes.
   Never write under the home directory (quota) and never under a path listed in
   `$HVANTK_READONLY_MOUNTS`.
5. Egress: if `$HVANTK_EGRESS` is `blocked` or unset, treat compute nodes as offline.
   Downloads happen on the login node first (`hvantk download <source> ...`), and the
   build runs with `--skip-download --no-check-drift`.
6. Once per image and node type, prove the stack works:
   `bash skills/hpc-hail-job/templates/submit-hail-job.sh --name check --cpus 4 --mem-gb 16 --time 0-00:20:00 --wait -- hvantk utils check-install`
   must print the Hail version and pass its smoke test.

## 3. Procedure

1. Size the job with `references/sizing.md`. Starting points: a gene-level table
   4 cores, 16 GB, 30 minutes; a variant-level table 16 cores, 120 GB, 4 hours; joint
   genotyping the whole node, exclusive, 12 hours. Never request more memory than
   `$HVANTK_NODE_MEM_GB` or more cores than `$HVANTK_NODE_CPUS`.
2. Dry-check the submission: the same command with `--dry-run` runs `sbatch --test-only`
   and prints the resolved request without submitting.
   ```bash
   bash skills/hpc-hail-job/templates/submit-hail-job.sh --name <run>-<what> \
     --cpus 4 --mem-gb 16 --time 0-00:30:00 --dry-run -- \
     hvantk reprocess hgnc:lookup --raw-dir "$HVANTK_WORK/raw/hgnc" \
       --output "$HVANTK_WORK/out/hgnc.ht" --skip-download --no-check-drift
   ```
3. Submit: drop `--dry-run`. The script prints the job id; write `job: <id>` into your
   run notes at once, with `commit:` and `container:` from the preflight.
4. Wait on the job id, never on a log string:
   ```bash
   until [ -z "$(squeue -j <id> -h -o %T)" ]; do sleep 60; done
   sacct -j <id> -X -o JobID,JobName%30,State,ExitCode,Elapsed,ReqMem,AllocCPUS
   ```
   `--wait` on the submit script does exactly this and exits non-zero unless the job
   reached `COMPLETED`.
5. Read the log at `$HVANTK_LOG_DIR/<name>_<id>.out`. A job that failed goes to
   `references/triage.md`; fix the request or the input and resubmit under a new name.
   Do not rerun on the login node "to see the error faster".
6. Verify the output (section 6) and finish the run notes.

The template `templates/hail-job.sbatch` is the batch form of
`containers/hvantk_run.sh`: it sets up node-local scratch, binds it and the work root
into the container, caps Spark at the allocation, and runs the command inside the
image. Use the wrapper for an interactive `srun` shell, the template for jobs.

## 4. Scale-up

- Bigger input means more memory on one node, not more nodes. Go to the fat partition
  (`$HVANTK_FAT_PARTITION`) with `--exclusive`, all cores and all memory; the template
  keeps the Spark driver at 80 percent of the request so the kernel's memory cgroup
  does not kill it.
- Time: ask for the time the job needs plus half. A time-limit kill loses everything;
  Hail tables are not written incrementally.
- Node-local scratch holds Spark spill and Hail temporary files. Check its size once
  (`df -h "$SLURM_TMPDIR"` inside a job); a variant-scale build can spill hundreds of GB.
  If the node has no local scratch, set `HVANTK_SCRATCH_ROOT` to the fastest filesystem
  you have, never to home.
- Outputs on a parallel filesystem are directories of many part files. Expect it, do
  not copy them file by file across the network, and archive with `tar` when moving.
- Many independent builds: submit one array job over a task list rather than many
  jobs; `sbatch --array=1-N%K` with one line of the list per task. The
  dataset-build-at-scale skill covers this.
- Multi-node Spark through the cluster's own launcher is possible
  (`HVANTK_SKIP_HAIL_INIT` makes hvantk adopt an existing session) and almost never
  worth it: both the scheduler and Spark must agree on resources, and the failure modes
  are opaque. Stay on one fat node until the working set does not fit.

## 5. Traps

- `[Errno 111] Connection refused` at Hail start inside the container. Cause: Spark has
  no writable node-local scratch, and the real message, `DiskBlockManager: Failed to
  create any local dir`, is one line up. Fix: run through the template or the wrapper;
  they set `SPARK_LOCAL_DIRS`.
- The job ends with no Python traceback, `sacct` says `OUT_OF_MEMORY` or the log just
  stops. Cause: the Spark driver was allowed the whole request and the cgroup killed
  it. Fix: the template caps the driver at 80 percent; if you set it yourself, keep it
  there and raise `--mem-gb` instead.
- `Executor heartbeat timed out` on a shared node. Cause: Spark took every core on the
  node, not the allocation. Fix: the template pins `local[<cpus>]`; do not override it
  with `local[*]` unless the job is `--exclusive`.
- Hail warns about the Java version or dies with `UnsupportedClassVersionError`.
  Cause: a host Java module leaked into the environment. Fix: never `module load`
  anything around the container; the image carries Java 11.
- `seff` or `MaxRSS` reports far more memory than requested for a job that completed.
  Cause: some accounting plugins sum per-process RSS including shared pages. Fix:
  trust the request and the cgroup; size from failures, not from `MaxRSS`.
- `hvantk reprocess` fails inside the job with a download or HTTP error. Cause: the
  compute node has no internet. Fix: stage the raw file on the login node, then
  `--skip-download --no-check-drift`.
- A write fails with "read-only file system" on a path that is writable on the login
  node. Cause: that mount is read-only on compute nodes. Fix: write under
  `$HVANTK_WORK`; copy into the read-only tree from the login node afterwards.
- A pipeline such as `zcat f | head` exits 141 under `set -o pipefail`. Cause:
  SIGPIPE when the reader closes early. Fix: `{ zcat f || true; } | head`.
- An exported VCF has comma decimals (`AF=0,25`). Cause: the JVM took a locale with a
  comma decimal separator. Fix: the template exports `JAVA_TOOL_OPTIONS` with an
  English locale; keep it when you run the container by hand.
- The job sits `PENDING` for hours. Cause: the partition's queue, not your job. Fix:
  `squeue -p "$HVANTK_CPU_PARTITION" --start -j <id>` for the estimate; many sites have a
  short-job partition for anything under two hours, and its name is a site value: ask,
  do not guess.
- A wait loop that greps a log for a string never ends. Cause: the string went to
  another stream, or the job died. Fix: wait on the job id (section 3 step 4).
- Binding a parent directory into the container shows empty directories inside.
  Cause: nested automounted filesystems. Fix: bind the exact subpaths; `HVANTK_BIND`
  takes a comma-separated list.
- `--mail-user` in a `#SBATCH` line. Cause: Slurm does not expand variables in
  directives, so people write the address in. Fix: the submit script passes
  `$SLURM_EMAIL_NOTIFICATION` on the command line when it is set; never write an
  address into a script.
- The cluster clone is behind the commit you reviewed. Cause: pushing to one remote
  does not update the other. Fix: preflight step 2, every time.

## 6. Verification

1. `sacct -j <id> -X -n -o State` prints `COMPLETED`.
2. For a dataset build: the output directory exists and carries `_SUCCESS`, and the
   sidecar `<output>.provenance.json` exists and names the `schema_id` you expected.
3. Describe the table inside the image and write the receipt:
   ```bash
   bash skills/hpc-hail-job/templates/submit-hail-job.sh --name <run>-describe \
     --cpus 2 --mem-gb 8 --time 0-00:15:00 --wait -- \
     python skills/hpc-hail-job/templates/hail-describe.py "$HVANTK_WORK/out/hgnc.ht" \
       --receipt "$HVANTK_WORK/out/hgnc.receipt.json"
   ```
   The log shows the key, the row count and the fields; the receipt carries the row
   count, the key and the SHA-256 of the provenance sidecar, so a verdict about this
   output is tied to exactly this output.
4. Put `commit:`, `container:` and every `job:` line in the run notes next to the
   output. A result without them cannot be reproduced.
