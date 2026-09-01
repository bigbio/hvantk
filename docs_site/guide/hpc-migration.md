# Migrating hvantk to an HPC cluster

This guide describes how to move the hvantk toolkit, its environment, and its data
from a local workstation to a shared HPC cluster, so that dataset builds, pipelines,
and experiments run at scale on the cluster.

It is written for a cluster with:

- a **Slurm** scheduler (with a supported Spark-on-Slurm launcher),
- **Apptainer/Singularity** available for containers,
- **outbound internet on compute nodes** (if your compute nodes are firewalled, see
  the staging notes in §6 below),
- a standard three-tier filesystem (small backed-up **home**, large **project/work**,
  fast purged **scratch**).

Adjust the partition names, paths, and Spark-on-Slurm launcher to your site.

---

## 0. TL;DR — the recommended path

1. **Code** → `git clone` the repo on the cluster. The gitignored `local/` data tree
   moves separately (step 4).
2. **Environment** → ship **one Apptainer image** (`hvantk.sif`) built from the
   repo's `poetry.lock`, pinning **Java 11** + Hail 0.2.137 + pyspark 3.5.8 + all
   native libraries. This is the most important decision: it eliminates the
   native-library portability class of failure (e.g. an `RPATH`/`libz` break in a
   compiled wheel like `pysam`/`pyBigWig`).
3. **Hail** → default to **Spark local mode on one exclusive fat node**; escalate to
   the cluster's **Spark-on-Slurm** multi-node launcher (via `HVANTK_SKIP_HAIL_INIT`)
   only when the working set exceeds a fat node's RAM.
4. **Data** → move the data that must be transferred via **Globus** (through a Data
   Transfer Node), **re-fetch** the public bulk downloads on the login node,
   **regenerate** the derived Hail `.mt`s on-cluster, and **skip** redundant/raw
   files.
5. **Streaming experiments** → where compute nodes have internet, code that streams
   from VEP / GCS / UCSC / eQTL sources runs as-is. Still **stage the big/repeated
   pulls** (gnomAD, the ~9 GB phyloP bigWig) to project storage for speed and to
   respect rate limits.
6. **Validate** every new node with `hvantk utils check-install` before running
   pipelines.

---

## 1. What moves, and how

| Layer | What | Mechanism | Notes |
|---|---|---|---|
| Code | `hvantk/`, `pyproject.toml`, `poetry.lock`, `tests/` | `git clone` | tracked; tiny |
| Environment | Python 3.10 + Hail/JVM + native libs | **Apptainer `.sif`** | built from `poetry.lock`; §3 |
| Data | gitignored `local/` tree | **Globus / refetch / regenerate** | §5 split |
| Secrets/access | EGA gvcfs, UKB-PPP, other gated sources | manual copy | §5 |

**Why a container and not conda or `poetry install` on the node.** A conda env or a
node-level `poetry install` still links against the *host's* system libraries
(`libz`, `libbz2`, `libcurl`, `libhdf5`, `glibc`). A mismatch there is the canonical
cause of `dlopen`/`RPATH` failures in compiled wheels. A `.sif` carries its own
copies of those libraries, so the dynamic loader resolves against them regardless of
the node. A container is also **one file** (vs. a conda env's ~100k files, which
exhausts inode quotas and stresses the shared filesystem). hvantk's stack is
unusually native-lib-heavy — Hail's JVM/Spark, `pysam` (htslib), `pyarrow` (Arrow
C++), `h5py` (HDF5), `scipy`/`scikit-learn` (BLAS/LAPACK), `duckdb` — so the
container payoff is large.

### 1.1 Keeping the cluster's clone current

The cluster has its own bare repo, and the working repo pushes to it via a **second
remote** called `hpc`. Local branches track `origin` (GitHub), so a bare `git push`
**does not** update the cluster:

| Remote | URL | Role |
|---|---|---|
| `origin` | `git@github.com:bigbio/hvantk.git` | GitHub; PRs and CI |
| `hpc` | `hpc:/user/<user>/git/pyvatk.git` | bare repo on the cluster; the node-side clone pulls from here |

After anything lands on `main` or `dev`:

```bash
# confirm it is a fast-forward, so no cluster-side commits are clobbered
git fetch hpc
git merge-base --is-ancestor hpc/dev origin/dev && echo "fast-forward, safe"

git push hpc main dev
```

**This matters more than ordinary repo hygiene**, because the `.sif` is built from
`poetry.lock` (§3.2). A stale cluster clone rebuilds the *previous* dependency set
without any error — the build succeeds, it is just the wrong environment. Rebuild the
image whenever `poetry.lock` changes, not only when `hvantk/` does.

---

## 2. Cluster facts to confirm on first login

Run on the login node to fill in the exact paths/flags used below:

```bash
sinfo -o "%P %l %c %m %G"          # partitions: time limit, cores, mem, gres (lscratch?)
module avail 2>&1 | grep -iE "java|apptainer|singularity|spark|globus"
echo $HOME; df -h $HOME            # home quota (code only)
# project/work + scratch paths & quotas are site-specific (check docs / `quota`)
which apptainer || which singularity
apptainer build --fakeroot --help 2>&1 | head -1   # is --fakeroot allowed?
java -version                      # default Java — MUST be 8 or 11 for Hail
```

Note the **scratch path** (often `/scratch/$USER` or `/lscratch/$SLURM_JOB_ID`),
whether **`--fakeroot`** is permitted (if not, build the image elsewhere and copy
it), and the **Spark-on-Slurm** launcher's name/usage.

---

## 3. Environment: build the Apptainer image

### 3.1 Critical version pins (do not drift)

The repo's `poetry.lock` encodes a tested triad: Hail **0.2.137** + pyspark
**3.5.8** + py4j 0.10.9.9 + **numpy 2.2.5** (with the `np.bool = np.bool_` shim in
`hvantk/core/utils/hail_context.py` that Hail 0.2.x needs under NumPy 2). Build from
the lock; do not regenerate or loosen it. And:

> **Java 11, not 17.** Hail 0.2.137 / Spark 3.5 (Scala 2.12) support **Java 8 or 11
> only**. HPC `module load java` frequently defaults to 17/21, which breaks Hail.
> Bake Java 11 into the container so the host's Java module is irrelevant.

### 3.2 `hvantk.def` (build from the locked environment)

The definition file lives in the repo at **`containers/hvantk.def`**
— build from that file rather than copying a snippet, so the base image and extras
cannot drift from what was last built and validated. It installs the **exact** locked
dependency set, not an unpinned `pip install hail`.

Two choices in it are load-bearing and must not be "modernised" without re-validating:

- **Base image is `python:3.10-slim-bullseye`, not bookworm.** Debian 12 (bookworm) has
  no `openjdk-11` package at all — `apt-get install openjdk-11-jdk-headless` fails with
  *"Package 'openjdk-11-jdk-headless' has no installation candidate … the following
  packages replace it: openjdk-17-jre-headless"*. Hail 0.2.x / Spark 3.5 support Java 8
  or 11 **only**, so a bookworm base either fails the build or silently gives you
  Java 17. Debian 11 (bullseye) ships `openjdk-11-jdk-headless` (11.0.32.1).
- **Extras include `expression`.** The set
  `hgc ptm ancestry psroc constraint enrichex cohort viz duckdb expression`
  is checked against `[project.optional-dependencies]` in `pyproject.toml` and covers
  every unique package across all extras. Omitting `expression` builds an image with no
  scanpy, so `hvantk expression …` cannot run. `ml` and `interactive` are redundant:
  scikit-learn arrives via `ancestry`/`psroc`, plotly via `viz`.

The `%post` block ends by printing `java -version` and importing `hail`, so a broken
image fails at **build** time rather than at first use.

### 3.3 Build (rootless) and validate

```bash
# Keep build cache/tmp OFF the small home quota
export SINGULARITY_CACHEDIR=$WORK/containers/cache     # APPTAINER_* if you have apptainer
export SINGULARITY_TMPDIR=$WORK/containers/tmp
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

# Build from the REPO ROOT: %files paths in the def are relative to the build CWD.
cd <the cluster clone>
singularity build $WORK/containers/hvantk.sif containers/hvantk.def
```

**`--fakeroot` is not usable on every cluster.** It needs a `/etc/subuid` entry for your
account; without one the build fails immediately with
`could not use fakeroot: no valid mapping entry found for <user>`, and a plain
unprivileged build is refused outright with
`--remote, --fakeroot, or the proot command are required to build this source as a
non-root user`. Check before assuming:

```bash
grep "^$USER:" /etc/subuid    # empty output => --fakeroot will NOT work here
```

When there is no `subuid` entry, SingularityCE 3.11+/4.x accepts a static **`proot`**
instead, which needs no privileges at all:

```bash
mkdir -p ~/bin && curl -fsSL -o ~/bin/proot https://proot.gitlab.io/proot/bin/proot
chmod +x ~/bin/proot
export PATH="$HOME/bin:$PATH"       # singularity picks proot up from PATH
singularity build $WORK/containers/hvantk.sif containers/hvantk.def
```

Builds take roughly 10 minutes and the image is ~2.2 GB.

Validate the full native + JVM stack on a **compute** node (not login) — and go through
the run wrapper, because Hail cannot initialise inside the container without Spark
scratch (see §4.1):

```bash
srun -c 4 --mem 16g bash containers/hvantk_run.sh utils check-install
#   expects: Hail version prints, balding_nichols_model smoke test passes.
```

`hvantk utils check-install` is the canonical go/no-go for a node — it initializes
Hail, prints `hl.version()`, runs a Hail smoke test, and diagnoses proxy problems.

## 4. Running Hail on the cluster

### 4.1 Default: Spark local mode on one exclusive fat node

hvantk runs Hail in **local Spark mode** — every `init_hail()` call uses local
defaults (no master/memory config in the repo). On HPC that maps to **one exclusive
node, many cores, high memory**, with Spark temp on **node-local scratch**.

> **`SPARK_LOCAL_DIRS` is mandatory when running from the container, and its absence
> is misdiagnosed.** Without a writable node-local scratch bound into the image, Hail
> init dies with
> `DiskBlockManager: ERROR: Failed to create any local dir`, followed by
> `Hail initialisation failed: [Errno 111] Connection refused` from py4j. The second
> message is what you see first and it reads like a network or proxy fault; it is not.
> `containers/hvantk_run.sh` sets the scratch dir up and binds it, so prefer:
>
> ```bash
> bash containers/hvantk_run.sh utils check-install
> HVANTK_BIND="$WORK:$WORK" bash containers/hvantk_run.sh reprocess clinvar:variants ...
> ```

```bash
#!/bin/bash
#SBATCH --job-name=hvantk-hail
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out

# node-local fast scratch for Spark spill + Hail tmp (NEVER home/Lustre/GPFS)
export SCRATCH_DIR=${SLURM_TMPDIR:-/scratch/$USER/$SLURM_JOB_ID}
mkdir -p "$SCRATCH_DIR"
export SPARK_LOCAL_DIRS="$SCRATCH_DIR"
export NO_PROXY=localhost,127.0.0.1,0.0.0.0,::1   # protect py4j localhost socket

apptainer exec \
  --bind /project/$USER/hvantk:/data \
  --bind "$SCRATCH_DIR":/scratch \
  hvantk.sif \
  hvantk hgc pipeline -i /data/gvcfs -o /data/out/cohort.vds
```

The three settings that prevent the most common failures:

1. **Driver memory must sit *below* the Slurm `--mem` cgroup limit.** In local mode
   the driver JVM *is* the executor, and Hail uses large off-heap/native memory
   (BLAS, LZ4) plus Python workers that live outside the JVM heap. Setting
   `spark.driver.memory == --mem` triggers a silent **cgroup OOM-kill**. Use
   **~75–85% of `--mem`** (e.g. `--mem=240G` → `spark.driver.memory=192g`).
2. **Pin cores; don't let Spark grab the machine.** On an `--exclusive` node
   `local[*]` is fine; on a shared node use `local[<cpus-per-task>]`. Uncapped
   `local[*]` on a shared node causes thread oversubscription and "Executor
   heartbeat timed out" errors.
3. **Spark/Hail temp on node-local scratch**, never the parallel filesystem.

To set these explicitly, wrap your script's `init_hail()` → `hl.init()`. A drop-in
Hail init for a 240 GB node:

```python
import os, hail as hl
scratch = os.environ["SCRATCH_DIR"]
hl.init(
    master="local[*]",                 # exclusive node → all cores OK
    tmp_dir=scratch,
    spark_conf={
        "spark.driver.memory": "192g",  # ~80% of --mem=240G
        "spark.local.dir": scratch,
        "spark.driver.maxResultSize": "16g",
        "spark.serializer": "org.apache.spark.serializer.KryoSerializer",
        "spark.kryo.registrator": "is.hail.kryo.HailKryoRegistrator",
    },
)
```

**Avoid small-file storms.** Each Hail `.ht`/`.mt` partition is a part file; too many
tiny parts saturate Lustre/GPFS metadata servers and slow the whole cluster.
Coalesce before writing big outputs, and prefer scratch then archive:

```python
mt = mt.naive_coalesce(64)        # ~hundreds-of-MB parts, not thousands of tiny ones
mt.write("/scratch/.../cohort.mt", overwrite=True)
# then tar the .mt directory to project storage rather than leaving 1000s of parts on Lustre
```

### 4.2 Escalation: multi-node via the cluster's Spark-on-Slurm launcher

Use this only when one fat node cannot hold the working set. hvantk supports
**adopting an externally-initialized Spark/Hail session** via the
**`HVANTK_SKIP_HAIL_INIT`** environment variable (checked in
`hvantk/core/utils/hail_context.py`): set it non-empty and `init_hail()` will not
call `hl.init()`, deferring to a session you created. Pattern:

1. Launch the cluster's Spark standalone cluster inside a Slurm allocation (per your
   center's docs — e.g. a `spark start ...` that yields a `spark://host:7077`
   master).
2. In a thin wrapper, call `hl.init(sc=<existing SparkContext>)` (or
   `hl.init(master="spark://host:7077", ...)`) **before** importing hvantk
   pipelines, then `export HVANTK_SKIP_HAIL_INIT=1` so hvantk adopts that backend.
3. Run inside the container with the Hail jar on the Spark classpath (`--jars
   $HAIL_HOME/hail-all-spark.jar`, Kryo serializer + `HailKryoRegistrator`).

Spark-on-Slurm launchers are cluster-specific and double-allocate resources (Slurm
*and* spark-submit flags must agree); expect operational friction. Stay in §4.1 local
mode unless data forces this.

---

## 5. Migrating the data

The gitignored `local/` tree is **not** part of the git clone — it must be moved
explicitly. Split it by action rather than copying it wholesale.

### 5.1 Action split

| Action | Examples | How |
|---|---|---|
| **Transfer** (experiment state / gated / source-specific) | notebooks + outputs, single-cell `.h5ad` atlases, planning docs, EGA gvcfs, interactome BED, GTEx proteomics, UKB-PPP, per-study cardiac data | **Globus** via DTN |
| **Refetch** on login node | raw UCSC Cell Browser dumps, GTEx/iPSC-CM eQTL stats, ArrayExpress, 1000 Genomes, GWAS Catalog, GenCC, ClinVar, Ensembl GTF, gnomAD lof-metrics | `hvantk download …` / skill downloaders |
| **Regenerate** on-cluster | derived Hail `.mt`s (e.g. Expression Atlas, GTEx TPM matrices) | hvantk build steps |
| **Skip** (redundant) | uncompressed `.tsv` where a `.bgz` exists; raw dumps that a built `.h5ad` supersedes; duplicate `.gz` next to `.bgz` | don't copy |

Net effect: physically move the experiment state and gated sources, refetch public
bulk on the cluster, regenerate derived artifacts, and avoid copying redundancy.
Globus is fast enough that transferring what you already have is often simpler than
re-deriving it — decide per-subtree.

### 5.2 Filesystem placement

- **`$HOME`** (small, backed up): the cloned repo + `hvantk.sif` + sbatch scripts.
  Never bulk data.
- **`/project` (or work)**: the `local/` data tree, e.g. `/project/$USER/hvantk/local/`,
  **bound into the container** as `/data`. Keep the same `local/...` substructure so
  experiment scripts' relative paths work.
- **scratch (purged)**: Spark/Hail tmp, per-job intermediates, container build cache.
  Copy results back to project before the purge window (often 30–60 days).

### 5.3 Moving the data

```bash
# Preferred: Globus (resumable, parallel, auto-checksum), endpoint = the cluster DTN
module load globus-cli 2>/dev/null || true
globus transfer --recursive --verify-checksum \
  "$LAPTOP_ENDPOINT:/path/to/pyvatk/local" \
  "$CLUSTER_ENDPOINT:/project/$USER/hvantk/local"

# Fallback: rsync through the Data Transfer Node (NOT the login node)
rsync -avP --partial \
  /path/to/pyvatk/local/ \
  $USER@dtn.cluster.edu:/project/$USER/hvantk/local/
# verify (checksum compare, dry-run)
rsync -avnc /path/to/pyvatk/local/ \
  $USER@dtn.cluster.edu:/project/$USER/hvantk/local/
```

Add `--exclude` rules for the §5.1 "skip" rows (e.g. uncompressed `.tsv` files that
duplicate a `.bgz`) to avoid moving redundant gigabytes.

---

## 6. External / streaming data

Several experiment paths stream from external sources at runtime. Where compute nodes
have internet they run in-job; the concern is volume and rate limits, so stage the
heavy/repeated pulls:

| Source | Used by | Recommendation |
|---|---|---|
| gnomAD v4.1 (GCS / googleapis) | gnomAD-MAPS / constraint experiments | **Stage once** to `/project/.../gnomad/` (per-chr sites or a filtered HT); repoint scripts at the local path |
| UCSC phyloP100way bigWig (~9 GB) | phyloP per-codon queries | **Download once** to project; point the query at the local `.bw` (also removes remote-query latency) |
| Ensembl VEP REST | VEP annotation experiments | keep batch ≤ 200 + honor `Retry-After`; for big sets use a local VEP install. (Public VEP does **not** serve AlphaMissense — use the standalone TSV if needed) |
| eQTL Catalogue / Pan-UKB (EBI FTP, AWS S3) | coloc / QTL experiments | transfer any existing local `data_cache` so re-streaming is unnecessary |
| Ensembl GTF (ftp.ensembl.org) | PTM pipeline | transfer the cached GTF; pass the local path to skip the download step |

Two operational notes:

- **Drift checks need internet by design.** Run `hvantk drift` / `hvantk download`
  on the login node; pass `--no-check-drift` to `hvantk reprocess` inside batch jobs
  if a compute partition is ever firewalled.
- **Proxy:** if the cluster routes outbound traffic through an HTTP proxy, set
  `http_proxy`/`https_proxy` **and** keep `NO_PROXY=localhost,127.0.0.1,0.0.0.0,::1`
  so Hail's py4j localhost socket is not intercepted.

---

## 7. Running experiments

### 7.1 A `reprocess` build (single node)

```bash
srun --partition=standard -c 8 --mem 48g -t 4:00:00 \
  apptainer exec --bind /project/$USER/hvantk/local:/data hvantk.sif \
    hvantk reprocess clinvar:variants \
      --raw-dir /data/data/clinvar_ts --output /data/out/clinvar.ht
```

### 7.2 Embarrassingly-parallel fan-out (array job)

Ideal for per-gene/per-site sweeps. Build a task list and index it with
`$SLURM_ARRAY_TASK_ID`:

```bash
#!/bin/bash
#SBATCH --job-name=hvantk-fanout
#SBATCH --array=1-500%50            # 500 genes, ≤50 concurrent
#SBATCH --cpus-per-task=4 --mem=16G --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
GENE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" genes.txt)
apptainer exec --bind /project/$USER/hvantk/local:/data hvantk.sif \
  hvantk qtlcascade cascade --gene "$GENE" --out /data/out/${GENE}.ht
```

For experiments that stream remote data (e.g. phyloP bigWig, gnomAD), pre-stage those
inputs to project storage (§6) and repoint the scripts at the local paths before
launching at scale.

---

## 8. Reproducibility and provenance

- **Pin the container, not just the code.** Record the `.sif`'s build date and the
  `poetry.lock` hash; a given `hvantk.sif` reproduces the exact Hail/numpy/native
  stack. Re-derive the image only by re-running `apptainer build` against the same
  lock.
- **hvantk stamps provenance** on every artifact (`ctx.provenance(...)`) and has
  `hvantk drift` for source fingerprints — keep using both on-cluster so results
  carry their input lineage.
- Treat scratch as ephemeral; the **project copy of outputs + the `.sif` + the git
  SHA** are the reproducible record.

---

## 9. Runbook (ordered)

1. **Login node:** `git clone` the repo into `$HOME`. Note partitions, scratch path,
   `--fakeroot` availability, Spark-on-Slurm launcher (§2).
2. **Build the image:** write `hvantk.def` (§3.2), set `APPTAINER_*` to scratch,
   `apptainer build --fakeroot hvantk.sif hvantk.def` (or build on a laptop + copy).
3. **Validate:** `srun ... apptainer exec hvantk.sif hvantk utils check-install`.
4. **Move data (§5):** Globus the "transfer" set into `/project/$USER/hvantk/local/`;
   refetch the public bulk on the login node; regenerate `.mt`s on-cluster; skip the
   redundant files.
5. **Stage streaming inputs (§6):** download the phyloP bigWig + a gnomAD slice to
   project; repoint the relevant scripts at local paths.
6. **First real run:** submit the §4.1 single-fat-node Hail job for a known pipeline;
   confirm it writes to project and stays under the cgroup memory limit.
7. **Scale out:** use array jobs (§7.2) for fan-out; reserve the Spark-on-Slurm path
   (§4.2) for working sets that exceed a fat node.

---

## 10. Gotchas quick-reference

- **Java 17 default** → Hail breaks. Force **Java 11** (baked into the container).
- **A bookworm base image** → there is no `openjdk-11` in Debian 12 at all; apt offers
  `openjdk-17-jre-headless` instead. Use `python:3.10-slim-bullseye`.
- **`--fakeroot` with no `/etc/subuid` entry** → `no valid mapping entry found`, and an
  unprivileged build is refused. Put a static `proot` on `PATH` instead (§3.3).
- **Container Hail init: `[Errno 111] Connection refused`** → almost never the network.
  Look one line up for `DiskBlockManager: Failed to create any local dir`: Spark has no
  writable scratch. Use `containers/hvantk_run.sh` (§4.1).
- **`spark.driver.memory == --mem`** → silent cgroup OOM-kill. Use ~80%.
- **`local[*]` on a shared node** → heartbeat timeouts. Pin to `--cpus-per-task`, or
  use `--exclusive`.
- **Spark/Hail tmp on Lustre/GPFS or `$HOME`** → slow + metadata storms. Use
  node-local scratch.
- **Over-partitioned `.mt`/`.ht`** → millions of tiny part files choke the parallel
  filesystem. `naive_coalesce` before writing; archive off the parallel FS.
- **HTTP proxy intercepting localhost** → Hail/py4j init fails. Set `NO_PROXY` for
  localhost.
- **Stray system/conda Python** below the 3.10 floor → use the container's Python;
  never run hvantk against an unmanaged interpreter.
- **Core install ≠ full toolkit** — install the right **extras** (`hgc ptm
  ancestry psroc constraint enrichex cohort viz duckdb expression`) in the image.
  Dropping `expression` yields an image with no scanpy, so `hvantk expression …` fails.
- **Scratch is purged** — copy results to project/home before the window.
