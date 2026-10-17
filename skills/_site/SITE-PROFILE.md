# The site profile

One file holds every cluster-specific value the workflow skills need. Fill it in on
the cluster, in the `$HVANTK_REPO` checkout: copy `site-profile.example.sh` to
`site-profile.local.sh`, fill it in, and export `HVANTK_SITE_PROFILE` pointing at it.
The filled copy is gitignored and must stay out of every repository: it names your
cluster.

Templates read values from it with `source "$HVANTK_SITE_PROFILE"`, and
`check-site-profile.sh` refuses to continue while a required value is empty. A skill
that meets an empty value stops and asks the user. It never guesses. Before the first
image exists, `HVANTK_SIF` is set but points nowhere; pass `check-site-profile.sh
--allow-missing-image` (the `hpc-container` templates already do), or run
`hvantk:hpc-container` to build the image first.

## Required

| variable | meaning | how to find it |
| --- | --- | --- |
| `HVANTK_REPO` | the hvantk checkout the jobs run from | where you cloned it on the cluster |
| `HVANTK_SIF` | the container image built from the repo's lock file | wherever you built or copied it; a symlink to the current image is fine |
| `HVANTK_CPU_PARTITION` | partition for single-node CPU jobs | `sinfo -s` lists partitions; your site docs say which is for general CPU work |
| `HVANTK_NODE_CPUS` | cores on one node of that partition | `sinfo -p <partition> -o "%c %m %l"` |
| `HVANTK_NODE_MEM_GB` | memory on one node of that partition, in GB | same command, second column, which Slurm prints in MiB: divide by 1024 and round down |
| `HVANTK_WORK` | large project or work storage, writable from compute nodes | your site docs; often a `$WORK` or `$PROJECT` variable |

## Optional

| variable | meaning | default |
| --- | --- | --- |
| `HVANTK_SIF_SHA256` | a cached digest of the image; skills recompute the digest from the resolved file and treat this as a cache that must match | recomputed when empty |
| `HVANTK_FAT_PARTITION` | partition for exclusive high-memory jobs | the CPU partition |
| `HVANTK_ACCOUNT`, `HVANTK_QOS` | passed to `sbatch` when set | unset |
| `HVANTK_MAX_TIME` | the partition's time limit, `D-HH:MM:SS` | unknown; `sinfo` shows it |
| `HVANTK_CONTAINER_BIN` | `singularity` or `apptainer` | `singularity` (apptainer installs that name too) |
| `HVANTK_BIND` | extra bind specs for the container, comma-separated | none; the work root is bound by the templates |
| `HVANTK_SCRATCH_ROOT` | node-local scratch root | `$SLURM_TMPDIR`, else `/tmp` |
| `HVANTK_LOG_DIR` | where job logs go | `$HVANTK_WORK/hvantk-logs` |
| `HVANTK_READONLY_MOUNTS` | paths that are read-only on compute nodes | none known |
| `HVANTK_EGRESS` | `open` or `blocked`: compute-node internet | unknown; treat as blocked |
| `HVANTK_EVAL_DATA` | where eval slices are staged | `$HVANTK_WORK/hvantk-eval` |
| `HVANTK_CLUSTER_REMOTE` | the name of the cluster remote in the clone you push from (`git remote -v`) | unset; the sync skill asks |

## Things that never go in a profile

- Your e-mail. Export `SLURM_EMAIL_NOTIFICATION=you@example.org` in your shell startup;
  the submit scripts pass it on the command line, never in a `#SBATCH` line and never
  in a file.
- Credentials of any kind.
- Anything that would also be in a committed file.

## Discovering values on a new cluster

```bash
sinfo -s                                   # partitions and their node counts
sinfo -p <partition> -o "%c %m %l %N"      # cores, MB, time limit, nodes
which singularity apptainer                # container runtime
df -h "$HOME" ; quota 2>/dev/null          # home is small: never put data there
srun -p <partition> -c 1 --mem 1G -t 2 bash -c 'touch /path/to/work/.probe && echo rw || echo ro'   # is the work root writable from a node?
srun -p <partition> -c 1 --mem 1G -t 2 curl -sI https://github.com | head -1   # egress?
```
