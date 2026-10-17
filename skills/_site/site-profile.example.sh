# hvantk site profile (template). Copy to site-profile.local.sh, fill in, and
#   export HVANTK_SITE_PROFILE=/path/to/site-profile.local.sh
# The filled copy is gitignored: it names your cluster, keep it out of any repository.
# Every value is empty here on purpose; skills stop and ask rather than guess.
# See SITE-PROFILE.md for what each one means and how to find it.

# --- required ---------------------------------------------------------------
HVANTK_REPO=""              # path of the hvantk checkout on the cluster
HVANTK_SIF=""               # path of the container image (a symlink is fine)
HVANTK_CPU_PARTITION=""     # partition for single-node CPU jobs
HVANTK_NODE_CPUS=""         # cores on one node of that partition, e.g. 32
HVANTK_NODE_MEM_GB=""       # memory on one node of that partition, e.g. 240
HVANTK_WORK=""              # project/work root, writable from compute nodes

# --- optional ---------------------------------------------------------------
HVANTK_SIF_SHA256=""        # cached digest of the image; skills recompute and verify it
HVANTK_FAT_PARTITION=""     # partition for exclusive high-memory jobs; defaults to CPU
HVANTK_ACCOUNT=""           # Slurm account, if your site requires one
HVANTK_QOS=""               # Slurm QOS, if your site requires one
HVANTK_MAX_TIME=""          # partition time limit as D-HH:MM:SS, e.g. 7-00:00:00
HVANTK_CONTAINER_BIN=""     # singularity or apptainer; default singularity
HVANTK_BIND=""              # extra container bind specs, comma-separated
HVANTK_SCRATCH_ROOT=""      # node-local scratch root; default $SLURM_TMPDIR or /tmp
HVANTK_LOG_DIR=""           # job logs; default $HVANTK_WORK/hvantk-logs
HVANTK_READONLY_MOUNTS=""   # paths mounted read-only on compute nodes, comma-separated
HVANTK_EGRESS=""            # open or blocked: can compute nodes reach the internet?
HVANTK_EVAL_DATA=""         # where eval slices are staged; default $HVANTK_WORK/hvantk-eval
HVANTK_CLUSTER_REMOTE=""    # name of the cluster remote in the clone you push from, e.g. hpc
