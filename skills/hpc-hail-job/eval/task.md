# Task: build the HGNC lookup dataset as a Slurm job

Build hvantk's `hgnc:lookup` dataset on this cluster, inside the hvantk container
image, as a Slurm job. Do not run the build on the login node.

- The site profile is at `$HVANTK_SITE_PROFILE`. It names the checkout, the image, the
  partition and the work root.
- The raw HGNC file is already staged at
  `$HVANTK_EVAL_DATA/hpc-hail-job/raw/hgnc/hgnc_complete_set.txt`, where
  `HVANTK_EVAL_DATA` is the site profile's value (default `$HVANTK_WORK/hvantk-eval`).
  Do not download anything; compute nodes have no internet access.
- Write the dataset to `$HVANTK_EVAL_OUT/hgnc.ht`.
- Name every Slurm job you submit with the prefix `$HVANTK_RUN_ID`.
- Budget per job: one node, at most 8 cores, 32 GB of memory, 30 minutes.
- When the build has completed, write `$HVANTK_EVAL_OUT/run-notes.md` containing one
  line `commit: <git commit of the checkout>`, one line `container: <sha256 of the image
  file>`, and one line `job: <id>` for every job you submitted.

You are done when the dataset exists, the run notes are written, and every job you
submitted has left the queue.
