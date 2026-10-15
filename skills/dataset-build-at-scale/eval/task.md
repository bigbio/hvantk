# Task: build three hvantk datasets as one array job

Build three of hvantk's small gene-level datasets on this cluster, inside the hvantk
container image, as a single Slurm array job with one task per dataset. Do not run any
build on the login node, and do not download anything: the raw files are already
staged.

- The site profile is at `$HVANTK_SITE_PROFILE`. It names the checkout, the image, the
  partition and the work root.
- The raw files are already staged under `$HVANTK_EVAL_DATA/dataset-build-at-scale/raw/`,
  where `HVANTK_EVAL_DATA` is the site profile's value (default
  `$HVANTK_WORK/hvantk-eval`):
  - `raw/hgnc/hgnc_complete_set.txt` for `hgnc:lookup`
  - `raw/clingen/` (one dated `Clingen-Gene-Disease-Summary-<date>.csv`) for
    `clingen:gene-disease`
  - `raw/gencc/` (one dated `gencc-submissions-<date>.tsv`) for `gencc:submissions`
- Build all three into `$HVANTK_EVAL_OUT/`:
  - `hgnc:lookup` → `$HVANTK_EVAL_OUT/hgnc-lookup.ht`
  - `clingen:gene-disease` → `$HVANTK_EVAL_OUT/clingen-gene-disease.ht`
  - `gencc:submissions` → `$HVANTK_EVAL_OUT/gencc-submissions.ht`
- Run the three builds as **one Slurm array job** (one task per dataset), not as three
  separate job submissions.
- Name every Slurm job you submit with the prefix `$HVANTK_RUN_ID`. The harness sets
  `HVANTK_RUN_ID` and `HVANTK_EVAL_OUT` (a directory under `HVANTK_EVAL_DATA`, so jobs
  can see it) in your environment before you start.
- Budget per array task: one node, at most 4 cores, 16 GB of memory, 30 minutes, at most
  three tasks running at once.
- When every build has completed, write `$HVANTK_EVAL_OUT/run-notes.md` containing one
  line `commit: <git commit of the checkout>`, one line `container: <sha256 of the
  image file>`, and one line `job: <id>` for every job you submitted (the array
  submission is one job, however many tasks it fans out to).

You are done when all three datasets exist, the run notes are written, and every job
you submitted has left the queue.
