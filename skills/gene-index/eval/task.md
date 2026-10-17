# Task: build the gene index

Build hvantk's gene index on this cluster, inside the hvantk container image, as Slurm
jobs. Do not run any build on the login node, and do not download anything: the raw
files are already staged.

There is no network access from this environment: nothing can be downloaded or probed;
the staged files are complete.

- The site profile is at `$HVANTK_SITE_PROFILE`. It names the checkout, the image, the
  partition and the work root.
- The raw files are already staged under `$HVANTK_EVAL_DATA/gene-index/raw/`, where
  `HVANTK_EVAL_DATA` is the site profile's value (default `$HVANTK_WORK/hvantk-eval`):
  - `raw/ensembl-gene/` (one pinned-release `.gtf.gz`) for `ensembl-gene:structure`
  - `raw/hgnc/hgnc_complete_set.txt` for `hgnc:lookup`
- Build, as Slurm jobs, inside the image:
  1. `ensembl-gene:structure` and `hgnc:lookup` from the staged raw directories. Both
     builds, with their `.provenance.json` sidecars, MUST be written somewhere under
     `$HVANTK_EVAL_OUT/` (the exact path and filename are yours to choose; the check
     finds each one by its sidecar's own `schema_id`, searching up to 3 directories
     below `$HVANTK_EVAL_OUT`).
  2. The protein-coding gene index (`hvantk annotate spine`, its default `--biotype`)
     from those two builds, into `$HVANTK_EVAL_OUT/gene-index.ht`.
- Name every Slurm job you submit with the prefix `$HVANTK_RUN_ID`. The harness sets
  `HVANTK_RUN_ID` and `HVANTK_EVAL_OUT` (a directory under `HVANTK_EVAL_DATA`, so jobs can
  see it) in your environment before you start.
- Budget per job: one node, at most 4 cores, 16 GB of memory, 30 minutes.
- When every build has completed, write `$HVANTK_EVAL_OUT/run-notes.md` containing:
  - one line `commit: <git commit of the checkout>`
  - one line `container: <sha256 of the image file>`
  - one line `job: <id>` for every job you submitted
  - one line `genes: <count>`, the gene count the gene-index build command printed
  - one line `mapping_rate: <rate>`, the HGNC mapping rate the gene-index build command
    printed (its own percentage form, e.g. `98.50%`, is fine)

You are done when the gene index exists, the run notes are written, and every job you
submitted has left the queue.
