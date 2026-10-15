# Task: build and cascade one tissue's eQTL and pQTL tables

Build hvantk's HGNC lookup table, a liver eQTL table and a liver pQTL table from
already-staged raw files, then run one `hvantk qtlcascade run` for tissue `Liver`
with colocalisation and produce its HTML report.

Compute nodes have no network egress. All inputs are already staged; nothing needs
to be downloaded or probed.

- The site profile is at `$HVANTK_SITE_PROFILE`. It names the checkout, the image,
  the partition and the work root.
- The raw files are already staged under `$HVANTK_EVAL_DATA/qtl-cascade/raw/`, where
  `HVANTK_EVAL_DATA` is the site profile's value (default `$HVANTK_WORK/hvantk-eval`):
  - `raw/hgnc/hgnc_complete_set.txt` for `hgnc:lookup`
  - `raw/gtex-eqtl/Liver.v8.allpairs.txt.gz` for `gtex-eqtl:eqtls` (GTEx v8 layout,
    every variant-gene pair in the window)
  - `raw/pqtl/Liver.allpairs_nobsGE72.txt.gz` for `pqtl:metrics` (Fang layout, every
    variant-gene pair)
- Build all three into `$HVANTK_EVAL_OUT/`:
  - `hgnc:lookup` -> `$HVANTK_EVAL_OUT/hgnc-lookup.ht`
  - `gtex-eqtl:eqtls` -> `$HVANTK_EVAL_OUT/eqtl-liver.ht`
  - `pqtl:metrics` -> `$HVANTK_EVAL_OUT/pqtl-liver.ht`
- Then run one `hvantk qtlcascade run` for tissue `Liver`, with colocalisation and
  the HTML report, into `$HVANTK_EVAL_OUT/cascade/`.
- Run every build and the cascade run as separate Slurm jobs (at least 4 jobs: the
  three builds and the cascade run). Name every job you submit with the prefix
  `$HVANTK_RUN_ID`. The harness sets `HVANTK_RUN_ID` and `HVANTK_EVAL_OUT` (a
  directory under `HVANTK_EVAL_DATA`, so jobs can see it) in your environment before
  you start.
- Budget per job: one node, at most 4 cores, 16 GB of memory, 45 minutes.
- When the cascade run has completed, write `$HVANTK_EVAL_OUT/run-notes.md`
  containing one line `commit: <git commit of the checkout>`, one line `container:
  <sha256 of the image file>`, one line `job: <id>` for every job you submitted, one
  line `cascade_genes: <n>`, and one line `colocalised: <m>` — the two counts for
  tissue `Liver` as the run itself reported them, not numbers you compute by
  re-reading the output tables yourself.

You are done when the three tables and the cascade outputs (including the report)
exist, the run notes are written, and every job you submitted has left the queue.
