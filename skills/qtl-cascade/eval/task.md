# Task: build and cascade one tissue's eQTL and pQTL tables

Build hvantk's HGNC lookup table, a liver eQTL table and a liver pQTL table from
already-staged raw files, then run the QTL cascade with colocalisation for tissue
`Liver` and produce its HTML report.

There is no network access from this environment: nothing can be downloaded or
probed; the staged files are complete.

- The site profile is at `$HVANTK_SITE_PROFILE`. It names the checkout, the image,
  the partition and the work root.
- The raw files are already staged under `$HVANTK_EVAL_DATA/qtl-cascade/raw/`, where
  `HVANTK_EVAL_DATA` is the site profile's value (default `$HVANTK_WORK/hvantk-eval`):
  - `raw/hgnc/hgnc_complete_set.txt` for `hgnc:lookup`
  - `raw/gtex-eqtl/Liver.v8.allpairs.txt.gz` for `gtex-eqtl:eqtls` (GTEx v8 layout,
    every variant-gene pair in the window — this is already an all-pairs file; do not
    filter it further at build time)
  - `raw/pqtl/Liver.allpairs_nobsGE72.txt.gz` for `pqtl:metrics` (Fang layout, every
    pair — also already all-pairs)
- Build all three into `$HVANTK_EVAL_OUT/`:
  - `hgnc:lookup` -> `$HVANTK_EVAL_OUT/hgnc-lookup.ht`
  - `gtex-eqtl:eqtls` (`source=gtex_v8`, `tissue=Liver`, `p_threshold=0`) ->
    `$HVANTK_EVAL_OUT/eqtl-liver.ht`
  - `pqtl:metrics` (`tissue=Liver`, `hgnc_ht=` the HGNC table above) ->
    `$HVANTK_EVAL_OUT/pqtl-liver.ht`
- Then run the QTL cascade for tissue `Liver`, with colocalisation, the default
  p-value thresholds and the default window, and the HTML report, into
  `$HVANTK_EVAL_OUT/cascade/`. Use the eQTL and pQTL tables above as both the cascade
  inputs and the coloc (allpairs) inputs.
- Run every build and the cascade run as separate Slurm jobs (at least 4 jobs: the
  three builds and the cascade run). Name every job you submit with the prefix
  `$HVANTK_RUN_ID`. The harness sets `HVANTK_RUN_ID` and `HVANTK_EVAL_OUT` (a
  directory under `HVANTK_EVAL_DATA`, so jobs can see it) in your environment before
  you start.
- Budget per job: one node, at most 4 cores, 16 GB of memory, 45 minutes.
- When the cascade run has completed, write `$HVANTK_EVAL_OUT/run-notes.md`
  containing one line `commit: <git commit of the checkout>`, one line `container:
  <sha256 of the image file>`, one line `job: <id>` for every job you submitted, one
  line `cascade_genes: <n>`, and one line `colocalised: <m>` — the two counts the
  cascade run itself reported for tissue `Liver` (its own printed summary line), not
  numbers you compute by re-reading the output tables yourself.

You are done when the three tables and the cascade outputs (including the report)
exist, the run notes are written, and every job you submitted has left the queue.
