# Task: infer the ancestry of held-out 1000 Genomes samples

A labelled reference and an unlabelled query were prepared from one chromosome of the
1000 Genomes high-coverage release. Run hvantk's ancestry pipeline on this cluster, as a
Slurm job inside the container image, and produce predictions for every query sample.

- The site profile is at `$HVANTK_SITE_PROFILE`.
- Inputs under `$HVANTK_EVAL_DATA/ancestry-inference/` (where `HVANTK_EVAL_DATA` is the
  site profile's value, default `$HVANTK_WORK/hvantk-eval`): `reference.mt` (the label
  is the flat column `super_pop`) and `query.mt` (no labels). Both are GRCh38 and share
  their variants. Do not download anything.
- Write the outputs under `$HVANTK_EVAL_OUT/`: at least `predictions.ht` and
  `predictions.tsv`, plus the pipeline's statistics file.
- The reference and query share about 20,000 variants, fewer than the pipeline's default
  minimum; choose the pipeline's options accordingly.
- The harness sets `HVANTK_EVAL_OUT` and `HVANTK_RUN_ID` in your environment before you
  start. Name every Slurm job you submit with the prefix `$HVANTK_RUN_ID`.
- Budget per job: one node, at most 16 cores, 64 GB of memory, 60 minutes.
- When done, write `$HVANTK_EVAL_OUT/run-notes.md` with one line `commit: <git commit
  of the checkout>`, one line `container: <sha256 of the image file>`, one line
  `job: <id>` per job you submitted, and the shared-variant count, the cross-validated
  accuracy and the assigned fraction the pipeline reported.

You are done when the predictions exist for every query sample, the notes are written,
and every job you submitted has left the queue.
