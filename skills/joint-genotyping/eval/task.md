# Task: joint-genotype three gVCFs into a cohort callset

Run hvantk's HGC pipeline on this cluster, as a Slurm job inside the container image,
on the three staged gVCFs, and produce the cohort callset with its QC tables and VCF.

- The site profile is at `$HVANTK_SITE_PROFILE`.
- The gVCFs (bgzipped, indexed, GRCh38, `chr`-named contigs) are staged under
  `$HVANTK_EVAL_DATA/joint-genotyping/gvcfs/`, where `HVANTK_EVAL_DATA` is the site
  profile's value (default `$HVANTK_WORK/hvantk-eval`). Do not download anything.
- Write the outputs under `$HVANTK_EVAL_OUT/cohort/` with the output prefix `cohort`,
  so that the directory holds `cohort.vds`, `cohort.mt`, `qc/` and `cohort.vcf.bgz`.
- The harness sets `HVANTK_EVAL_OUT` and `HVANTK_RUN_ID` in your environment before you
  start. Name every Slurm job you submit with the prefix `$HVANTK_RUN_ID`.
- Budget per job: one node, at most 8 cores, 32 GB of memory, 30 minutes.
- When done, write `$HVANTK_EVAL_OUT/run-notes.md` containing one line `commit: <git
  commit of the checkout>`, one line `container: <sha256 of the image file>`, one line
  `job: <id>` for every job you submitted, and the three counts you verified: the sample
  set, the MatrixTable row count and the VCF record count.

You are done when the callset exists with its QC tables and VCF, the notes are written,
and every job you submitted has left the queue.
