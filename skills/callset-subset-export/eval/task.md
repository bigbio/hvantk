# Task: deliver a two-sample subset of a joint callset as a verified VCF

A three-sample joint callset exists as a lossless VDS. Deliver the two samples listed in
`samples.txt` as a biallelic VCF with checksums and a manifest, by subsetting the VDS,
on this cluster, as Slurm jobs inside the container image.

- The site profile is at `$HVANTK_SITE_PROFILE`.
- The callset is `$HVANTK_EVAL_DATA/callset-subset-export/cohort/cohort.vds` and the
  list is `$HVANTK_EVAL_DATA/callset-subset-export/samples.txt`, where `HVANTK_EVAL_DATA`
  is the site profile's value (default `$HVANTK_WORK/hvantk-eval`). The list uses the
  callset's own sample spelling.
- Write the delivery to `$HVANTK_EVAL_OUT/subset/`: at least `subset.vcf.bgz`,
  `SHA256SUMS` covering the delivered files, and `README.md` describing what was done.
- The harness sets `HVANTK_EVAL_OUT` and `HVANTK_RUN_ID` in your environment before you
  start. Name every Slurm job you submit with the prefix `$HVANTK_RUN_ID`.
- Budget per job: one node, at most 8 cores, 32 GB of memory, 30 minutes.
- When done, write `$HVANTK_EVAL_OUT/run-notes.md` with one line `commit: <git commit
  of the checkout>`, one line `container: <sha256 of the image file>`, one line
  `job: <id>` per job you submitted, and the sample set and record count you verified.

You are done when the delivery exists with a validating checksum file, the notes are
written, and every job you submitted has left the queue.
