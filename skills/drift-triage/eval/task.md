# Task: triage a drifted ClinGen signal and repair the plugin

`hvantk drift clingen:gene-disease` has reported `drifted`. Triage the signal and get
the dataset building again.

- The site profile is at `$HVANTK_SITE_PROFILE`. It names the checkout, the image, the
  partition and the work root.
- The drift report is already staged at
  `$HVANTK_EVAL_DATA/drift-triage/drift-report.json`, where `HVANTK_EVAL_DATA` is the
  site profile's value (default `$HVANTK_WORK/hvantk-eval`). The upstream file it was
  compared against is staged at `$HVANTK_EVAL_DATA/drift-triage/raw/clingen/` (one CSV).
  There is no network access from here, so that staged file is what upstream now serves
  -- do not try to reach it over the network.
- The checkout to repair is `$HVANTK_REPO`, in place. Make whatever code change the
  signal calls for.
- Build `clingen:gene-disease` from the staged raw directory into
  `$HVANTK_EVAL_OUT/clingen-gene-disease.ht`, as a Slurm job inside the hvantk container
  image. Do not build on the login node.
- The plugin's committed schema snapshot (`hvantk/skills/clingen/tests/snapshots/schema.json`)
  is the contract: the built table must carry its key and every field it declares. Do
  not edit the snapshot files to make the build match.
- Record the outcome the way hvantk's drift workflow requires: the dataset's
  fingerprint file must reflect the header the staged file actually has, the rebuild
  ledger must carry a row for this change, and that row must say the dataset was
  rebuilt.
- Name every Slurm job you submit with the prefix `$HVANTK_RUN_ID`. The harness sets
  `HVANTK_RUN_ID` and `HVANTK_EVAL_OUT` (a directory under `HVANTK_EVAL_DATA`, so jobs
  can see it) in your environment before you start.
- Budget per job: one node, at most 4 cores, 16 GB of memory, 30 minutes.
- When you are done, write `$HVANTK_EVAL_OUT/run-notes.md` containing one line
  `commit: <git commit of the checkout>`, one line `container: <sha256 of the image
  file>`, and one line `job: <id>` for every job you submitted.

You are done when the dataset exists, the checkout carries the repair and the records,
the run notes are written, and every job you submitted has left the queue.

Leave the repair and the records in the working tree; do not create commits.
