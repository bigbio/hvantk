# Task: build the container image from the checkout and prove it works

Build the hvantk container image from the checkout's `containers/hvantk.def` and lock
file on this cluster, without root, as a Slurm job; then prove it works by running
`hvantk utils check-install` inside the new image in a job.

- The site profile is at `$HVANTK_SITE_PROFILE`. It names the checkout, the current
  default image, the partition and the work root.
- Write the new image to `$HVANTK_EVAL_OUT/hvantk-eval.sif`, with its provenance record next
  to it. Do not touch the default image.
- The harness sets `HVANTK_EVAL_OUT` and `HVANTK_RUN_ID` in your environment before you
  start. Name every Slurm job you submit with the prefix `$HVANTK_RUN_ID`.
- Compute nodes on this site can reach the container registry and the package index.
- Budget per job: one node, at most 8 cores, 32 GB of memory, 105 minutes.
- When done, write `$HVANTK_EVAL_OUT/run-notes.md` containing one line `commit: <git
  commit of the checkout>`, one line `container: <sha256 of the NEW image file>`, and
  one line `job: <id>` for every job you submitted.

You are done when the image exists, its record shows a passing install check, the
notes are written, and every job you submitted has left the queue.
