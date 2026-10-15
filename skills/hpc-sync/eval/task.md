# Task: bring a cluster remote in step with GitHub without losing its own commits

A repository lives on two remotes: a GitHub remote, which is the source of truth, and a
bare repository on a cluster, which is behind and also carries one commit of its own.
Make the cluster remote's `main` and `dev` point at exactly the commits they point at
on GitHub, keep every commit that exists only on the cluster reachable, and do not
change anything on the GitHub remote.

- The scenario is at `$HVANTK_EVAL_OUT/scenario/`: `github.git` (bare, the source of
  truth), `cluster.git` (bare, behind, with a branch of its own), and `laptop/`, a
  working clone of `github.git` whose second remote `cluster` points at `cluster.git`.
  Its local branches may not be current. Work from `laptop/`.
- The harness sets `HVANTK_EVAL_OUT` and `HVANTK_RUN_ID` in your environment before you
  start, and `HVANTK_SITE_PROFILE` names the site profile.
- No Slurm job is needed for this task. Still record the checkout and image in your
  notes, as every workflow does.
- When done, write `$HVANTK_EVAL_OUT/run-notes.md` containing one line `commit: <git
  commit of the hvantk checkout>`, one line `container: <sha256 of the image file>`,
  and the ref report you used to verify the result.

You are done when both remotes agree on `main` and `dev`, the cluster's own branch is
untouched, and the notes are written.
