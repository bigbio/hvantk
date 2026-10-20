# The image provenance record

`templates/build-image.sbatch` writes `<image>.provenance.json` next to every image it
builds, before the final pass-or-fail gate, so a failed install check is recorded too.

```json
{
  "sif": "<absolute path of the image>",
  "built_utc": "2026-10-15T07:42:10Z",
  "slurm_job": "<job id, empty when built outside a job>",
  "host": "<node that built it>",
  "repo": "<checkout the files were copied from>",
  "git_commit": "<full sha>",
  "git_branch": "<branch>",
  "worktree": "clean | dirty",
  "poetry_lock_md5": "<md5 of poetry.lock at build time>",
  "container_tool": "singularity-ce version 4.x | apptainer version 1.x",
  "build_method": "fakeroot | unprivileged + proot",
  "hail_version": "0.2.137-<hash>",
  "java_version": "11.0.x",
  "check_install_rc": 0
}
```

How to read it:

| field | what it answers | what to do when it is wrong |
| --- | --- | --- |
| `git_commit`, `git_branch`, `worktree` | which code's lock file and package went in; `dirty` means uncommitted edits went in too | rebuild from a clean tree at the intended commit |
| `poetry_lock_md5` | whether the image matches the checkout's dependencies: compare with `md5sum "$HVANTK_REPO/poetry.lock"` | the clone was stale at build time, or the lock changed since: sync, rebuild |
| `hail_version`, `java_version` | the two pins Hail is most sensitive to; `hail_version` starts with the locked version and carries a build hash after it | a different Hail or a Java other than 8 or 11 means the def file or the lock drifted |
| `check_install_rc` | whether Hail started and ran its smoke test inside the image on a compute node | anything but 0 is not a usable image; read the job log |
| `build_method`, `container_tool` | how it was built, for a site that changes its runtime | informational |

The record is a sidecar, like the `.provenance.json` hvantk writes next to every dataset:
a job's run notes name the image digest, the image's record names the commit, and the
dataset's record names the source fingerprint. Inside an image built from
`templates/build-image.sbatch`, a dataset's own `builder_commit` field records that same
commit too: the def file's `%arguments` section stamps it into `HVANTK_BUILDER_COMMIT`
at build time, and `reprocess` reads it there once its own git check (correctly) finds no
checkout to ask (see `dataset-build-at-scale`'s verification reference). An image built
without `--build-arg HVANTK_COMMIT=...` -- a manual `singularity build` with no flag --
falls back to recording `hvantk <version>` instead. Together they answer "what produced
this file" without a conversation.

Promotion keeps records with images: `promote-image.sh` moves the old default's record
with the old default and points the default's record at the new one, so
`source "$HVANTK_SITE_PROFILE" && cat "$HVANTK_SIF.provenance.json"` always describes the
image the profile's default resolves to.
