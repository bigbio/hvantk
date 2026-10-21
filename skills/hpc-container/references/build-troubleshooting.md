# Build troubleshooting

## Minimum container tool

The def file's `%arguments` section and the template's `--build-arg` need SingularityCE
4.0 or later, or Apptainer 1.2 or later. Check with `"$HVANTK_CONTAINER_BIN" --version`
before the first build on a site; an older tool fails at once (see the failure table).

## Choosing the unprivileged method

| `grep -cE "^($USER|$(id -u)):" /etc/subuid` | method | flag |
| --- | --- | --- |
| 1 or more | fakeroot (user namespaces mapped for you) | `--fakeroot` |
| 0, or the file does not exist | a static `proot` on `PATH` with SingularityCE (4.0 or later, as above); Apptainer releases may build without either, so check `apptainer build --help` and your site's notes | none |

Getting `proot` on a cluster without it: download the static binary from the proot
project's release page into a directory on your `PATH` (for example `~/bin`), make it
executable, and check `proot --version`. Record the version in the run notes. If the site
forbids downloading binaries, ask the administrators for `/etc/subuid` entries instead.

## Failure signatures

| text in the job log | cause | fix |
| --- | --- | --- |
| `unknown flag: --build-arg` | SingularityCE older than 4.0 or Apptainer older than 1.2 | use a newer container tool (a site module, or build elsewhere and copy the image) |
| `no valid mapping entry found` | fakeroot requested without a subuid entry | use proot (table above); the template chooses automatically |
| `--remote, --fakeroot, or the proot command are required` | neither method available | install proot on `PATH` |
| `cp: cannot stat 'pyproject.toml'` or similar during `%files` | the build did not start from the repository root | the template `cd`s there; by hand, `cd "$HVANTK_REPO"` first |
| `Package 'openjdk-11-jdk-headless' has no installation candidate` | a base image without Java 11 | keep the def file's base; do not edit it to a newer Debian |
| `apt-get update` succeeds and every install 404s | a base whose distribution left support | the def file pins a supported base; an old def must be updated in the repository, not patched locally |
| `FATAL: python3 is not 3.10` | the base image's Python moved | the lock file was resolved for 3.10; a base change needs a lock review |
| `Failed to create any local dir` then `Connection refused` in the check step | no node-local scratch for Spark | the template sets `SPARK_LOCAL_DIRS`; by hand use `containers/hvantk_run.sh` on a compute node |
| the pull hangs, then times out | no egress from the node | the template probes the registry first; build where egress exists and copy the file |
| `No space left on device` under the home directory | cache or tmp on the home quota | the template uses `$HVANTK_WORK/containers/{cache,tmp}`; by hand export `SINGULARITY_CACHEDIR` and `SINGULARITY_TMPDIR` (or the `APPTAINER_` pair) |
| `module 'attr.setters' has no attribute 'pipe'` during pip | pip's vendored `rich` racing a concurrent attrs install | the def file pins pip's upgrade before poetry; a retry usually passes |

## How long it takes

A build is dominated by the base image pull and the wheel downloads: around ten to
twenty minutes on a shared partition with warm caches, longer the first time. The
install check inside the job takes a minute or two (JVM start, a tiny Hail model). The
template requests 1h45 so a slow registry does not kill a build that would have
finished.

## Checking an image by hand

On a compute node (never the login node), through the wrapper that sets up scratch:

```bash
source "$HVANTK_SITE_PROFILE" && HVANTK_SIF="<image>" bash "$HVANTK_REPO/containers/hvantk_run.sh" utils check-install
```

Expect the Hail version and `All checks passed.` Then `singularity inspect <image>` for
the labels the def file sets (Hail and pyspark versions, base image).
