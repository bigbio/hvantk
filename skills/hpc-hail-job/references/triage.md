# Triage for a Slurm job

## Where things are

- Job log: `${HVANTK_LOG_DIR:-$HVANTK_WORK/hvantk-logs}/<name>_<id>.out` (stdout and stderr
  together). Hail's own `hail-<timestamp>.log` is written to the same directory --
  except for an array task, which gets its own subdirectory,
  `<name>_<array-job-id>_<task-id>/`, so the task's Hail log does not collide with
  its siblings'.
- Accounting: `sacct -j <id> -X -o JobID,JobName%30,State,ExitCode,Elapsed,Timelimit,ReqMem,AllocCPUS,NodeList`.
- While it runs: `squeue -j <id> -o "%i %j %T %M %L %R"` (state, elapsed, time left, reason).
- Details of the request: `scontrol show job <id>`.

## States and what to do

| state | meaning | action |
| --- | --- | --- |
| `PENDING (Priority)` | waiting in the queue | wait; `squeue --start -j <id>` estimates |
| `PENDING (Resources)` | no node fits the request | check the request against `sinfo -p <partition> -o "%c %m"` |
| `PENDING (QOSMaxJobsPerUserLimit)` or similar | your own job cap | wait for your other jobs |
| `PENDING (ReqNodeNotAvail)` | nodes down or reserved | ask the site; do not resubmit in a loop |
| `RUNNING` | running | read the log banner: scratch, cores, driver memory |
| `COMPLETED` | exit 0 | verify the output; completion is not correctness |
| `FAILED` | non-zero exit | the Python traceback is at the end of the log |
| `OUT_OF_MEMORY` | the cgroup killed it | raise memory by half |
| `TIMEOUT` | the time limit hit | ask for elapsed plus half |
| `CANCELLED` | someone cancelled it | `sacct` shows by whom; a site cancels jobs that write to home quotas |
| `NODE_FAIL` | the node died | resubmit unchanged |

## Failure signatures in the log

| text in the log | cause | fix |
| --- | --- | --- |
| `Failed to create any local dir` then `Connection refused` | no writable scratch for Spark | run through the template |
| `java.lang.OutOfMemoryError: Java heap space` | driver heap too small for the input | raise `--mem-gb` |
| `Executor heartbeat timed out` | Spark oversubscribed the node | pin cores (the template does) |
| `UnsupportedClassVersionError` or a Java version warning from Hail | host Java leaked in | no modules around the container |
| `No such file or directory` on an input that exists on the login node | path not visible on compute nodes, or a bind is missing | check the mount from a job; set `HVANTK_BIND` in the site profile (a shell export is overwritten when the job sources the profile) |
| `Read-only file system` | the output mount is read-only on compute nodes | write under `$HVANTK_WORK` |
| `HTTP` or `URLError` during `reprocess` | compute node has no internet | stage on the login node; `--skip-download --no-check-drift` |
| `is already initialized` from Hail | a second Hail start in one process | run one hvantk command per job |

## Cancelling

`scancel <id>` for one job; `scancel --name=<prefix>...` does not take a prefix, so cancel
by id. Cancel your own stuck waiters too: a shell loop left on the login node is the
quieter cousin of a runaway job.
