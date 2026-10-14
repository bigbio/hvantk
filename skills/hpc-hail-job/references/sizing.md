# Sizing a single-node Hail job

Every number here is a starting point for one node of a general CPU partition. Measure
once, then adjust.

## The three rules

1. **Driver memory below the cgroup limit.** In local mode the Spark driver is the only
   JVM, and Hail also uses off-heap memory (native libraries, compression buffers) plus
   Python workers outside the JVM. If the driver heap equals the Slurm request, the
   kernel kills the job silently. The template sets the driver to 80 percent of
   `--mem`. Raise the request, not the percentage.
2. **Cores equal the allocation.** `local[<cpus-per-task>]`, never `local[*]` on a
   shared node. On an `--exclusive` node both are the same.
3. **Scratch is node-local.** Spark spill and Hail temporary files go to the node's
   local disk (`$SLURM_TMPDIR` on most sites), never to home or to the parallel
   filesystem.

## Starting points

| workload | cores | memory | time | notes |
| --- | --- | --- | --- | --- |
| install check, describe a table | 2 | 8 GB | 15 min | |
| gene-level table (tens of thousands of rows) | 4 | 16 GB | 30 min | |
| gene-level with heavy joins (gene index, cohort attach) | 8 | 32 GB | 1 h | |
| variant-level table from one provider | 16 | 120 GB | 4 h | input size drives memory |
| QTL cascade over one tissue | 8 | 64 GB | 2 h | |
| joint genotyping, hundreds of samples | whole node, exclusive | whole node | 12 h | one job per contig scales better than one job for all |

## Reading a failure

- Finished early with `OUT_OF_MEMORY` or a log that stops without a traceback: raise
  memory by half, keep cores.
- `TIMEOUT`: the elapsed time is in `sacct -j <id> -X -o Elapsed,Timelimit`; ask for
  what it needed plus half. Hail does not write partial tables; a timeout wastes the
  whole run.
- Slow but finishing: check the scratch was node-local (`SPARK_LOCAL_DIRS` in the log
  banner) before adding cores.

## What the template exports

| variable | value | why |
| --- | --- | --- |
| `SPARK_LOCAL_DIRS`, `TMPDIR`, `HAIL_TMPDIR` | the node-local scratch directory | rule 3 |
| `PYSPARK_SUBMIT_ARGS` | `--master local[<cpus>] --driver-memory <80% of mem> --conf spark.local.dir=<scratch> pyspark-shell` | rules 1 and 2; hvantk starts Hail with defaults and PySpark reads this variable |
| `OMP_NUM_THREADS` and the BLAS equivalents | `<cpus>` | libraries otherwise size their pools to the node |
| `JAVA_TOOL_OPTIONS` | `-Duser.language=en -Duser.country=US` | decimal separators in exported text |
| `PYTHONPATH` | the checkout (`HVANTK_REPO`) | the checkout's code runs inside the image's environment |
