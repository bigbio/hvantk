# Sizing a build, by dataset class

Read `hpc-hail-job`'s `references/sizing.md` first for the three rules (driver memory
below the cgroup limit, cores equal the allocation, scratch is node-local) — they apply
unchanged to every task in a build array. In particular, every array task runs through
`hail-job.sbatch` exactly as a single `hpc-hail-job` submission would: Spark's driver is
capped at 80 percent of `--mem-gb`, and `local[<cpus>]` is pinned to `--cpus`, both per
task, not per array. This page adds what is specific to `hvantk reprocess`: what drives
memory for each dataset class, and when a dataset does not belong in the array at all.

## What drives memory in a build

- **Import of the raw file.** `hl.import_table` / `hl.import_vcf` reads the whole file
  into Spark partitions before anything else happens; a wide file (many columns, HGNC's
  54) costs more than a narrow one of the same row count.
- **Joins.** A builder that joins against a gene catalog (ClinGen and GenCC both key on
  `hgnc_id`) holds both sides in memory during the shuffle.
- **Hail's key sort.** `key_by` on a multi-field key (ClinGen's `(hgnc_id, mondo_id)`,
  GenCC's `(hgnc_id, mondo_id, submitter)`) sorts the whole table; a wider key costs
  more than HGNC's single-field one for a similar row count.

None of the three eval datasets (HGNC, ClinGen, GenCC) need more than a gene-level
request; they are the small end of this table on purpose.

## Starting points per class

| class | example datasets | cores | memory | time | notes |
| --- | --- | --- | --- | --- | --- |
| gene-level table | `hgnc:lookup`, `clingen:gene-disease`, `gencc:submissions`, `ensembl-gene:structure` | 4 | 16 GB | 30 min | tens of thousands of rows; one array request covers several of these at once |
| variant-level table | `clinvar:variants` | 8 | 64 GB | 2 h | millions of rows; import dominates |
| variant-level, very large | `dbnsfp:variants` | whole fat node | whole fat node | a day | does not belong in an array; its own `hpc-hail-job` submission |
| expression matrix (`.h5ad`) | `expression-atlas:dataset`, `cptac:expression` | 8 | 32 GB | 1 h | the matrix is held in memory whole; size for the largest study, not the average one |
| QTL summary stats, per tissue | `gtex-eqtl:eqtls` | 8 | 64 GB | 2 h | one array task per tissue is a reasonable use of this skill |

## One request for the whole array

`submit-build-array.sh` sizes every task in the array identically — there is one
`--cpus`/`--mem-gb`/`--time`, not one per line of the task list. Size for the list's
largest dataset. A dataset that would need a request order-of-magnitude larger than the
rest (dbNSFP next to HGNC, say) does not belong in the same array: give it its own
`hpc-hail-job` submission and keep the array for datasets that share a size class.

## `--max-concurrent` against the node share

`--max-concurrent` (the `%K` in `--array=1-N%K`) caps how many tasks run at once, which
caps how many nodes (or how much of a shared node pool) the array can occupy
simultaneously. On a partition you share with other users, set it to a fraction of what
you could plausibly claim without queueing everyone else out — not to `N`, even when the
site's node count would allow it.

## Node-local scratch, per task

Each array task gets its own node-local scratch directory from `hail-job.sbatch`
(`$SLURM_TMPDIR`-derived, named with the task's own `$SLURM_JOB_ID`), so concurrent
tasks never collide on scratch even when `--max-concurrent` lands two of them on the
same node. Nothing extra to configure here; it falls out of running each task through
the same per-job template `hpc-hail-job` uses.

## Outputs and sidecars

Each dataset's output (a `.ht` directory on a parallel filesystem) is many part files,
like any Hail table. Keep the dated raw file that produced it and the output's
`.provenance.json` sidecar next to the output, or at least discoverable from it: the
sidecar is the only record of which raw file, plugin version and commit built it, and it
does not live inside the output directory (moving one without the other loses that
link).

## Re-running only failed indices

A failed task in an array does not fail the array job as a whole, and does not stop the
other tasks. Read `sacct -j <A> -X -n -o JobID,State` to find which indices did not
reach `COMPLETED`, fix whatever caused it (usually a raw-dir problem; see the skill's
Traps section), and resubmit only those with `submit-build-array.sh --list <the same
list> --tasks i,j ...` — a fresh array job id, not a retry of the old one.
