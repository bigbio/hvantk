# Task: burden-test a cohort and attach it to the gene index

A synthetic case/control cohort was prepared from one 1000 Genomes chromosome: a clean,
split-biallelic, rare-variant MatrixTable, already annotated with a gene and a variant
route. Run hvantk's cohort burden test on it, write and validate the cohort manifest,
then attach the result onto the gene index -- each step on this cluster, inside the
container image.

There is no network access from this environment: nothing can be downloaded or probed;
the inputs are already staged.

- The site profile is at `$HVANTK_SITE_PROFILE`.
- Inputs under `$HVANTK_EVAL_DATA/cohort-attach-burden/` (`HVANTK_EVAL_DATA` is the
  site profile's value, default `$HVANTK_WORK/hvantk-eval`):
  - `cohort.mt`: the case/control MatrixTable. Row fields `gene_id` (Ensembl, string,
    unversioned) and `route` (string: `snv` or `indel`); column field `is_case`
    (boolean); entry field `GT`. Already split-biallelic and rare; no further
    cleaning needed.
  - `gene-index.ht`: the gene index to attach onto, keyed by `gene_id`.
- Produce, under `$HVANTK_EVAL_OUT/`:
  1. `gene_burden.tsv`: the gene-level burden test, counting a sample as a carrier
     when it has a **heterozygous** qualifying genotype, with a
     **Benjamini-Hochberg** multiple-testing correction applied to the result.
  2. `cohort.yaml`: a cohort manifest for that table, with:
     - `name: eval-cohort`;
     - the gene key in the **Ensembl gene id** space (`gene_id`);
     - the prior statistic the ranking is built on is the burden test's own min-p
       column, with the direction stated as "lower is better";
     - an axis named `architecture` holding `n_case_var`, `conc`, `driver_af`;
     - an **absolute** path for `table:` (a relative one is resolved against the
       working directory of whatever process reads it next, not against the
       manifest file, so it would break the moment validation or attach runs from
       a different directory than the one it was written in).
  3. Validated (`hvantk cohort validate`).
  4. `attached.ht`: that cohort attached onto `gene-index.ht`, with its report at the
     default path (`attached.ht.report.json`).
- Name every Slurm job you submit with the prefix `$HVANTK_RUN_ID`. The harness sets
  `HVANTK_RUN_ID` and `HVANTK_EVAL_OUT` (a directory under `HVANTK_EVAL_DATA`, so jobs
  can see it) in your environment before you start.
- Budget per job: one node, at most 8 cores, 32 GB of memory, 30 minutes.
- When every job has completed, write `$HVANTK_EVAL_OUT/run-notes.md` containing one
  line `commit: <git commit of the checkout>`, one line `container: <sha256 of the
  image file>`, one line `job: <id>` for every job you submitted, and one line
  `attached: <n>/<m>` copied verbatim from the attach job's own printed
  `attach eval-cohort: <n>/<m> genes tested -> ...` line (read it from that job's log).

You are done when `gene_burden.tsv`, `cohort.yaml` and `attached.ht` all exist, the run
notes are written, and every job you submitted has left the queue.
