# Task: rank a cohort's genes by gene structure and evaluate against ClinGen labels

Rank the genes a synthetic cohort tested, by one feature axis built from the gene
index's structural fields, and evaluate the ranking against ClinGen gene-disease
validity labels, on this cluster inside the hvantk container image.

- The site profile is at `$HVANTK_SITE_PROFILE`.
- Inputs are already staged under `$HVANTK_EVAL_DATA/rerank-evaluate/` (where
  `HVANTK_EVAL_DATA` is the site profile's value, default `$HVANTK_WORK/hvantk-eval`).
  Do not download anything; compute nodes have no internet access.
  - `cohort/gene_burden.tsv` -- the cohort's gene-level table, one row per gene
    tested, keyed by an unversioned Ensembl `gene_id` in its `gene` column.
  - `cohort/cohort.yaml` -- the matching cohort manifest (`key: gene_id`,
    `key_column: gene`, prior column `minp`, and an `architecture` axis of
    `n_case_var`, `conc`, `driver_af`).
  - `gene-index.ht` -- the gene index (keyed `gene_id`; carries `hgnc_id`,
    `cds_length`, `n_coding_exons`, `n_transcripts`, among other fields).
  - `raw/clingen/` -- one dated ClinGen gene-disease CSV
    (`Clingen-Gene-Disease-Summary-<date>.csv`).
  - `raw/hgnc/hgnc_complete_set.txt` -- the raw HGNC complete-set TSV.
- **The ranked universe is the genes `cohort/gene_burden.tsv` tested** -- every gene
  in that table's `gene` column, nothing more and nothing fewer.
- **Feature axis**, named `structure`: the gene index's `cds_length`,
  `n_coding_exons` and `n_transcripts`, for exactly the genes in the ranked universe
  above.
- **Labels**: ClinGen gene-disease validity at Moderate or stronger (any disease --
  do not restrict to one GCEP or one disease term), built from the staged ClinGen
  download, restricted to the ranked universe above.
- **Paralogue blocks**: `--blocks` from the staged HGNC complete-set TSV.
- **Seed**: 42.
- The cohort's architecture axis (`n_case_var`, `conc`, `driver_af`) must still be
  present and intact wherever your config's `cohort:` manifest ends up pointing, so
  the case/control architecture audit runs.
- Write the config to `$HVANTK_EVAL_OUT/config.yaml` and the ranked table to
  `$HVANTK_EVAL_OUT/ranked.tsv`. Evaluate the ranked table and keep its receipt
  somewhere under `$HVANTK_EVAL_OUT/`.
- Name every Slurm job you submit with the prefix `$HVANTK_RUN_ID`. The harness sets
  `HVANTK_RUN_ID` and `HVANTK_EVAL_OUT` (a directory under `HVANTK_EVAL_DATA`, so jobs
  can see it) in your environment before you start.
- Budget per job: one node, at most 4 cores, 16 GB of memory, 45 minutes.
- When everything has completed, write `$HVANTK_EVAL_OUT/run-notes.md` containing one
  line `commit: <git commit of the checkout>`, one line `container: <sha256 of the
  image file>`, one line `job: <id>` for every job you submitted, and one line
  `auc: <the ROC AUC of score against y over the ranked table>`.

You are done when `config.yaml` and `ranked.tsv` exist, the evaluation receipt
exists, the run notes are written, and every job you submitted has left the queue.
