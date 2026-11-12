# Task: build a labelled, leakage-tiered ClinVar training set for chromosome 22

Build an ML training set from a feature spec with `hvantk training-set build`, on this
cluster inside the hvantk container image: ClinVar missense SNVs on chromosome 22,
labelled pathogenic/benign, with dbNSFP scores, protein-site flags, gene-level features
and a gene-list flag, every feature with a leakage tier.

- The site profile is at `$HVANTK_SITE_PROFILE`.
- Inputs are already staged under `$HVANTK_EVAL_DATA/ml-features/` (where
  `HVANTK_EVAL_DATA` is the site profile's value, default `$HVANTK_WORK/hvantk-eval`).
  Do not download anything; compute nodes have no internet access.
  - `spec.yaml` -- the feature spec. Use it as it is.
  - `genes.txt` -- the gene list for the spec's `is_eval_gene` column.
  - `raw/clinvar/`, `raw/dbnsfp/`, `raw/insider_variants/`, `raw/insider_interfaces/`,
    `raw/uniprot_ptm/`, `raw/gnomad_metrics/`, `raw/gevir/` -- one raw file each, for
    the datasets `clinvar:variants`, `dbnsfp:variants`, `insider:variants`,
    `insider:interfaces`, `uniprot-ptm:sites`, `gnomad-metrics:metrics` and
    `gevir:metrics`. The ClinVar, dbNSFP and INSIDER-variants files hold chromosome 22
    only.
  - `upstream/gene-index.ht` -- the gene index (keyed `gene_id`).
  - `upstream/hgnc-lookup.ht` -- the `hgnc:lookup` dataset the gene index was built
    with.
- Build each of the seven datasets yourself from its raw file, and the gene matrix the
  spec's `layer1` section declares, on the gene index above. Write every output under
  `$HVANTK_EVAL_OUT`.
- Write the training set with the output prefix `$HVANTK_EVAL_OUT/clinvar_chr22`, so
  that `clinvar_chr22.ht`, `clinvar_chr22.tsv.gz`, `clinvar_chr22.dictionary.tsv` and
  `clinvar_chr22.manifest.json` exist in `$HVANTK_EVAL_OUT`.
- Name every Slurm job you submit with the prefix `$HVANTK_RUN_ID`. The harness sets
  `HVANTK_RUN_ID` and `HVANTK_EVAL_OUT` (a directory under `HVANTK_EVAL_DATA`, so jobs
  can see it) in your environment before you start.
- Budget per job: one node, at most 4 cores, 16 GB of memory, 45 minutes.
- When everything has completed, write `$HVANTK_EVAL_OUT/run-notes.md` containing one
  line `commit: <git commit of the checkout>`, one line `container: <sha256 of the
  image file>`, one line `job: <id>` for every job you submitted, one line
  `rows: <number of rows in the training set>` and one line
  `positives: <number of rows labelled 1>`.

You are done when the four training-set files exist, the run notes are written, and
every job you submitted has left the queue.
