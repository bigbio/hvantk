# Feature axes

A feature axis is any gene-keyed table (TSV, CSV, `.bgz`/`.gz`, or `.parquet`;
`hvantk/algorithms/rerank/catalog/builders.py` `table_axis`) with a column literally
named `gene` plus one or more feature columns. `table_axis`'s loader does nothing
clever: `pd.read_csv`/`pd.read_parquet`, then every non-`gene` column is coerced with
`pd.to_numeric(..., errors="coerce")` inside `engine.rerank()` -- a non-numeric cell
(or an entire text column) becomes `NaN`, silently, not an error. Check a new axis's
dtypes before trusting it.

## The baseline axis

The FIRST axis listed under `features:` in the config is the ablation baseline
(`engine.rerank()`: `baseline = next(iter(groups))`, and Python dict iteration order
follows insertion order, which follows the YAML list order). It stays in the full
model like every other axis -- "baseline" only means every other axis's ablation
delta is computed against it. A one-axis config is valid (the eval in this skill's
`eval/` runs exactly one, `structure`) and produces a trivial, single-row ablation
table, because there is nothing left to ablate against the baseline but itself.

## Building an axis from an hvantk dataset

- **Gene structure.** The gene index itself (`hvantk annotate spine`; see the
  `gene-index` skill) carries `cds_length`, `n_coding_exons`, `n_transcripts` per
  gene -- `prepare-rerank-inputs.py features --ht <gene-index.ht> --fields
  cds_length,n_coding_exons,n_transcripts ...` is the whole recipe. This is what this
  skill's own eval ranks on.
- **Constraint.** `gnomad-metrics:metrics` (gnomAD's per-gene constraint table,
  keyed `gene_id`) or `gevir:metrics`. Neither carries `hgnc_id` directly the way
  the gene index does, so point `prepare-rerank-inputs.py features`'s `--gene-index`
  at the gene index separately from `--ht` (which stays the constraint table) so the
  `gene_id -> symbol` lookup still has something to join through.
- **Expression.** An expression `.h5ad` summarised per tissue or cell type with
  `hvantk expression summarize` (wide output, one row per gene) -- pick the
  tissues/cell types that actually fit the phenotype under study; an unrelated
  tissue is not wrong, exactly, but it adds noise the ablation table will then have
  to explain.

## Naming columns uniquely across axes

`FeatureAssembler.assemble` raises if two axes declare the same non-`gene` column
name (`rerank feature axis '{ax}' declares column(s) [...] already provided by axis
'{other}'`) -- the outer join would otherwise silently suffix both to `_x`/`_y`,
feeding the model the same logical quantity twice under two names and breaking the
per-axis grouping the ablation and the audit both rely on (they match columns by
exact name). Give every feature table's columns a name unique across every axis you
plan to combine, even axes you have not written yet.

## Circularity

An axis's own values must not have been derived, even partially, from the SAME
curated resource the labels come from. The canonical trap: a "known-disease-gene"
feature axis built from ClinGen or GenCC, scored against ClinGen- or GenCC-derived
labels -- the axis is then measuring "was this gene already in the database", which
trivially predicts "is this gene in the database" and inflates AUC with nothing
learned. `references/audit-and-leakage.md` covers the one formal tool for this
(`feature_provenance`/`label_provenance`, Python API only) and the informal check
(read what each axis is actually built from before adding it) that this skill's own
eval relies on, since its one axis (`structure`, built from Ensembl/HGNC structural
fields) shares no provenance with its labels (ClinGen gene-disease validity).
