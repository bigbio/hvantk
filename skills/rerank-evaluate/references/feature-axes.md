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
delta is computed against it. A one-axis config is valid and produces a trivial,
single-row ablation table, because there is nothing left to ablate against the
baseline but itself -- that is a property of running one axis, not a sign the run
is broken.

## Building an axis from an hvantk dataset

- **Gene structure.** The gene index itself (`hvantk annotate spine`; see the
  `gene-index` skill) carries `cds_length`, `n_coding_exons`, `n_transcripts` per
  gene -- `prepare-rerank-inputs.py features --ht <gene-index.ht> --fields
  cds_length,n_coding_exons,n_transcripts ...` is the whole recipe.
- **Constraint.** `gnomad-metrics:metrics`'s DEFAULT build (gnomAD v2.1.1's
  `by_gene` table, `hvantk download gnomad-metrics`'s own default) is keyed
  `gene_id`, one row per gene, usable with `prepare-rerank-inputs.py features` the
  same way as gene structure. It does not carry `hgnc_id` directly the way the gene
  index does, so point `--gene-index` at the gene index separately from `--ht`
  (which stays the constraint table) so the `gene_id -> symbol` lookup still has
  something to join through. gnomAD's v4.0 release (`--version v4.0`,
  `constraint_metrics`) is PER-TRANSCRIPT, not per-gene (`hvantk/skills/
  gnomad_metrics/builder.py`'s own comment: it has no `gene_id` column at all, and
  the builder's `key` parameter has to be set to `transcript` to read it) -- this
  skill's templates have no recipe for reducing that to one row per gene (e.g. the
  MANE/canonical transcript); build v2.1.1 instead, or write that reduction
  yourself before using v4.0 as an axis. `gevir:metrics` is the alternative,
  per-gene constraint source.
- **Expression.** `hvantk expression summarize` writes a per-group-by-per-gene
  AnnData (`.h5ad`, shape `(n_groups, n_genes)`) with layers `mean`, `sum`,
  `count_nonzero`, `fraction_expressed` -- genes are columns (`.var`), groups are
  rows (`.obs`), so one tissue/cell-type axis is one ROW of one layer, which still
  needs transposing into a one-row-per-gene table before it is a feature axis in
  this sense. `prepare-rerank-inputs.py` does not read `.h5ad` at all -- its
  `features` subcommand takes a Hail table (`--ht`) only. Building an expression
  axis therefore needs its own small script (pick the layer and the group, write
  one row per `var_names` entry with its value in that layer), outside what this
  skill ships, and that script still has to map whatever identifier `var_names`
  holds for the chosen expression source onto the same symbol space as everything
  else (`references/inputs-and-identifiers.md`) -- it is not guaranteed to be a
  symbol, or even consistent, across providers.

## Naming columns uniquely across axes

`FeatureAssembler.assemble` raises if two axes declare the same non-`gene` column
name (`rerank feature axis '{ax}' declares column(s) [...] already provided by axis
'{other}'`) -- the outer join would otherwise silently suffix both to `_x`/`_y`,
feeding the model the same logical quantity twice under two names and breaking the
per-axis grouping the ablation and the within-axis feature selection both rely on
(they match columns by exact name; the case/control architecture audit does NOT --
it reads `n_case_var`/`conc`/`driver_af` by name off the full merged table,
regardless of which axis, if any, declared them). Give every feature table's
columns a name unique across every axis you plan to combine, even axes you have
not written yet.

## Circularity

An axis's own values must not have been derived, even partially, from the SAME
curated resource the labels come from. The canonical trap: a "known-disease-gene"
feature axis built from ClinGen or GenCC, scored against ClinGen- or GenCC-derived
labels -- the axis is then measuring "was this gene already in the database", which
trivially predicts "is this gene in the database" and inflates AUC with nothing
learned. `references/audit-and-leakage.md` covers the one formal tool for this
(`feature_provenance`/`label_provenance`, Python API only) and the informal check
every axis still needs regardless -- read what it is actually built from before
adding it, since the formal tool is opt-in and nothing stops a CLI config from
pairing a circular axis with its label source silently.
