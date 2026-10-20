# PTM constraint analysis

`hvantk ptm constraint` compares the population allele frequencies (AF) of PTM variants
with those of non-PTM variants. It stratifies the comparison by the group in which each
gene is most expressed (a tissue, a cell type, or any categorical field of an expression
dataset), by tissue specificity and by gene constraint. It reports group-level
statistics, not per-variant scores.

PTM variants are the site and proximal variants defined by `hvantk ptm annotate`; see
[How variants are classified](ptm.md#how-variants-are-classified). To fit mixed models
on a table you have already split into strata, use `hvantk ptm test`
([Per-stratum mixed models](ptm.md#per-stratum-mixed-models)).

## Quick start

```bash
hvantk ptm constraint \
  --variants-ht data/variants_ptm.ht \
  --expression-source anndata \
  --expression-path data/atlas.h5ad \
  --grouping cell_type \
  --output-dir results/ptm-constraint/
```

The command needs the `constraint` extra: tspex computes tissue specificity, and
matplotlib and seaborn draw the plots ([Installation](../getting-started/installation.md)).
If plotting fails, the tables and `summary.json` are still written.

## Inputs

### Variant table

`--variants-ht` takes the Hail Table written by `hvantk ptm annotate`, or the same rows
as a `.parquet`, `.pkl`, `.csv` or `.tsv` file; a `.csv` or `.tsv` may be compressed
(`.gz` or `.bgz`).

| Column | Set with (default) | Required |
| --- | --- | --- |
| `is_ptm_site`, `is_ptm_proximal` | fixed names | Yes; a tabular file may omit `is_ptm_site` |
| Gene identifier | `--gene-field` (`gene_symbol`) | Yes |
| AF | `--af-field` (`gnomad_af_genomes`) | Yes, unless `--include-zero-af`; then a missing column counts as AF 0 |
| Label | `--label-field` (`rf_label`) | Yes, unless `--label-filter all` |
| LOEUF | `--loeuf-field` (`loeuf`) | No; without it the factorial table is empty |
| PTM categories | `--ptm-category-field` (`ptm_types`) | No; without it the category table is empty |

`annotate` adds only the PTM fields. Add the gene, AF, LOEUF and label columns to the
variant table yourself and pass their names. `--label-filter` keeps the rows whose label
is exactly `TN` (benign, the default) or `TP` (pathogenic), the values the ClinVar
training-set labelling writes to `rf_label`; pass `all` to keep every row.

### Expression data

| `--expression-source` | `--expression-path` | `--grouping` | `--expression-metric` | `--min-cells-per-group` counts |
| --- | --- | --- | --- | --- |
| `hail-mt` | MatrixTable directory; values in the entry field `x`, genes in the first row key | A column field | `mean` only | Samples |
| `anndata` | `.h5ad` file with gene identifiers as `var_names`; values from `X` | An `obs` column | `median`, `mean` or `median_nonzero` | Cells |
| `tabular` | Gene × group matrix as `.parquet` or `.pq`, `.pkl` or `.pickle`, `.tsv` or `.txt`, or `.csv`; genes in the index (the first column of a text file) | Required but ignored: the columns are the groups | Ignored | Ignored |

- With `hail-mt`, pass `--expression-metric mean`. The default, `median`, stops with an
  error.
- With `anndata`, `mean` aggregates through scanpy, which the `expression` extra
  installs. `median` and `median_nonzero` (the median of the positive values) do not
  need it. An AnnData keyed by transcript (`var.index.name` is `transcript_id`) is
  rejected. `hvantk expression describe -m data/atlas.h5ad` lists the `obs` fields you
  can group by.
- Values are used as stored, without a log transform. Non-numeric values become missing,
  genes with no value are dropped, the remaining missing values become 0, and negative
  values are clipped to 0. The matrix must have between 2 and 200 groups.

### Matching genes and groups

- Gene identifiers are compared after removing version suffixes (`ENSG00000141510.17`
  becomes `ENSG00000141510`) on both sides. When the two inputs use different
  identifiers, pass `--gene-id-mapping`: a TSV with a header row, the expression
  identifier without its version suffix in the first column and the variant table's
  identifier in the second.
  Identifiers missing from the file are kept as they are. Genes that end up with the
  same identifier are merged, keeping the highest value per group.
- `--group-mapping` collapses group labels, for example GTEx tissues into organ systems:
  a TSV with a header row, the raw label in the first column and the analysis group in
  the second. It changes only each gene's primary group; tissue specificity and the
  expression filter still use the original groups. Labels missing from the file are
  kept.

## Method

1. Load the variants and apply `--label-filter`. A variant is a PTM variant if
   `is_ptm_site` or `is_ptm_proximal` is true. A missing AF counts as 0.
2. Aggregate expression into a gene × group matrix. Groups with fewer than
   `--min-cells-per-group` samples or cells are dropped.
3. For each gene, compute τ, the tissue-specificity index of Yanai et al. (2005), with
   [tspex](https://pypi.org/project/tspex/); its primary group, the group with the
   highest value (after `--group-mapping`); and that highest value.
4. Keep the variants whose gene is in the matrix with a highest value of at least
   `--expressed-threshold`. With `--af-observed-only` (the default), keep only variants
   with an AF above 0.
5. Each variant takes its gene's primary group, τ and LOEUF, so it belongs to exactly
   one group.

The first four analyses compare a set of PTM variants with a set of non-PTM variants in
the same way:

- `ratio` is (median non-PTM AF + 1e-8) / (median PTM AF + 1e-8), and `log2_ratio` its
  base-2 logarithm. A ratio above 1 (log2 above 0) means the PTM variants have lower AF.
  `mean_ratio` and `log2_mean_ratio` repeat this with means.
- `p_value` comes from a one-sided Mann-Whitney U test that non-PTM AFs are greater than
  PTM AFs. It needs at least 3 variants on each side and is empty otherwise.
- `n_ptm`, `n_non`, `med_ptm`, `med_non`, `mean_ptm` and `mean_non` give the counts,
  medians and means.

| Analysis (file) | PTM and non-PTM variants compared within | A row is skipped when |
| --- | --- | --- |
| Per-group ranking (`per_group_ranking.tsv`), sorted by `log2_ratio`, largest first | Each primary group | The group has fewer than `--min-variants-per-group` variants |
| τ quartiles (`tau_quartile.tsv`) | Each quartile of τ over the variants: `Q1_housekeeping`, `Q2`, `Q3`, `Q4_tissue_specific` | The quartile has fewer than `--min-variants-per-group` variants. The table is empty if τ does not split into four bins |
| τ × LOEUF factorial (`loeuf_factorial.tsv`) | Each cell of two median splits over the variants: τ `broad` or `tissue_specific`, and LOEUF `constrained` (lower half) or `unconstrained` | The cell has fewer than `--min-variants-per-group` variants |
| PTM category × group (`category_group_heatmap.tsv`) | Each group, comparing the PTM variants of one category with all non-PTM variants of the group | The group has fewer than 3 non-PTM variants, or fewer than half of `--min-variants-per-group`; or the category has fewer than 3 PTM variants |

In the category analysis, categories are lowercased, a variant counts in every category
of its set, and PTM variants without a category are left out. In a tabular input,
set-like strings such as `{'phosphorylation'}` or `phosphorylation;acetylation` are
split into categories.

The fifth analysis, the within-gene comparison (`within_gene_paired.tsv`), takes every
gene with at least one PTM and one non-PTM variant and computes `delta`, the median
non-PTM AF minus the median PTM AF. A one-sided Wilcoxon signed-rank test asks whether
the deltas are above 0; it needs at least 3 genes. Its result is in `summary.json` under
`within_gene_paired`: `n_genes`, `statistic`, `p_value`, `pct_genes_depleted` (the
fraction of genes with a `delta` above 0, between 0 and 1) and `median_delta`.

No p-value is adjusted for multiple testing. Quartiles and median splits are taken over
variants, not genes, so genes with many variants weigh more.

## Outputs

```text
<output-dir>/
├── summary.json                    # counts, top 10 groups, rows of the other tables, within-gene test, configuration
├── per_group_ranking.tsv
├── tau_quartile.tsv
├── loeuf_factorial.tsv
├── category_group_heatmap.tsv      # long form: one row per group and category
├── within_gene_paired.tsv          # gene, n_ptm, n_non, med_ptm, med_non, delta
├── plots/
│   ├── per_group_ranking.png       # log2 ratio of up to 25 groups
│   ├── tau_quartile.png
│   ├── loeuf_factorial.png         # needs seaborn
│   └── category_group_heatmap.png  # needs seaborn
└── report.html                     # plots embedded, with the top rows of each table
```

The first four tables start with their stratum columns (`group`; `tau_quartile`;
`tau_bin` and `loeuf_bin`; `group` and `category`), followed by the comparison columns
above. A test with no rows writes an empty TSV, and its plot is skipped. The command
prints one summary line: the number of variants, groups and genes, and the top group's
ratio.

## Options

| Option | Default | Description |
| --- | --- | --- |
| `--variants-ht` | required | Annotated variants: a Hail Table directory or a tabular file |
| `--expression-source` | required | `hail-mt`, `anndata` or `tabular` |
| `--expression-path` | required | Expression data |
| `--grouping` | required | Group field: a column field (`hail-mt`) or an `obs` column (`anndata`); ignored for `tabular` |
| `-o`, `--output-dir` | required | Output directory |
| `--label-filter` | `TN` | `TN`, `TP` or `all` |
| `--label-field` | `rf_label` | Label column |
| `--gene-field` | `gene_symbol` | Gene identifier column |
| `--af-field` | `gnomad_af_genomes` | AF column |
| `--loeuf-field` | `loeuf` | LOEUF column |
| `--ptm-category-field` | `ptm_types` | PTM category column |
| `--gene-id-mapping` | none | TSV mapping expression gene identifiers to the variant table's |
| `--group-mapping` | none | TSV mapping group labels to analysis groups |
| `--expression-metric` | `median` | `median`, `mean` or `median_nonzero` |
| `--min-cells-per-group` | 50 | Minimum samples or cells per group (`hail-mt` and `anndata`) |
| `--min-variants-per-group` | 20 | Minimum variants per group, quartile or factorial cell |
| `--expressed-threshold` | 1.0 | Minimum highest group value for a gene to be kept |
| `--af-observed-only` / `--include-zero-af` | `--af-observed-only` | Keep only variants with an AF above 0, or also those with an AF of 0 |
| `--flanking-codons` | 7 | Shown in the report header only; the window was set when you ran `annotate` |
| `--overwrite` | off | Has no effect; existing outputs are always replaced |

The choices are case-insensitive.

## Examples

```bash
# Bulk tissues from a GTEx-style MatrixTable (hail-mt aggregates by mean only)
hvantk ptm constraint \
  --variants-ht data/variants_ptm.ht \
  --expression-source hail-mt \
  --expression-path data/gtex_tpm.mt \
  --grouping SMTSD \
  --expression-metric mean \
  --output-dir results/ptm-gtex/

# A pre-aggregated gene x group matrix keyed by Ensembl ID; --grouping is not used
hvantk ptm constraint \
  --variants-ht data/variants_ptm.parquet \
  --expression-source tabular \
  --expression-path data/gene_by_celltype.tsv \
  --grouping cell_type \
  --gene-id-mapping data/ensembl_to_symbol.tsv \
  --output-dir results/ptm-custom/
```

## Caveats

- Nothing adjusts for mutation rate. Codons of modified residues may mutate at different
  baseline rates from other codons.
- Each gene belongs to one group, the one where it is most expressed, so a broadly
  expressed gene contributes to that group only.
- The gene identifiers of the two inputs must match. A join that matches no gene stops
  with an error that suggests checking symbols against Ensembl IDs.

## Python API

```python
from hvantk.algorithms.ptm import PTMConstraintConfig, run_ptm_constraint

config = PTMConstraintConfig(
    variants_ht_path="data/variants_ptm.ht",
    expression_source="anndata",  # lower case, as on the command line
    expression_path="data/atlas.h5ad",
    grouping="cell_type",
    output_dir="results/ptm-constraint/",
    label_filter="TN",  # exactly "TN", "TP" or "all"
)
result = run_ptm_constraint(config)  # validates the config and writes every output
print(result.summary())
```

## See also

- [PTM](ptm.md): the sites table, the classification and `hvantk ptm test`
- [PTM example](../examples/ptm.md): building and annotating the inputs
- [Expression analysis](expression.md): inspecting and summarizing AnnData files
- [Installation](../getting-started/installation.md): the `constraint` and `expression` extras
