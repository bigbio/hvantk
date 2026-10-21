# Expression analysis

`hvantk expression` inspects, summarizes, and extracts marker genes from expression
`AnnData` (`.h5ad`) files -- built by the Expression Atlas, UCSC Cell Browser, or CPTAC
plugins, or any other gene-keyed AnnData.

## When to use it

- `describe` -- a quick look at an AnnData's observation metadata (what you can group
  or filter by) before running anything else.
- `summarize` -- collapse cells/samples into a per-group x per-gene pseudo-bulk
  AnnData, e.g. as a `matrix:` source for [`hvantk annotate prepare`](annotate.md).
- `markers` -- per-group differential expression, written as a `GeneSetCollection`
  (see [Gene sets](genesets.md)) for `hvantk psroc` / `hvantk enrichex`.
- `summarize-ucsc` -- the same output as `summarize`, fused with UCSC's own streaming
  reader so a raw `exprMatrix.tsv(.gz)` + metadata TSV never materializes an
  intermediate cells x genes `.h5ad`.

## Inputs

An expression AnnData with `obs` = cells/samples and `var` = genes, and `--group-by` /
`--filter-by` naming columns that exist in `obs`. `summarize` and `markers` require
`var` to be gene-keyed (`adata.var.index.name` must not be `"transcript_id"`); `describe`
only inspects, so it never refuses -- it just labels the variable axis "transcripts"
instead of "genes" for a transcript-keyed file. `summarize-ucsc` takes the UCSC raw
files directly instead: `-e` (`exprMatrix.tsv[.gz]`) and `-m` (the cell metadata TSV).

### The `expression` extra

`describe` and `summarize-ucsc` work on a base install. `summarize` and `markers` need
scanpy -- install the `expression` extra (see
[Installation](../getting-started/installation.md#optional-features-extras)) -- and
**it does not install on Intel macOS** (llvmlite ships no x86_64 macOS wheel); run them
on Linux, Apple Silicon, or the cluster instead.

## Worked example

```bash
# Inspect what you can group or filter by
hvantk expression describe -m data/heart_sc.h5ad

# Collapse into a per-cell-type summary
hvantk expression summarize \
  -m data/heart_sc.h5ad \
  --group-by cell_type \
  -o out/heart_celltype_summary.h5ad

# Extract Wilcoxon marker genes per cell type
hvantk expression markers \
  -m data/heart_sc.h5ad \
  --group-by cell_type \
  --top-n 200 \
  -o out/heart_markers.json

# Same shape as `summarize`, streamed straight from UCSC's raw files
hvantk expression summarize-ucsc \
  -e data/exprMatrix.tsv.gz -m data/meta.tsv \
  --group-by Class \
  -o out/class_summary.h5ad
```

## Options reference

### `expression describe`

| Option | Required | Description |
| --- | --- | --- |
| `-m, --matrix-table PATH` | yes | Path to an expression AnnData (`.h5ad`). |

Prints `n_obs`/`n_vars` and, per `obs` column, its dtype plus (for a categorical
column) its level count or (for a numeric column) its min/max. Uses a pre-computed
`adata.uns["column_summary"]` when present (instant); otherwise computes it on the fly.

### `expression summarize`

| Option | Required | Description |
| --- | --- | --- |
| `-m, --matrix-table PATH` | yes | Path to an expression AnnData (`.h5ad`). |
| `--group-by TEXT` | yes, repeatable | `obs` field(s) to group by. |
| `--filter-by TEXT` | no, repeatable | `FIELD=VALUE`; a cell/sample is kept only if `obs[FIELD]` equals (one of) the given value(s) -- exact match, not substring. |
| `--min-cells INTEGER` | no | Drop groups with fewer cells. Default `10`. |
| `-o, --output PATH` | yes | Must end in `.h5ad`. |
| `--overwrite` | no | Without it, an existing `--output` is a hard error. |

### `expression markers`

| Option | Required | Description |
| --- | --- | --- |
| `-m, --matrix-table PATH` | yes | Path to an expression AnnData (`.h5ad`). |
| `--group-by TEXT` | yes | A single `obs` field for differential expression. |
| `--filter-by TEXT` | no, repeatable | `FIELD=VALUE`, as above. |
| `--method [wilcoxon\|t-test\|t-test_overestim_var\|logreg]` | no | Passed to `scanpy.tl.rank_genes_groups`. Default `wilcoxon`. |
| `--top-n INTEGER` | no | Maximum markers per group. Default `200`. |
| `-o, --output PATH` | yes | `.json` or `.gmt`. |
| `--results-tsv PATH` | no | Also save the full per-gene results table. |
| `--overwrite` | no | Without it, an existing `--output` is a hard error. |

### `expression summarize-ucsc`

| Option | Required | Description |
| --- | --- | --- |
| `-e, --expression-matrix PATH` | yes | UCSC expression TSV (plain or gzipped). |
| `-m, --metadata PATH` | yes | UCSC cell metadata TSV. |
| `--group-by TEXT` | yes, repeatable | Metadata column(s) to group by. |
| `--filter-by TEXT` | no, repeatable | `FIELD=VALUE`, as above (against the metadata). |
| `--min-cells INTEGER` | no | Drop groups with fewer cells. Default `10`. |
| `--gene-column TEXT` | no | Default `gene`. |
| `--delimiter TEXT` | no | Default tab. |
| `--split-gene-field` / `--no-split-gene-field` | no | Default on: split a pipe-delimited UCSC gene label (e.g. `A\|B`) into one row per gene. |
| `-o, --output PATH` | yes | Must end in `.h5ad`. |
| `--overwrite` | no | Without it, an existing `--output` is a hard error. |

## Outputs

**`summarize` / `summarize-ucsc`** write an AnnData shaped `(n_groups, n_genes)` with
layers `mean`, `sum`, `count_nonzero`, `fraction_expressed`, and `obs["n_cells"]`
holding the per-group cell count. Multi-field `--group-by` joins group labels with
`_`. `summarize-ucsc`'s output has the same schema, so downstream consumers don't need
to branch on which command produced it.

**`markers`** writes a `GeneSetCollection` (one `GeneSet` per group, the top `--top-n`
genes by the chosen method) as `.json` or `.gmt` -- see
[Gene sets](genesets.md#the-genesetcollection-json-format) for the format. With
`--results-tsv`, the full `scanpy.get.rank_genes_groups_df` table is also written.

## Python API

```python
from hvantk.core.io.anndata_io import load_anndata
from hvantk.algorithms.expression.matrix_utils import describe_expression_ad, summarize_expression_ad

adata = load_anndata("data/heart_sc.h5ad")
info = describe_expression_ad(adata)  # {"n_obs", "n_vars", "fields": [...]}

summary = summarize_expression_ad(
    adata, group_by=["cell_type"], filter_by={"region": "LV"}, min_cells_per_group=10,
)  # raises TranscriptKeyedError / ImportError exactly as the CLI does
```

`filter_by_metadata_ad(adata, filters)` (same module) is the standalone pre-filter
`markers` applies before ranking. `require_scanpy()` returns the `scanpy` module or
raises `ImportError` naming the `expression` extra; `summarize_expression_ad` and the
CLI both call it, so the message is identical either way.

## Caveats

- `--filter-by`/`--group-by` values are matched exactly against `obs`/metadata values
  as strings; there is no substring or regex mode.
- `summarize`/`summarize-ucsc` drop any cell with a missing value in a `--group-by`
  column before aggregating (logged, not an error) -- otherwise it would form a
  literal `"nan"` group.
- A transcript-keyed AnnData (from an Expression Atlas build without
  `--plugin-arg level=gene`) is refused by `summarize` and `markers`, with the rebuild
  command in the error message; `describe` still works, labeling the axis
  "transcripts". See [Data Sources](../guide/data-sources.md) for the
  `expression-atlas:dataset` acquisition details.
- `--min-cells` is applied after aggregation, not as a pre-filter on raw cells.

## See also

- [Building datasets](../guide/usage.md) -- build the expression datasets this reads
- [Data Sources](../guide/data-sources.md) -- Expression Atlas, UCSC Cell Browser, CPTAC acquisition
- [Gene sets](genesets.md) -- the `GeneSetCollection` format `markers` writes
- [Installation](../getting-started/installation.md#optional-features-extras) -- the `expression` extra and its Intel-macOS caveat
