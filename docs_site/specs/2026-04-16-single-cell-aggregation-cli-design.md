# Design spec: single-cell aggregation CLI (`hvantk expression summarize`)

**Date:** 2026-04-16
**Status:** Design approved — ready to plan implementation
**Related planning doc:** `local/planning/2026-04-16-single-cell-aggregation-cli.md`

## Goal

Provide a first-class library + CLI primitive for aggregating a single-cell
expression AnnData into a per-group × per-gene summary, and make the upstream
`mkmatrix ucsc` builder usable on atlases that don't fit in RAM. Analysis
notebooks must consume a pre-aggregated artifact; all streaming/aggregation
optimization lives in hvantk.

## Scope

- Stream-capable `mkmatrix ucsc` (fix `build_ucsc_ad` to handle multi-GB
  gzipped TSV without OOM).
- Reworked `hvantk expression summarize` backed by `scanpy.get.aggregate`.
- Canonical `.h5ad` output for aggregates (groups × genes, aggregations as
  `layers`).
- Retrofit notebook K.1 Cell 3 to consume the aggregated `.h5ad` directly.

**Out of scope for v1:** median aggregation, pickle/parquet/Hail-Table
exporters, multi-modal AnnData, cross-cluster normalization, automatic gene
symbol harmonization, `--backed` mode for atlases > RAM (flagged for v2).

## Architectural principle

Two-step pipeline, mapped onto hvantk's existing Builder/Streamer split:

```
[Builder]                           [Streamer / aggregator]
mkmatrix ucsc   ──►   atlas.h5ad   ──►   expression summarize   ──►   class_mean.h5ad
(source-specific:     (canonical,       (source-agnostic:             (small,
 streams raw UCSC     reusable for       sc.get.aggregate on           reusable across
 TSV → sparse CSR)    describe /         any .h5ad)                    analysis runs)
                      markers /
                      summarize)
```

`atlas.h5ad` is built **once per atlas**. Multiple `summarize` passes (different
`--group-by`, different `--filter-by`, etc.) reuse it. `describe` and `markers`
also consume the same `atlas.h5ad`.

## Why `.h5ad`-only output (not pickle / parquet / HT)

- `sc.get.aggregate` returns an AnnData natively. Writing that directly is
  scverse-native and zero re-implementation.
- All aggregations (`mean`, `sum`, `count_nonzero`, derived `fraction_expressed`)
  travel as `layers` in one file.
- Wide DataFrames and Hail Tables can be rederived in one line each; exporters
  can be added later when a concrete consumer demands them (not v1).

## Why no native scverse reader for UCSC

Research confirmed no scverse tool ships a UCSC Cell Browser reader:

- No `read_ucsc`/`read_cellbrowser` in `scanpy`, `scanpy.external`, or
  `anndata`.
- UCSC's own `cellbrowser` PyPI package writes UCSC format (authoring tool);
  it does not reverse-read.
- `anndata.read_text` / `anndata.read_csv` handle `.gz` but materialize the
  full matrix **dense** via `np.array(list_of_row_arrays, dtype=dtype)` — for
  a 520k × 20k UCSC matrix that is ~41 GB dense in RAM. Strictly worse than
  the current `pd.read_csv` path.
- `anndata.experimental` (`CSRDataset`, `sparse_dataset`, `read_elem`,
  `concat_on_disk`) operates on **existing** h5ad/zarr stores. Useful
  downstream, not at the pre-`.h5ad` boundary.

The streaming gz → CSR → h5ad step is ~30 lines we own. Everything downstream
delegates to `scanpy.get.aggregate`.

## Components

### C1. `build_ucsc_ad` (streaming rewrite)

**File:** `hvantk/tables/ucsc.py::create_anndata_from_ucsc_matrix`
**Caller:** `hvantk/tables/matrix_builders.py::build_ucsc_ad` (one new kwarg)

Replace the single `pd.read_csv` body with a streamed builder:

1. Read header once (via `gzip.open`) to get cell IDs.
2. Loop `pd.read_csv(..., sep=delimiter, chunksize=chunk_size, dtype=np.float32)`
   (works transparently on `.gz`). Default `chunk_size = 500` genes.
3. Each chunk: split off the gene-ID column, optionally split on `|` per
   `split_gene_field`, convert the numeric block to `scipy.sparse.csr_matrix`,
   append to list. Accumulate gene names separately.
4. After loop: `scipy.sparse.vstack([chunks])`, `.T.tocsr()` to get
   `(n_cells × n_genes)` CSR.
5. Wrap in `AnnData(X=csr, obs=..., var=...)`. Existing metadata-join block
   below is unchanged.

**New kwarg on `build_ucsc_ad`:** `chunk_size: int = 500`. Plumb through to
`create_anndata_from_ucsc_matrix`.

**Peak RAM:** one chunk (~500 × 520k × 4B ≈ 1 GB dense) + growing sparse
vstack (sparse; typically < 5 GB for realistic sc atlases).

**Error:** if a chunk fails float32 conversion, raise `ValueError` naming the
offending gene row.

### C2. `summarize_expression_ad` (rewrite on `sc.get.aggregate`)

**File:** `hvantk/utils/matrix_utils.py::summarize_expression_ad`

Replace the hand-rolled per-group loop with:

```python
def summarize_expression_ad(
    adata,
    group_by,                       # str or list[str] — multi-column kept
    filter_by=None,
    min_cells_per_group=10,
):
    if filter_by:
        adata = filter_by_metadata_ad(adata, filter_by)

    by = [group_by] if isinstance(group_by, str) else list(group_by)

    agg = sc.get.aggregate(
        adata,
        by=by,                           # sc.get.aggregate accepts list[str] natively
        func=["mean", "sum", "count_nonzero"],
    )
    # agg: AnnData of shape (n_groups, n_genes), layers["mean"|"sum"|"count_nonzero"].
    # Multi-column groupings produce obs_names joined with "_" (e.g. "IT_LV").

    # sc.get.aggregate does NOT stash per-group cell counts — verified against
    # scanpy 1.10.3. Compute them from the pre-aggregate adata, joining the
    # group-by columns with "_" so the resulting index matches agg.obs_names.
    group_labels = adata.obs[by].astype(str).agg("_".join, axis=1)
    n_cells = group_labels.value_counts().reindex(agg.obs_names).astype(int)

    agg.obs["n_cells"] = n_cells.values
    agg.layers["fraction_expressed"] = (
        agg.layers["count_nonzero"] / n_cells.values[:, None]
    )

    keep = agg.obs["n_cells"].values >= min_cells_per_group
    return agg[keep].copy()
```

**Return type change:** `AnnData` (was `pd.DataFrame`). This is a breaking
change to the library function.

**`count_nonzero` is the scanpy-native primitive**; `fraction_expressed`
is one division. `mean` is the canonical v1 agg. `sum` / `var` /
`count_nonzero` come along for free (one `sc.get.aggregate` pass) and are
available as layers without extra cost.

### C3. `summarize_expression_cmd` (CLI — output change only)

**File:** `hvantk/commands/summarize_expression_cli.py`

Deltas vs. the existing command:

| Parameter | Change |
|---|---|
| `-o, --output` | Suffix must be `.h5ad`. Reject others with `click.BadParameter`. |
| `--min-cells` | Default changed from `50` → `10` (planning-doc alignment). |
| `--group-by` | Keep `multiple=True` (list-of-obs-columns joined to a single grouping key). |
| `--filter-by` | Unchanged. |
| `--overwrite` | Unchanged. |

Output handling: `summary_adata.write_h5ad(output_path)` replaces
`summary_df.to_parquet(output_path)`.

Report card prints: `n_groups`, `n_genes`, min/max cells per group, layer
names, and a table of group labels with cell counts (first 10 + "… and N
more").

### C4. Notebook K.1 retrofit

**File:** `local/notebooks/ptm-eda/notebook_k1_brain_metaatlas_ptm.ipynb`,
Cell 3 (the "Build per-gene × per-Class mean-expression matrix (cached)"
cell).

Collapses from ~45 lines of streaming logic to:

```python
import anndata as ad, subprocess
subprocess.run([
    "hvantk", "expression", "summarize",
    "-m", BRAIN_ATLAS_H5AD,
    "--group-by", CELLTYPE_COL,
    "--min-cells", "10",
    "-o", BRAIN_CACHE_H5AD,
    "--overwrite",
], check=True)
adata = ad.read_h5ad(BRAIN_CACHE_H5AD)
brain_wide = (
    pd.DataFrame(adata.layers["mean"].T, index=adata.var_names, columns=adata.obs_names)
      .rename_axis("gene_symbol")
)
```

Pre-requisite: `BRAIN_ATLAS_H5AD` must exist (one-time build via
`hvantk mkmatrix ucsc -e exprMatrix.tsv.gz -m meta.tsv -o BRAIN_ATLAS_H5AD`).
The old `BRAIN_CACHE` pickle is replaced by `BRAIN_CACHE_H5AD`. Downstream
cells that reference `brain_wide` continue to work unchanged (still a wide
pandas DataFrame with `gene_symbol` index and class names as columns).

## CLI reference (final)

```
hvantk mkmatrix ucsc \
    -e exprMatrix.tsv.gz \
    -m meta.tsv \
    -o atlas.h5ad \
    --chunk-size 500          # new; gene rows per streaming chunk
```

```
hvantk expression summarize \
    -m          atlas.h5ad \                     # required
    --group-by  Class [--group-by Subclass ...] \# required, repeatable
    --filter-by FIELD=VALUE ... \                # optional, repeatable
    --min-cells 10 \                             # optional (default 10)
    -o          class_mean.h5ad \                # required, must end .h5ad
    [--overwrite]
```

## Error handling

| Condition | Behavior |
|---|---|
| `--group-by` column missing from `obs` | `click.BadParameter` listing available obs columns. |
| `--group-by` column has NaN cells | Drop silently; log count. |
| `--filter-by FIELD=VALUE` where FIELD missing | `click.BadParameter`. |
| Filter produces zero cells | `SystemExit(1)` with explicit message. |
| All groups drop out via `--min-cells` | `SystemExit(1)`; report group sizes. |
| Output `.h5ad` exists, no `--overwrite` | Existing behavior (error + exit 1). |
| Output suffix is not `.h5ad` | `click.BadParameter`. |
| Streamed builder hits a non-numeric row | `ValueError` naming the offending gene ID. |

## Testing

Per `MEMORY.md → feedback_no_more_tests.md`, no new tests.

**Update existing tests in place:** `hvantk/tests/test_matrix_utils_anndata.py`
contains `TestSummarizeExpressionAd` which currently asserts a `pd.DataFrame`
return with columns `{gene_id, group, mean, fraction_expressed, n_cells}`.
Migrate those assertions to the new `AnnData` return contract — check
`isinstance(result, ad.AnnData)`, `result.shape == (n_groups, n_genes)`,
`result.obs["n_cells"]`, and the required layers (`mean`, `sum`,
`count_nonzero`, `fraction_expressed`). Any CLI-level tests that asserted
parquet output must migrate similarly.

**Manual validation:** build `atlas.h5ad` from the brain meta-atlas
(520k cells) via `mkmatrix ucsc`; run `expression summarize --group-by Class`;
verify end-to-end completion under 30 min on a laptop; compare the resulting
`layers["mean"]` against the existing K.1 pickle at a handful of (gene,
class) cells for numerical agreement.

If manual validation passes, ship.

## Open implementation details

These are resolved during implementation, not by design:

1. Chunked-build logging cadence (every N genes processed) — cosmetic.
2. Whether to also expose `var` as a first-class per-group aggregation in
   the CLI help (it's already computed into `layers["count_nonzero"]` and
   derivable from `mean` + `sum`; cheap to promote later).

### Verified against scanpy 1.10.3 during design

- `sc.get.aggregate` returns an AnnData with `obs` containing only the
  grouping column(s), not a cell-count column. Per-group cell counts must
  be computed separately from the pre-aggregate `adata.obs`.
- With a list `by=[...]`, `obs_names` are joined with `_`
  (e.g. `"IT_LV"`), which matches the existing hvantk convention.
- Group ordering in the returned AnnData is deterministic (sorted by
  category levels).

## Deliverables (this PR)

1. `create_anndata_from_ucsc_matrix` streaming rewrite + `build_ucsc_ad`
   `chunk_size` kwarg.
2. `summarize_expression_ad` rewrite on `sc.get.aggregate`, AnnData return.
3. `summarize_expression_cmd` output-format change (parquet → `.h5ad`) and
   `--min-cells` default bump.
4. `test_matrix_utils_anndata.py::TestSummarizeExpressionAd` updated to the
   new output contract.
5. Notebook K.1 Cell 3 retrofit.

## Deferred (v2)

- `--backed` mode for atlases > RAM (requires `sc.get.aggregate` verification
  against backed `X`).
- Exporters from aggregated `.h5ad` to wide pickle / parquet / Hail Table
  (add when a concrete consumer requests).
- Retrofit of notebooks F / G / H (optional; they work as-is).
- `sum` / `log2_mean` as first-class `--agg` selectors (the data is already
  present in `layers`, just not surfaced via the CLI).
