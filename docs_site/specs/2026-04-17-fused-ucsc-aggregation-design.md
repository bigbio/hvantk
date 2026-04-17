# Design spec: backed-write `mkmatrix ucsc` + fused `expression summarize-ucsc`

**Date:** 2026-04-17
**Status:** Design approved — ready to plan implementation
**Supersedes parts of:** `docs_site/specs/2026-04-16-single-cell-aggregation-cli-design.md`
(The prior spec assumed the atlas was ~5 GB in RAM. Manual validation on the
real 122k × 520k adult-cortex meta-atlas showed the in-memory `vstack` +
`.T.tocsr()` path is a wall-clock bottleneck (~2–3 h) and the atlas artifact
itself reaches ~20–40 GB on disk. This spec adjusts both stages and adds a
fused shortcut for the common `matrix → aggregate → summary` workflow.)

## Goal

Make single-cell atlas ingestion and aggregation scale to 500k+-cell UCSC
atlases on a laptop. Two complementary execution paths, user-selected per
invocation:

1. **Build-and-keep** — incremental disk write of a canonical `atlas.h5ad`
   (`obs=cells`, `var=genes`) for downstream reuse (multiple aggregations,
   `describe`/`markers`, collaborator share).
2. **Fused shortcut** — one-pass stream + aggregate that emits only the tiny
   `summary.h5ad`, never materializing the cells × genes matrix.

## Context

- Issue #93 (closed by this branch) replaced `pd.read_csv(chunksize=...)`
  with a row-oriented `gzip + np.fromstring` streamer inside
  `create_anndata_from_ucsc_matrix`. That fixed the parser failure on large
  gzipped inputs and gave ~6× speedup over pandas/pyarrow for wide (520k-col)
  atlases. It did **not** change the fact that the function still builds the
  full sparse matrix in RAM before writing `.h5ad`.
- Manual validation run on `adult-ctx-meta-atlas`
  (`exprMatrix.tsv.gz`, 8.7 GB, 520,014 cells × 122,697 gene/transcript rows):
  - Streaming parser: ~78 ms/row steady-state → projected ~160 min to parse
  - In-memory `vstack(chunks).T.tocsr()`: unknown, projected multi-hour, RAM
    ceiling unclear (accumulated sparse chunks + transpose copy ≳ 40 GB)
  - User aborted at ~40,000 rows processed due to projected runtime.
- User workflow confirmed: atlases of this size will be routinely processed.
  Common pipeline is `full expression matrix → aggregate → summary stats /
  joins`. Keeping `atlas.h5ad` is valuable when (a) running multiple
  aggregations from the same input, or (b) sharing the raw matrix with a
  collaborator.

## Scope

**In scope (this PR)**

- `mkmatrix ucsc` rewrite: stream rows and write `X` as a backed CSC dataset
  column-by-column (one gene = one column). Output `atlas.h5ad` shape
  `(n_cells, n_genes)`, AnnData standard convention.
- New `expression summarize-ucsc` CLI: takes raw UCSC files + `--group-by`,
  streams + aggregates in a single pass via `np.bincount`, emits
  `summary.h5ad` directly. Schema matches what the existing
  `expression summarize` emits, so downstream consumers cannot tell which
  path produced a given summary.
- Refactor: extract `_iter_ucsc_rows()` from `create_anndata_from_ucsc_matrix`
  as the single row-streaming primitive shared by both paths. The in-memory
  `create_anndata_from_ucsc_matrix` remains for small inputs and tests.
- `build_ucsc_ad` gains a `backed: bool` kwarg (default: auto-select by
  input size).
- `mkmatrix ucsc` CLI gains `--backed/--no-backed` with default
  `--backed=auto` (backed if input `.gz` > 1 GB, else in-memory).

**Out of scope (deferred)**

- `expression summarize` gaining a `--backed` mode for reading large
  `atlas.h5ad` files on disk. Required eventually — otherwise running stage-2
  alone on a 40 GB atlas OOMs the same way. Tracked for a follow-up PR; for
  now, users of large atlases should prefer `summarize-ucsc` or use the fused
  path.
- Retrofit of `describe` / `markers` to handle backed atlases.
- Retrofit of non-UCSC expression sources (Expression Atlas, CPTAC).

## Architectural principle

One low-level streamer, two accumulation strategies, three CLI entry points.

```
                      UCSC exprMatrix.tsv.gz + meta.tsv
                                   │
                                   ▼
              ┌──────────────────────────────────────────┐
              │ _iter_ucsc_rows()                        │  ← row streamer
              │ yields (gene_name, float32[n_cells])     │    (already written
              └───────────┬──────────────────────┬───────┘     for #93)
                          │                      │
          ┌───────────────┘                      └──────────────┐
          ▼                                                     ▼
 mkmatrix ucsc                              expression summarize-ucsc
 ─────────────                              ────────────────────────
 append CSC column                          np.bincount(group_idx, weights=row)
 (one gene per step)                        into (groups × genes) accumulators
          │                                                     │
          ▼                                                     ▼
     atlas.h5ad                                            summary.h5ad
 (cells × genes sparse CSC)                         (groups × genes, layers)

         atlas.h5ad ──► expression summarize (existing) ──► summary.h5ad
```

The CSC orientation is not incidental: it matches the parse order (one gene
per append) so the backed builder never has to buffer rows or transpose.

## Components

### C1. `_iter_ucsc_rows` (extract as helper)

**File:** `hvantk/tables/ucsc.py`
**New function:**

```python
def _iter_ucsc_rows(
    expression_matrix_path: str,
    delimiter: str = "\t",
    split_gene_field: bool = True,
) -> tuple[list[str], Iterator[tuple[str, np.ndarray]]]:
    """Open a UCSC expression TSV (plain or gzipped) and return
    ``(cell_ids, row_iterator)``. ``cell_ids`` is the header's cell-column
    list (computed once, returned eagerly). ``row_iterator`` yields
    ``(gene_name, row_values_float32)`` per data line; the underlying file
    handle is closed when the iterator is exhausted.
    """
```

Behavior already implemented in the current `create_anndata_from_ucsc_matrix`
(gzip + `np.fromstring` + shape validation). This refactor just lifts the
parse loop into its own generator so the backed builder and the fused
aggregator share it byte-for-byte. Returning `cell_ids` separately (rather
than on every yield) avoids passing the same 520k-element list through every
row.

The existing `create_anndata_from_ucsc_matrix` becomes a thin wrapper:
iterate `_iter_ucsc_rows`, accumulate sparse chunks, vstack, transpose, wrap
in AnnData. Contract unchanged.

### C2. `build_ucsc_atlas_backed` (new backed builder)

**File:** `hvantk/tables/ucsc.py`
**Function:**

```python
def build_ucsc_atlas_backed(
    expression_matrix_path: str,
    output_path: str,
    metadata_df: pd.DataFrame | None = None,
    gene_column: str = UCSC_GENE_COLUMN,
    delimiter: str = "\t",
    split_gene_field: bool = True,
    column_batch: int = 64,
    overwrite: bool = False,
) -> str:
    """Stream-build an AnnData .h5ad file on disk, appending one gene
    (or batch of genes) at a time as CSC columns.

    Never materializes the full cells × genes matrix in RAM. Peak RSS is
    roughly ``column_batch × n_cells × 4 B`` for the current parse buffer
    plus HDF5 I/O overhead.
    """
```

Implementation outline:

1. Open `output_path` as HDF5 with `h5py.File(..., "w")`.
2. Write `obs`, `var` placeholders (shape from header + metadata).
3. Create the `X` group with
   `CSCDataset(..., dtype=np.float32, indptr_dtype=np.int64)` —
   `indptr_dtype=np.int64` is required: at ~5% density on 520k × 122k we hit
   ~3.2B nonzeros, exceeding the int32 default.
4. Iterate `_iter_ucsc_rows`; every `column_batch` rows, convert the batch
   into a `scipy.sparse.csc_matrix` of shape `(n_cells, column_batch)`,
   append via `csc_dataset.append(batch)`.
5. After stream end, write final `var` (gene names), `obs` (cells + metadata
   join), and `uns["hvantk_metadata"]` provenance.

Atlas shape on disk: `(n_cells, n_genes)`, `X` stored as CSC, standard
AnnData axis convention. Scanpy and downstream scverse tools handle CSC
transparently; operations that prefer CSR convert on-demand at read time.

### C3. `summarize_ucsc_streaming` (new fused aggregator)

**File:** `hvantk/tables/ucsc.py` (or a new `hvantk/tables/ucsc_aggregate.py`
if the file gets too large)
**Function:**

```python
def summarize_ucsc_streaming(
    expression_matrix_path: str,
    metadata_df: pd.DataFrame,
    group_by: str | list[str],
    filter_by: dict[str, str | list[str]] | None = None,
    min_cells_per_group: int = 10,
    delimiter: str = "\t",
    split_gene_field: bool = True,
) -> ad.AnnData:
    """Stream the UCSC expression matrix and aggregate per-group per-gene
    statistics in one pass. Returns a small AnnData (groups × genes) with
    layers ``mean``, ``sum``, ``count_nonzero``, ``fraction_expressed``.
    """
```

Implementation outline:

1. Apply `filter_by` to `metadata_df`. The existing helper
   `filter_by_metadata_ad` operates on an `AnnData`; the fused path only has
   a DataFrame at this point, so either generalize that helper to accept a
   DataFrame or factor out a small `_filter_metadata_df(df, filter_by)` that
   both callers use. Decide during implementation — contract is the same
   either way (`FIELD=VALUE` equality; repeated keys form an OR; multiple
   keys form an AND).
2. Build composite group label per cell: `metadata_df[group_by].astype(str)
   .agg("_".join, axis=1)` — matches `scanpy.get.aggregate`'s `obs_names`
   convention.
3. Factorize the labels → `group_idx: np.ndarray[int32]`, `n_groups`.
4. Parse the expression header via `_iter_ucsc_rows`, reindex `group_idx` and
   `metadata_df` to the expression-matrix cell order.
5. Pre-allocate accumulators:
   - `sum_matrix: float64[n_groups, n_genes]`
   - `count_nz_matrix: int64[n_groups, n_genes]`
   - `gene_names: list[str]` (grown lazily — not a tight loop)
   - `n_cells_per_group: int64[n_groups]` (from `np.bincount(group_idx)`)
6. For each `(gene_name, row)` from `_iter_ucsc_rows`, compute per-gene column
   `j`:
   ```python
   sum_matrix[:, j] = np.bincount(group_idx, weights=row, minlength=n_groups)
   count_nz_matrix[:, j] = np.bincount(
       group_idx, weights=(row != 0).astype(np.float64), minlength=n_groups
   )
   ```
7. After stream end:
   ```python
   mean            = sum_matrix / n_cells_per_group[:, None]
   fraction_exprd  = count_nz_matrix / n_cells_per_group[:, None]
   ```
8. Drop groups with `n_cells_per_group < min_cells_per_group`.
9. Build output AnnData:
   - `obs_names = composite labels`
   - `obs` holds one column per field in `group_by`, plus `n_cells`
   - `var_names = gene_names`
   - `X = mean` (float32)
   - `layers = {sum, count_nonzero, fraction_expressed}`

Memory footprint: `(n_groups + 2·n_groups) × n_genes × 8 B` accumulators +
one `n_cells × 4 B` row buffer. At 20 groups × 122k genes × 520k cells:
~60 MB accumulators + 2 MB per-row buffer. Independent of `n_cells`.

### C4. `build_ucsc_ad` (wire through)

**File:** `hvantk/tables/matrix_builders.py`

Add a `backed: bool | None = None` kwarg. When `None`, auto-select: backed
if input file size > `BACKED_BUILDER_THRESHOLD_BYTES` (constant, default
1 GB); else in-memory. Dispatches to `build_ucsc_atlas_backed` or
`create_anndata_from_ucsc_matrix` accordingly.

### C5. `mkmatrix ucsc` CLI (existing)

**File:** `hvantk/commands/make_matrix_cli.py`

Add `--backed / --no-backed` toggle (default omitted = auto).

### C6. `expression summarize-ucsc` CLI (new)

**File:** `hvantk/commands/summarize_expression_cli.py`

New subcommand under the existing `expression_group`. Flags:

```
-e, --expression-matrix PATH    # required
-m, --metadata          PATH    # required
--group-by              STR     # required, repeatable
--filter-by             STR     # optional, repeatable, FIELD=VALUE
--min-cells             INT     # default 10
-o, --output            PATH    # required, must end .h5ad
--overwrite                     # flag
--gene-column           STR     # default "gene"
--delimiter             STR     # default "\t"
--split-gene-field      FLAG    # default True
```

Body:

1. Validate output suffix; check `--overwrite` / existence.
2. Load `metadata_df` via `load_ucsc_metadata`.
3. Parse `--filter-by FIELD=VALUE` pairs; apply to `metadata_df`.
4. Validate every `--group-by` column is in `metadata_df.columns`;
   `click.BadParameter` with available columns if not.
5. Call `summarize_ucsc_streaming(...)`.
6. Write result via `adata.write_h5ad(output)`.
7. Print the same report card as `expression summarize` for UX parity.

## CLI reference (final)

```
# Backed build: user will reuse atlas.h5ad downstream
hvantk mkmatrix ucsc \
    -e exprMatrix.tsv.gz \
    -m meta.tsv \
    -o atlas.h5ad \
    [--backed | --no-backed]       # default: auto (backed if .gz > 1 GB)
    [--chunk-size 500]

# Existing two-step aggregator (unchanged)
hvantk expression summarize \
    -m atlas.h5ad \
    --group-by Class [--group-by Subclass ...] \
    [--filter-by FIELD=VALUE ...] \
    [--min-cells 10] \
    -o summary.h5ad [--overwrite]

# Fused one-pass aggregator (new)
hvantk expression summarize-ucsc \
    -e exprMatrix.tsv.gz \
    -m meta.tsv \
    --group-by Class [--group-by Subclass ...] \
    [--filter-by FIELD=VALUE ...] \
    [--min-cells 10] \
    -o summary.h5ad [--overwrite]
```

Both `summarize` commands emit identical `.h5ad` schema. Downstream consumers
(notebook K.1, etc.) do not branch on provenance.

## Data flow — backed builder

```
exprMatrix.tsv.gz ── _iter_ucsc_rows ──► (gene_name, float32[n_cells])
                                              │
                                              ▼  batch of column_batch genes
                          csc_matrix of shape (n_cells, column_batch)
                                              │
                                              ▼  .append(block)
                               atlas.h5ad :: X  (CSCDataset, indptr=int64)
                                              │
                    meta.tsv ──► reindex to cell_ids ──► atlas.h5ad :: obs
                                                         atlas.h5ad :: var
                                                         atlas.h5ad :: uns
```

## Data flow — fused aggregator

```
meta.tsv ──► obs DataFrame ──► pick --group-by column(s), apply --filter-by
                                              │
                                              ▼
                       group_idx: np.ndarray[int32], shape=(n_cells,)
                       n_groups from np.unique(group_idx)

                       pre-allocate:
                         sum[n_groups, n_genes]        float64
                         count_nz[n_groups, n_genes]   int64
                         n_cells_per_group[n_groups]   int64

exprMatrix.tsv.gz ── _iter_ucsc_rows ──► (gene_name, float32[n_cells])
                                              │
                                              ▼  per row (gene j):
             sum[:, j]      = np.bincount(group_idx, weights=row, minlength=n_groups)
             count_nz[:, j] = np.bincount(group_idx, weights=(row != 0), minlength=n_groups)
                                              │
                                              ▼  after last row:
                        mean            = sum / n_cells_per_group[:, None]
                        fraction_exprd  = count_nz / n_cells_per_group[:, None]
                        drop groups with n_cells < --min-cells
                                              │
                                              ▼
                        AnnData(X=mean, layers={sum, count_nonzero,
                                                 fraction_expressed},
                                obs includes n_cells + group_by columns)
                                              │
                                              ▼
                                         summary.h5ad
```

## Error handling

| Condition | Path | Behavior |
|---|---|---|
| Input `.tsv.gz` / `.tsv` missing | both | `FileNotFoundError` naming the path. |
| Metadata `.tsv` missing | both | `FileNotFoundError` naming the path. |
| Row has wrong column count or non-numeric token | both | `ValueError` naming the offending gene (inherited from #93 fix). |
| `--group-by` column not in metadata columns | fused | `click.BadParameter` listing available obs columns. |
| `--filter-by FIELD=VALUE` with FIELD missing | fused | `click.BadParameter` listing available obs columns. |
| `--filter-by` produces zero cells | fused | `SystemExit(1)` with "no cells match filter". |
| Every group falls below `--min-cells` | fused | `SystemExit(1)` with a table of observed group sizes. |
| Output `.h5ad` exists, no `--overwrite` | both | Error + `SystemExit(1)`. |
| Output suffix ≠ `.h5ad` | both | `click.BadParameter`. |
| Metadata cell ids don't overlap expression header | both | **`ValueError`** — flip from the current warn-only behavior. A backed build or stream-aggregate that would produce a zero-obs output should fail at step 0. |
| Cell with NaN in `--group-by` column | fused | Drop silently; log `count` plus one example cell id. No imputation (belongs upstream — it is a domain-modeling choice, not a reduce operation). |
| Append would exceed CSC `indptr` int32 cap | builder | Non-issue: builder passes `indptr_dtype=np.int64` up front. If a future change drops that, `anndata` raises `OverflowError` with a clear message. |

## Testing

Per `MEMORY.md → feedback_no_more_tests.md`: **no new test files**.

**Regression harness (unchanged tests that must keep passing):**

| Test | Covers |
|---|---|
| `hvantk/tests/test_expression_builders_anndata.py::TestBuildUcscAd` | `create_anndata_from_ucsc_matrix` contract. Guards the in-memory path and (via the extracted helper) the parse loop shared with backed/fused paths. |
| `hvantk/tests/test_matrix_utils_anndata.py::TestSummarizeExpressionAd` | `summarize_expression_ad` contract. Pins the `.h5ad` schema that `summarize-ucsc` must match. |
| `hvantk/tests/test_expression_builders_anndata.py::TestMkmatrixCli::test_ucsc_produces_h5ad` | CLI wiring for `mkmatrix ucsc`. |

**In-place extension**: add a `backed` parameter to `TestBuildUcscAd` via
`pytest.mark.parametrize`. Same assertions, now run in both
`backed=False` and `backed=True` modes. No new test function, no new file.

**Manual validation gate** (documented in the implementation plan; no
auto-coverage; must pass before PR):

| # | Check | Pass criterion |
|---|---|---|
| 1 | `mkmatrix ucsc --backed` on `adult-ctx-meta-atlas` | Finishes; peak RSS ≤ 3 GB; `adata.shape == (520014, 122697)`; `atlas.h5ad` reopens cleanly. |
| 2 | `expression summarize-ucsc --group-by Class` on same input | Finishes in ≤ 45 min; `summary.h5ad.obs["n_cells"].sum()` equals the cell count of metadata after filter. |
| 3 | Numerical agreement: fused `mean` layer vs. two-step `mean` layer on a toy fixture | `max_abs_diff < 1e-5`. |
| 4 | CLI smoke: `summarize-ucsc --group-by NonExistent` | Exits 2 with a `click.BadParameter` message listing the available obs columns. |

## Deferred (v2, not this PR)

- `expression summarize` gains a `--backed` mode that streams `atlas.h5ad`
  from disk. Required before anyone can aggregate from a 40 GB atlas without
  OOM; for now, users hitting that should use `summarize-ucsc` directly.
- `describe` / `markers` retrofits for backed atlases.
- Apply the same streaming pattern to Expression Atlas and other
  wide-TSV expression sources.
- Promote additional aggregations (`var`, `std`, `log2_mean`) to first-class
  output layers.

## Open implementation details

These are decided during implementation, not by this design:

1. `column_batch` default for the backed builder — likely 32 or 64; tune
   against real-run wall-clock and `atlas.h5ad` file size.
2. Whether `_iter_ucsc_rows` should take a pre-opened file handle (so
   callers can peek the header first without re-opening). Cosmetic.
3. How to surface progress logging in both new paths so the user does not
   have to pass `-v` (tracked as issue #92; likely a click-level `err=True`
   heartbeat every N rows). Deferred to #92's own PR.

## Verified against the ecosystem during design

- `anndata.experimental` / `anndata.io`:
  `CSRDataset` and `CSCDataset` both expose `.append(sparse_matrix)` for
  incremental on-disk writes (source:
  `anndata/src/anndata/_core/sparse_dataset.py::BaseCompressedSparseDataset.append`).
  CSR append requires matching column count; CSC append requires matching
  row count. Parse order in UCSC is gene-by-gene → CSC matches naturally.
- `indptr_dtype=np.int64` option on `CSCDataset` prevents `OverflowError`
  for `> 2^31` nonzeros (required at our scale).
- `scanpy.get.aggregate` composite group labels are `"_"`-joined across the
  `by=[...]` columns; the fused path replicates that convention so the
  emitted `.h5ad` shares `obs_names` format with the existing two-step
  output.

## Deliverables (this PR)

1. Extract `_iter_ucsc_rows` from `create_anndata_from_ucsc_matrix`.
2. New `build_ucsc_atlas_backed` with CSC incremental writes.
3. New `summarize_ucsc_streaming` with `np.bincount` aggregation.
4. `build_ucsc_ad` gains `backed: bool | None` (auto-select default).
5. `mkmatrix ucsc` gains `--backed / --no-backed`.
6. New `expression summarize-ucsc` CLI subcommand.
7. Parametrized extension of `TestBuildUcscAd` covering `backed=True`.
8. Manual validation notes recorded against `adult-ctx-meta-atlas`.
