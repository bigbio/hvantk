# Single-cell aggregation CLI — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship a scverse-native `hvantk expression summarize` that emits a canonical `.h5ad` (groups × genes, aggregations as layers) via `scanpy.get.aggregate`, plus a stream-capable `mkmatrix ucsc` that handles multi-GB gzipped expression matrices without OOM.

**Architecture:** Two-step Builder/Streamer pipeline. `build_ucsc_ad` streams raw UCSC gz → CSR-sparse AnnData (`.h5ad`); `summarize_expression_ad` wraps `sc.get.aggregate` to produce a small aggregated `.h5ad`. Notebook K.1 consumes the aggregated artifact directly.

**Tech Stack:** Python 3.9+, scanpy ≥ 1.10 (`sc.get.aggregate`), anndata ≥ 0.10, scipy sparse, pandas, Click CLI.

**Spec:** `docs_site/specs/2026-04-16-single-cell-aggregation-cli-design.md`

**Git workflow:** Project policy forbids direct commits to `dev` or `main`. Work on a feature branch (e.g., `feat/sc-aggregation-cli`). Open a PR into `dev` at the end. Consider using `superpowers:using-git-worktrees` for isolation.

---

## File map

| File | Action | Responsibility |
|---|---|---|
| `hvantk/tables/ucsc.py` | Modify | Stream chunks from gz → CSR → AnnData. |
| `hvantk/tables/matrix_builders.py` | Modify | Plumb `chunk_size` kwarg through `build_ucsc_ad`. |
| `hvantk/commands/make_matrix_cli.py` | Modify | Expose `--chunk-size` on `mkmatrix ucsc`. |
| `hvantk/utils/matrix_utils.py` | Modify | Rewrite `summarize_expression_ad` on `sc.get.aggregate`. |
| `hvantk/commands/summarize_expression_cli.py` | Modify | `.h5ad` output, `.h5ad` suffix validation, `--min-cells` default 10, new report card. |
| `hvantk/tests/test_matrix_utils_anndata.py` | Modify | Migrate `TestSummarizeExpressionAd` to new AnnData return contract. |
| `local/notebooks/ptm-eda/notebook_k1_brain_metaatlas_ptm.ipynb` | Modify | Replace Cell 3 streaming loop with CLI call + AnnData load. |

No new files. No new tests (`MEMORY.md → feedback_no_more_tests.md`).

---

## Task 1: Stream-capable `create_anndata_from_ucsc_matrix`

**Files:**
- Modify: `hvantk/tables/ucsc.py` (function `create_anndata_from_ucsc_matrix`, lines 69–142)
- Test: `hvantk/tests/test_expression_builders_anndata.py::TestBuildUcscAd` (existing — must still pass, unchanged)

Contract stays the same: input (gz or plain TSV) → `AnnData(cells × genes)` sparse CSR. Only the body changes. Because contract is unchanged, existing tests are the regression harness.

- [ ] **Step 1: Verify existing ucsc-builder tests pass against current code**

Run: `pytest hvantk/tests/test_expression_builders_anndata.py::TestBuildUcscAd -v`

Expected: 4 tests PASS.

- [ ] **Step 2: Replace the body of `create_anndata_from_ucsc_matrix` with a streamed builder**

In `hvantk/tables/ucsc.py`, replace the function body starting at line 106 (the `if not os.path.exists(...)` block) with the streamed implementation. Add a `chunk_size` kwarg to the signature. The full replacement function:

```python
def create_anndata_from_ucsc_matrix(
    expression_matrix_path: str,
    metadata_df: pd.DataFrame = None,
    gene_column: str = UCSC_GENE_COLUMN,
    delimiter: str = "\t",
    split_gene_field: bool = True,
    chunk_size: int = 500,
) -> ad.AnnData:
    """Create an AnnData object from a UCSC Cell Browser expression matrix.

    The input file has genes as rows and cells as columns. This function
    streams the file in chunks of ``chunk_size`` genes, converts each chunk
    to a ``scipy.sparse.csr_matrix``, and vstacks them, so peak RAM stays
    roughly ``chunk_size × n_cells × 4 B`` plus the growing sparse result.
    The stacked matrix is transposed to the AnnData convention
    (obs = cells, var = genes) and stored as ``float32`` CSR.

    Parameters
    ----------
    expression_matrix_path : str
        Path to the expression TSV or `.tsv.gz` (genes x cells).
    metadata_df : pd.DataFrame, optional
        Cell metadata to join into ``adata.obs``. Index must match cell ids.
    gene_column : str
        Name of the first column containing gene identifiers.
    delimiter : str
        Column delimiter (default tab).
    split_gene_field : bool
        If True, split gene names on ``|`` and keep only the first element.
    chunk_size : int
        Gene rows per streaming chunk. Lower ⇒ less peak RAM, slower.

    Returns
    -------
    ad.AnnData
        AnnData with shape (n_cells, n_genes), sparse CSR X.

    Raises
    ------
    FileNotFoundError
        If *expression_matrix_path* does not exist.
    ValueError
        If a chunk fails float32 conversion (names the offending gene row).
    """
    if not os.path.exists(expression_matrix_path):
        raise FileNotFoundError(
            f"Expression matrix file not found: {expression_matrix_path}"
        )

    logger.info(
        "Streaming expression matrix from %s (chunk_size=%d)",
        expression_matrix_path,
        chunk_size,
    )

    gene_name_chunks = []
    sparse_chunks = []
    cell_ids = None
    n_seen = 0

    reader = pd.read_csv(
        expression_matrix_path,
        sep=delimiter,
        header=0,
        index_col=0,
        chunksize=chunk_size,
    )
    for chunk in reader:
        if cell_ids is None:
            cell_ids = list(chunk.columns)

        chunk_genes = chunk.index.astype(str)
        if split_gene_field:
            chunk_genes = chunk_genes.str.split("|").str[0]
        gene_name_chunks.append(np.asarray(chunk_genes))

        try:
            chunk_values = chunk.to_numpy(dtype=np.float32, copy=False)
        except (ValueError, TypeError) as exc:
            bad_gene = chunk_genes[0] if len(chunk_genes) else "<unknown>"
            raise ValueError(
                f"Non-numeric value in chunk starting at gene {bad_gene!r}: {exc}"
            ) from exc

        sparse_chunks.append(sparse.csr_matrix(chunk_values))
        n_seen += len(chunk_genes)
        if n_seen % (chunk_size * 10) == 0:
            logger.info("  streamed %d genes", n_seen)

    if not sparse_chunks:
        raise ValueError(
            f"Expression matrix {expression_matrix_path!r} contained no data rows."
        )

    # Genes-x-cells sparse → cells-x-genes AnnData convention.
    X = sparse.vstack(sparse_chunks, format="csr").T.tocsr()
    gene_names = np.concatenate(gene_name_chunks)

    var = pd.DataFrame(index=pd.Index(gene_names, name=gene_column))
    obs = pd.DataFrame(index=pd.Index(cell_ids, name=UCSC_CELL_ID_COLUMN))

    adata = ad.AnnData(X=X, obs=obs, var=var)

    if metadata_df is not None:
        common = adata.obs.index.intersection(metadata_df.index)
        if len(common) == 0:
            logger.warning(
                "No overlapping cell ids between expression matrix and metadata"
            )
        meta_aligned = metadata_df.reindex(adata.obs.index)
        for col in meta_aligned.columns:
            adata.obs[col] = meta_aligned[col].values

    logger.info("Built AnnData: %d cells × %d genes", adata.shape[0], adata.shape[1])
    return adata
```

- [ ] **Step 3: Run existing ucsc-builder tests to verify no regression**

Run: `pytest hvantk/tests/test_expression_builders_anndata.py::TestBuildUcscAd -v`

Expected: all 4 tests PASS (no test changes; same contract). If any fail, the streaming body diverged from the previous contract — fix before continuing.

- [ ] **Step 4: Commit**

```bash
git add hvantk/tables/ucsc.py
git commit -m "$(cat <<'EOF'
feat(ucsc): stream exprMatrix in chunks to avoid dense OOM on large atlases

create_anndata_from_ucsc_matrix now streams the gz/tsv in chunks of
`chunk_size` genes, converting each chunk to a CSR sparse block before
vstacking. Peak RAM ≈ chunk_size × n_cells × 4 B + growing sparse result,
making 500k-cell UCSC atlases feasible on a laptop.
EOF
)"
```

---

## Task 2: Plumb `chunk_size` through `build_ucsc_ad`

**Files:**
- Modify: `hvantk/tables/matrix_builders.py` (function `build_ucsc_ad`, lines 24–84)
- Test: `hvantk/tests/test_expression_builders_anndata.py::TestBuildUcscAd` (unchanged — regression harness).

- [ ] **Step 1: Add `chunk_size` kwarg to `build_ucsc_ad` and pass it through**

In `hvantk/tables/matrix_builders.py`, update the `build_ucsc_ad` signature and its call to `create_anndata_from_ucsc_matrix`:

```python
def build_ucsc_ad(
    expression_matrix_path: str,
    metadata_path: str,
    output_path: Optional[str] = None,
    gene_column: str = "gene",
    delimiter: str = "\t",
    split_gene_field: bool = True,
    overwrite: bool = False,
    chunk_size: int = 500,
) -> "ad.AnnData":
    """Build an AnnData object from UCSC Cell Browser expression + metadata.

    Parameters
    ----------
    expression_matrix_path : str
        Path to expression TSV (genes x cells).
    metadata_path : str
        Path to metadata TSV.
    output_path : str, optional
        If provided, save the AnnData as ``.h5ad``.
    gene_column : str
        Name of the gene identifier column (default ``"gene"``).
    delimiter : str
        Column delimiter (default tab).
    split_gene_field : bool
        Split pipe-separated gene names, keeping the first element.
    overwrite : bool
        Allow overwriting *output_path* if it exists.
    chunk_size : int
        Gene rows per streaming chunk (default 500).

    Returns
    -------
    ad.AnnData
        Expression AnnData with metadata in ``obs`` and provenance in ``uns``.
    """
    from hvantk.tables.ucsc import load_ucsc_metadata, create_anndata_from_ucsc_matrix
    from hvantk.core.anndata_utils import (
        build_anndata_metadata,
        annotate_column_summary_ad,
        save_anndata,
    )

    logger.info("Loading UCSC metadata from %s", metadata_path)
    metadata_df = load_ucsc_metadata(metadata_path, sep=delimiter)

    logger.info("Creating AnnData from UCSC expression matrix")
    adata = create_anndata_from_ucsc_matrix(
        expression_matrix_path=expression_matrix_path,
        metadata_df=metadata_df,
        gene_column=gene_column,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
        chunk_size=chunk_size,
    )

    adata.uns["hvantk_metadata"] = build_anndata_metadata(
        "UCSC", expression_matrix_path
    )
    annotate_column_summary_ad(adata)

    if output_path:
        save_anndata(adata, output_path, overwrite=overwrite)

    return adata
```

- [ ] **Step 2: Verify existing tests still pass**

Run: `pytest hvantk/tests/test_expression_builders_anndata.py::TestBuildUcscAd -v`

Expected: 4 tests PASS.

- [ ] **Step 3: Commit**

```bash
git add hvantk/tables/matrix_builders.py
git commit -m "feat(ucsc): expose chunk_size on build_ucsc_ad for streaming tune"
```

---

## Task 3: Expose `--chunk-size` on `mkmatrix ucsc`

**Files:**
- Modify: `hvantk/commands/make_matrix_cli.py` (function `mkmatrix_ucsc`, lines 49–93)
- Test: `hvantk/tests/test_expression_builders_anndata.py::TestMkmatrixCli::test_ucsc_produces_h5ad` (unchanged).

- [ ] **Step 1: Add `--chunk-size` option to the `ucsc` subcommand**

In `hvantk/commands/make_matrix_cli.py`, add one click option and one kwarg pass-through. The final `mkmatrix_ucsc` command:

```python
@mkmatrix_group.command("ucsc")
@click.option(
    "-e",
    "--expression-matrix",
    "expression_matrix",
    required=True,
    type=click.Path(exists=True),
)
@click.option("-m", "--metadata", required=True, type=click.Path(exists=True))
@click.option(
    "-o",
    "--output",
    "output_path",
    required=True,
    type=click.Path(),
    help="Output path for AnnData (.h5ad).",
)
@click.option("-g", "--gene-column", default="gene", show_default=True)
@click.option("-d", "--delimiter", default="\t", show_default=True)
@click.option(
    "--split-gene-field/--no-split-gene-field", default=True, show_default=True
)
@click.option(
    "--chunk-size",
    type=int,
    default=500,
    show_default=True,
    help="Gene rows per streaming chunk. Lower ⇒ less peak RAM, slower.",
)
@click.option("-w", "--overwrite", is_flag=True)
def mkmatrix_ucsc(
    expression_matrix,
    metadata,
    output_path,
    gene_column,
    delimiter,
    split_gene_field,
    chunk_size,
    overwrite,
):
    """Build an AnnData object from UCSC Cell Browser expression + metadata files."""
    logger.info("Building UCSC AnnData")
    adata = _build_ucsc_ad(
        expression_matrix_path=expression_matrix,
        metadata_path=metadata,
        output_path=output_path,
        gene_column=gene_column,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
        chunk_size=chunk_size,
        overwrite=overwrite,
    )
    click.echo(f"AnnData created at {output_path}")
    click.echo(f"  Shape: {adata.shape[0]} obs x {adata.shape[1]} vars")
```

- [ ] **Step 2: Verify CLI test still passes**

Run: `pytest hvantk/tests/test_expression_builders_anndata.py::TestMkmatrixCli::test_ucsc_produces_h5ad -v`

Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add hvantk/commands/make_matrix_cli.py
git commit -m "feat(mkmatrix): add --chunk-size to mkmatrix ucsc"
```

---

## Task 4: Migrate `TestSummarizeExpressionAd` to new AnnData contract (fail-first)

**Files:**
- Modify: `hvantk/tests/test_matrix_utils_anndata.py` (class `TestSummarizeExpressionAd`, lines 66–78)

This is the TDD step: update the existing assertions to the new contract, verify they fail against the current DataFrame-returning implementation, then flip the implementation in Task 5 to make them pass.

- [ ] **Step 1: Replace `TestSummarizeExpressionAd` with the new contract**

Replace the class (lines 66–78) in `hvantk/tests/test_matrix_utils_anndata.py` with:

```python
class TestSummarizeExpressionAd:
    def test_returns_anndata(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        assert isinstance(result, ad.AnnData)

    def test_shape_is_groups_by_genes(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        expected_groups = len(set(test_adata.obs["cell_type"]))
        assert result.shape == (expected_groups, test_adata.n_vars)

    def test_has_required_layers(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        expected_layers = {"mean", "sum", "count_nonzero", "fraction_expressed"}
        assert expected_layers.issubset(set(result.layers.keys()))

    def test_obs_has_n_cells(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        assert "n_cells" in result.obs.columns
        assert (result.obs["n_cells"] > 0).all()
        assert int(result.obs["n_cells"].sum()) == test_adata.n_obs

    def test_groups_match_unique_cell_types(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        assert set(result.obs_names) == {"neuron", "astrocyte", "microglia"}

    def test_fraction_expressed_between_0_and_1(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        frac = np.asarray(result.layers["fraction_expressed"])
        assert frac.min() >= 0.0
        assert frac.max() <= 1.0

    def test_min_cells_filters_small_groups(self, test_adata):
        # all three cell types should be well above n=5 in the fixture (100 cells / 3),
        # but at min_cells=10_000 we should drop every group.
        result = summarize_expression_ad(
            test_adata, group_by="cell_type", min_cells_per_group=10_000
        )
        assert result.n_obs == 0
```

- [ ] **Step 2: Run the updated tests — they must fail against the current implementation**

Run: `pytest hvantk/tests/test_matrix_utils_anndata.py::TestSummarizeExpressionAd -v`

Expected: every test FAILS (current implementation returns a `pd.DataFrame`, so `isinstance(result, ad.AnnData)` trips). If any test unexpectedly passes, inspect it — the contract didn't actually change there.

- [ ] **Step 3: Commit**

```bash
git add hvantk/tests/test_matrix_utils_anndata.py
git commit -m "test(summarize): migrate TestSummarizeExpressionAd to AnnData return contract"
```

---

## Task 5: Rewrite `summarize_expression_ad` on `sc.get.aggregate`

**Files:**
- Modify: `hvantk/utils/matrix_utils.py` (function `summarize_expression_ad`, lines 94–162)
- Test: `hvantk/tests/test_matrix_utils_anndata.py::TestSummarizeExpressionAd` (from Task 4).

- [ ] **Step 1: Replace the function body**

In `hvantk/utils/matrix_utils.py`, replace the `summarize_expression_ad` function (and its imports block at the top if needed) with:

```python
def summarize_expression_ad(
    adata: ad.AnnData,
    group_by: Union[str, List[str]],
    filter_by: Optional[Dict[str, Union[str, List[str]]]] = None,
    min_cells_per_group: int = 10,
) -> ad.AnnData:
    """Collapse an AnnData expression matrix into a per-group, per-gene AnnData.

    Thin wrapper around :func:`scanpy.get.aggregate` that also derives
    ``fraction_expressed`` and attaches per-group cell counts.

    Parameters
    ----------
    adata
        Expression AnnData (obs = cells/samples, var = genes).
    group_by
        One or more obs columns to group by. Multi-column groupings produce
        ``obs_names`` joined with ``_`` (matches ``scanpy.get.aggregate``).
    filter_by
        Optional pre-filter passed to :func:`filter_by_metadata_ad`.
    min_cells_per_group
        Drop groups with fewer than this many cells.

    Returns
    -------
    ad.AnnData
        Shape ``(n_groups, n_genes)`` with layers ``mean``, ``sum``,
        ``count_nonzero``, ``fraction_expressed``. ``obs["n_cells"]`` stores
        the per-group cell count.
    """
    import scanpy as sc

    if filter_by:
        adata = filter_by_metadata_ad(adata, filter_by)

    by = [group_by] if isinstance(group_by, str) else list(group_by)

    agg = sc.get.aggregate(
        adata,
        by=by,
        func=["mean", "sum", "count_nonzero"],
    )

    # sc.get.aggregate does not stash per-group cell counts; compute from the
    # pre-aggregate adata, joining group_by columns with "_" to match
    # agg.obs_names.
    group_labels = adata.obs[by].astype(str).agg("_".join, axis=1)
    n_cells = group_labels.value_counts().reindex(agg.obs_names).astype(int)

    agg.obs["n_cells"] = n_cells.values
    agg.layers["fraction_expressed"] = (
        np.asarray(agg.layers["count_nonzero"]) / n_cells.values[:, None]
    )

    keep = agg.obs["n_cells"].values >= min_cells_per_group
    return agg[keep].copy()
```

Also update the imports at the top of the file (line 11) from:

```python
from typing import Any, Dict, List, Optional, Union
```

No typing changes needed — already imports what's used. Confirm `import anndata as ad` is present (line 13); it is.

- [ ] **Step 2: Run the migrated tests — they must pass**

Run: `pytest hvantk/tests/test_matrix_utils_anndata.py::TestSummarizeExpressionAd -v`

Expected: all 7 tests PASS.

- [ ] **Step 3: Run the full test_matrix_utils_anndata.py to check no collateral damage**

Run: `pytest hvantk/tests/test_matrix_utils_anndata.py -v`

Expected: all tests PASS (Describe + Filter tests are unaffected).

- [ ] **Step 4: Commit**

```bash
git add hvantk/utils/matrix_utils.py
git commit -m "$(cat <<'EOF'
feat(summarize): rewrite summarize_expression_ad on scanpy.get.aggregate

Replace the hand-rolled per-group loop with a thin wrapper around
sc.get.aggregate that returns an AnnData keyed by group (obs_names) with
mean / sum / count_nonzero / fraction_expressed in layers and per-group
cell counts in obs["n_cells"]. Multi-column group_by joins column values
with "_" to match scanpy's obs_names convention.

Breaking change: return type is now ad.AnnData (was pd.DataFrame).
EOF
)"
```

---

## Task 6: Update `summarize_expression_cmd` CLI for `.h5ad` output

**Files:**
- Modify: `hvantk/commands/summarize_expression_cli.py` (function `summarize_expression_cmd`, lines 70–199)

Deltas: enforce `.h5ad` suffix, lower `--min-cells` default from 50 to 10, replace `to_parquet` with `write_h5ad`, update the report card to reflect AnnData output.

- [ ] **Step 1: Replace the `summarize_expression_cmd` body**

In `hvantk/commands/summarize_expression_cli.py`, replace the function (starting at `@expression_group.command("summarize")` and running through the bottom of the function) with:

```python
@expression_group.command("summarize")
@click.option(
    "-m",
    "--matrix-table",
    "matrix_path",
    type=click.Path(exists=True),
    required=True,
    help="Path to an expression AnnData file (.h5ad).",
)
@click.option(
    "--group-by",
    multiple=True,
    required=True,
    help="Observation metadata field(s) to group by. Repeat for multi-field grouping.",
)
@click.option(
    "--filter-by",
    multiple=True,
    default=None,
    help="Pre-filter observations: FIELD=VALUE (repeatable). "
    "Example: --filter-by time_point=9wpc --filter-by region=LV",
)
@click.option(
    "--min-cells",
    type=int,
    default=10,
    show_default=True,
    help="Drop groups with fewer cells than this threshold.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    required=True,
    help="Output path for the aggregated AnnData (.h5ad).",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite existing output.",
)
def summarize_expression_cmd(
    matrix_path,
    group_by,
    filter_by,
    min_cells,
    output,
    overwrite,
):
    """Aggregate an expression AnnData into a per-group × per-gene AnnData.

    The output .h5ad has shape (n_groups, n_genes) with layers:
    ``mean``, ``sum``, ``count_nonzero``, ``fraction_expressed``. Per-group
    cell counts are in ``obs['n_cells']``.

    \b
    Examples:

      # Single-field grouping
      hvantk expression summarize \\
          -m data/heart_sc.h5ad \\
          --group-by cell_type \\
          -o data/heart_celltype_summary.h5ad

      # Multi-field grouping
      hvantk expression summarize \\
          -m data/heart_sc.h5ad \\
          --group-by cell_type --group-by region \\
          -o data/heart_celltype_region_summary.h5ad

      # With pre-filtering
      hvantk expression summarize \\
          -m data/heart_sc.h5ad \\
          --group-by cell_type \\
          --filter-by time_point=9wpc --filter-by region=LV \\
          --min-cells 10 \\
          -o data/heart_9wpc_LV_summary.h5ad
    """
    from pathlib import Path

    from hvantk.core.anndata_utils import load_anndata
    from hvantk.utils.matrix_utils import summarize_expression_ad

    output_path = Path(output)
    if output_path.suffix != ".h5ad":
        raise click.BadParameter(
            f"Output must end in '.h5ad' (got {output_path.suffix!r}).",
            param_hint="--output",
        )

    if output_path.exists() and not overwrite:
        click.echo(
            f"Error: Output file already exists: {output_path}\n"
            "Use --overwrite to replace it.",
            err=True,
        )
        raise SystemExit(1)

    filters = None
    if filter_by:
        filters = {}
        for item in filter_by:
            if "=" not in item:
                raise click.BadParameter(
                    f"Expected FIELD=VALUE format, got: '{item}'",
                    param_hint="--filter-by",
                )
            key, value = item.split("=", 1)
            filters[key.strip()] = value.strip()

    adata = load_anndata(matrix_path)

    missing = [col for col in group_by if col not in adata.obs.columns]
    if missing:
        available = sorted(adata.obs.columns)
        raise click.BadParameter(
            f"--group-by column(s) not in obs: {missing}. Available: {available}",
            param_hint="--group-by",
        )

    summary = summarize_expression_ad(
        adata,
        group_by=list(group_by),
        filter_by=filters,
        min_cells_per_group=min_cells,
    )

    if summary.n_obs == 0:
        click.echo(
            "Error: no groups passed --min-cells threshold. "
            "Lower --min-cells or check --filter-by / --group-by.",
            err=True,
        )
        raise SystemExit(1)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary.write_h5ad(str(output_path))

    n_cells = summary.obs["n_cells"].astype(int)
    click.echo(f"\nSummary AnnData written to: {output_path}")
    click.echo(f"  Shape:  {summary.n_obs} groups × {summary.n_vars} genes")
    click.echo(f"  Layers: {sorted(summary.layers.keys())}")
    click.echo(
        f"  Cells per group: min={int(n_cells.min())}, "
        f"max={int(n_cells.max())}, total={int(n_cells.sum())}"
    )
    click.echo("")
    click.echo(f"Group labels ({summary.n_obs}):")
    for label, n in list(zip(summary.obs_names, n_cells))[:10]:
        click.echo(f"  {label:<40s}  {int(n):>8,} cells")
    if summary.n_obs > 10:
        click.echo(f"  ... and {summary.n_obs - 10} more")
```

- [ ] **Step 2: Smoke-test the CLI end-to-end with a tiny fixture**

Run:

```bash
python - <<'PY'
import anndata as ad, numpy as np, pandas as pd, tempfile, subprocess, pathlib
rng = np.random.default_rng(0)
X = rng.poisson(0.5, size=(30, 10)).astype(np.float32)
obs = pd.DataFrame({"cell_type": ["A"]*12 + ["B"]*10 + ["C"]*8},
                    index=[f"c{i}" for i in range(30)])
var = pd.DataFrame(index=[f"g{i}" for i in range(10)])
adata = ad.AnnData(X=X, obs=obs, var=var)
with tempfile.TemporaryDirectory() as d:
    in_path  = pathlib.Path(d) / "in.h5ad"
    out_path = pathlib.Path(d) / "out.h5ad"
    adata.write_h5ad(in_path)
    r = subprocess.run(
        ["hvantk", "expression", "summarize",
         "-m", str(in_path),
         "--group-by", "cell_type",
         "--min-cells", "5",
         "-o", str(out_path)],
        capture_output=True, text=True,
    )
    print("stdout:", r.stdout)
    print("stderr:", r.stderr)
    print("exit:", r.returncode)
    assert r.returncode == 0
    loaded = ad.read_h5ad(out_path)
    print("loaded shape:", loaded.shape)
    print("layers:", sorted(loaded.layers.keys()))
    print("obs['n_cells']:", loaded.obs["n_cells"].tolist())
    assert set(loaded.layers.keys()) >= {"mean", "sum", "count_nonzero", "fraction_expressed"}
    assert loaded.shape == (3, 10)
print("OK")
PY
```

Expected: "OK" printed. Shape `(3, 10)`, all 4 layers present, `n_cells` matches the fixture.

- [ ] **Step 3: Run the fast test suite to confirm no regression**

Run: `pytest hvantk/tests/test_matrix_utils_anndata.py hvantk/tests/test_expression_builders_anndata.py -v`

Expected: all tests PASS.

- [ ] **Step 4: Commit**

```bash
git add hvantk/commands/summarize_expression_cli.py
git commit -m "$(cat <<'EOF'
feat(expression): summarize writes .h5ad (groups × genes with layers)

- Enforce --output suffix ==  .h5ad
- Default --min-cells lowered 50 → 10 (matches planning doc and typical sc atlases)
- Replace long-form parquet write with AnnData.write_h5ad
- Report card now surfaces n_groups, n_genes, layers, and per-group cell counts
- Bail early with a helpful message when a --group-by column is missing or
  every group drops below --min-cells
EOF
)"
```

---

## Task 7: Retrofit notebook K.1 Cell 3

**Files:**
- Modify: `local/notebooks/ptm-eda/notebook_k1_brain_metaatlas_ptm.ipynb` (Cell 7 in the ipynb JSON, labeled "Cell 3" in the markdown; currently the ~45-line streaming loop)

The atlas `.h5ad` must be produced once via `hvantk mkmatrix ucsc` before running the notebook (see Task 8 for the command). Cell 3 becomes a CLI call + AnnData load.

- [ ] **Step 1: Replace the source of the streaming cell via a reproducible script**

Run this script — it locates the streaming cell by its leading marker `# === Cell 3 ===` + `# Vectorized chunked build`, replaces its `source`, clears outputs, and rewrites the notebook:

```bash
python - <<'PY'
import json, pathlib

NB = pathlib.Path("local/notebooks/ptm-eda/notebook_k1_brain_metaatlas_ptm.ipynb")
MARKER = "# === Cell 3 ==="
SENTINEL = "Vectorized chunked build"

new_source = '''# === Cell 3 ===
# Aggregation is now done by the hvantk CLI, which leans on
# scanpy.get.aggregate over the atlas .h5ad built from mkmatrix ucsc.
# One-time preparation (run outside the notebook):
#   hvantk mkmatrix ucsc -e exprMatrix.tsv.gz -m meta.tsv -o atlas.h5ad
import subprocess
import anndata as ad

BRAIN_ATLAS_H5AD  = f"{META_ATLAS_DIR}/adult_ctx_meta_atlas.h5ad"
BRAIN_CACHE_H5AD  = f"{OUTPUT_DIR}/brain_metaatlas_gene_class_mean_expr.h5ad"

if not os.path.exists(BRAIN_CACHE_H5AD):
    log.info("Aggregating atlas per %s via hvantk expression summarize...", CELLTYPE_COL)
    subprocess.run(
        [
            "hvantk", "expression", "summarize",
            "-m", BRAIN_ATLAS_H5AD,
            "--group-by", CELLTYPE_COL,
            "--min-cells", "10",
            "-o", BRAIN_CACHE_H5AD,
            "--overwrite",
        ],
        check=True,
    )
else:
    log.info("Loading cached aggregated AnnData from %s", BRAIN_CACHE_H5AD)

brain_agg = ad.read_h5ad(BRAIN_CACHE_H5AD)
brain_wide = (
    pd.DataFrame(
        brain_agg.layers["mean"].T,
        index=brain_agg.var_names,
        columns=brain_agg.obs_names,
    )
    .rename_axis("gene_symbol")
)
# De-dup on gene symbol (meta-atlas occasionally has duplicate var index)
brain_wide = brain_wide[~brain_wide.index.duplicated(keep="first")]

print(f"Meta-atlas wide matrix: {brain_wide.shape[0]:,} genes \u00d7 {brain_wide.shape[1]} classes")
print(f"Classes: {list(brain_wide.columns)}")
print(brain_wide.describe().T)
'''

nb = json.loads(NB.read_text())

matches = []
for i, cell in enumerate(nb["cells"]):
    if cell.get("cell_type") != "code":
        continue
    src = "".join(cell.get("source", []))
    if src.startswith(MARKER) and SENTINEL in src:
        matches.append(i)

if len(matches) != 1:
    raise SystemExit(
        f"Expected exactly one cell starting with {MARKER!r} and containing "
        f"{SENTINEL!r}; found {len(matches)}. Inspect the notebook manually."
    )

idx = matches[0]
cell = nb["cells"][idx]
cell["source"] = [line + "\n" for line in new_source.splitlines()]
if cell["source"]:
    cell["source"][-1] = cell["source"][-1].rstrip("\n")
cell["outputs"] = []
cell["execution_count"] = None

NB.write_text(json.dumps(nb, indent=1) + "\n")
print(f"Replaced cell #{idx} in {NB}")
PY
```

Expected: prints `Replaced cell #N in local/notebooks/...` (exactly one match).

- [ ] **Step 2: Verify the notebook still parses**

Run:

```bash
python -c "import json; json.load(open('local/notebooks/ptm-eda/notebook_k1_brain_metaatlas_ptm.ipynb'))"
```

Expected: no output (valid JSON).

- [ ] **Step 3: Verify nbformat validates**

Run:

```bash
python -c "import nbformat; nbformat.read('local/notebooks/ptm-eda/notebook_k1_brain_metaatlas_ptm.ipynb', as_version=4)"
```

Expected: no error.

- [ ] **Step 4: Commit**

```bash
git add local/notebooks/ptm-eda/notebook_k1_brain_metaatlas_ptm.ipynb
git commit -m "$(cat <<'EOF'
refactor(notebook-k1): consume pre-aggregated AnnData from hvantk CLI

Cell 3 previously re-implemented a vectorized chunked streaming loop
against exprMatrix.tsv.gz. Replace with a subprocess call to
`hvantk expression summarize` and an AnnData load of the resulting .h5ad.
The wide pandas DataFrame (gene_symbol × Class) that the rest of the
notebook consumes is rederived from layers["mean"] in two lines.

Prereq: one-time build of atlas.h5ad via `hvantk mkmatrix ucsc`.
EOF
)"
```

---

## Task 8: Manual validation on the brain meta-atlas

No automated tests for this step. Purpose: verify end-to-end behavior on realistic data before opening the PR. If any check fails, pause and investigate before proceeding.

- [ ] **Step 1: Build the canonical atlas `.h5ad` via the streamed `mkmatrix ucsc`**

```bash
hvantk mkmatrix ucsc \
    -e /Users/enrique/projects/github/pyvatk/local/data/UCSC/adult-ctx-meta-atlas/adult-ctx-meta-atlas/exprMatrix.tsv.gz \
    -m /Users/enrique/projects/github/pyvatk/local/data/UCSC/adult-ctx-meta-atlas/adult-ctx-meta-atlas/meta.tsv \
    -o /Users/enrique/projects/github/pyvatk/local/data/UCSC/adult-ctx-meta-atlas/adult-ctx-meta-atlas/adult_ctx_meta_atlas.h5ad \
    --chunk-size 500
```

Expected: finishes without OOM. Log lines show streamed gene counts increasing. Final shape roughly ~520k cells × ~20k genes.

Acceptance: wall-clock under ~45 min on a laptop; no RSS spike beyond ~6 GB.

- [ ] **Step 2: Run `expression summarize` over the atlas**

```bash
hvantk expression summarize \
    -m /Users/enrique/projects/github/pyvatk/local/data/UCSC/adult-ctx-meta-atlas/adult-ctx-meta-atlas/adult_ctx_meta_atlas.h5ad \
    --group-by Class \
    --min-cells 10 \
    -o /Users/enrique/projects/github/pyvatk/output/ptm_assembly/brain_metaatlas_gene_class_mean_expr.h5ad \
    --overwrite
```

Expected: completes in 10–15 min. Report card shows ~10 groups (Jorstad Class labels) and ~20k genes with all four layers.

- [ ] **Step 3: Spot-check numerical agreement against the old pickle cache**

If the old pickle still exists at `output/ptm_assembly/brain_metaatlas_gene_class_mean_expr.pkl`, run:

```bash
python - <<'PY'
import pandas as pd, anndata as ad, numpy as np

old = pd.read_pickle("output/ptm_assembly/brain_metaatlas_gene_class_mean_expr.pkl")
new = ad.read_h5ad("output/ptm_assembly/brain_metaatlas_gene_class_mean_expr.h5ad")
new_wide = pd.DataFrame(
    new.layers["mean"].T, index=new.var_names, columns=new.obs_names
).rename_axis("gene_symbol")
new_wide = new_wide[~new_wide.index.duplicated(keep="first")]

common_genes = old.index.intersection(new_wide.index)
common_cls = old.columns.intersection(new_wide.columns)
o = old.loc[common_genes, common_cls].to_numpy(dtype=np.float64)
n = new_wide.loc[common_genes, common_cls].to_numpy(dtype=np.float64)
abs_err = np.nanmax(np.abs(o - n))
rel_err = np.nanmax(np.abs(o - n) / np.where(np.abs(o) > 1e-8, np.abs(o), 1.0))
print(f"overlap: {len(common_genes)} genes × {len(common_cls)} classes")
print(f"max abs diff: {abs_err:.3e}")
print(f"max rel diff: {rel_err:.3e}")
assert abs_err < 1e-3 or rel_err < 1e-3, "mean expression diverged — investigate"
print("OK — new pipeline agrees with old pickle within tolerance")
PY
```

Expected: "OK — new pipeline agrees with old pickle within tolerance".

- [ ] **Step 4: Run notebook K.1 end-to-end (or at least Cells 1–3)**

Either execute via `jupyter nbconvert --to notebook --execute ...`, or open in JupyterLab and run Cells 1 through 3. Verify `brain_wide.shape` matches the previous run and `brain_wide.describe()` looks sane.

- [ ] **Step 5: Clean up — delete the now-stale pickle if validation passed**

```bash
rm -f output/ptm_assembly/brain_metaatlas_gene_class_mean_expr.pkl
```

- [ ] **Step 6: Open the PR**

```bash
git push -u origin feat/sc-aggregation-cli
gh pr create --base dev --title "feat(expression): scverse-native summarize + streaming mkmatrix ucsc" --body "$(cat <<'EOF'
## Summary
- Streams the UCSC exprMatrix gz in chunks → CSR → AnnData, so `mkmatrix ucsc` handles 500k-cell atlases without OOM (new `--chunk-size`).
- Rewrites `expression summarize` on top of `scanpy.get.aggregate`; output is now a canonical `.h5ad` of shape (n_groups × n_genes) with `mean`, `sum`, `count_nonzero`, `fraction_expressed` as layers and per-group cell counts in `obs["n_cells"]`.
- Drops `median` / pickle / Hail-Table output (deferred to v2; no concrete consumer yet).
- Retrofits notebook K.1 Cell 3 to consume the aggregated `.h5ad` directly instead of streaming the 5 GB gz in-notebook.

## Test plan
- [x] `pytest hvantk/tests/test_expression_builders_anndata.py::TestBuildUcscAd -v`
- [x] `pytest hvantk/tests/test_matrix_utils_anndata.py -v`
- [x] Manual end-to-end build + summarize of the adult-cortex meta-atlas (~520k cells)
- [x] Numerical agreement with the prior pickle cache (max abs diff < 1e-3)
- [x] Notebook K.1 Cells 1–3 run cleanly against the new artifact
EOF
)"
```

---

## Deferred (NOT in this PR)

- `--backed` mode on `expression summarize` for atlases > RAM.
- Exporters from the aggregated `.h5ad` to wide pickle / parquet / Hail Table.
- Retrofit of notebooks F / G / H.
- Promote `sum` / `var` to first-class `--agg` selectors.
