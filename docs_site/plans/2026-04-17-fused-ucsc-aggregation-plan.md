# Fused UCSC aggregation + backed-write builder — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship two scalable paths for 500k+-cell UCSC atlases: (a) `mkmatrix ucsc` with backed CSC writes (no RAM materialization), (b) a new `expression summarize-ucsc` that streams + aggregates in one pass and emits a `groups × genes` `.h5ad` directly.

**Architecture:** A shared row-streaming helper (`_iter_ucsc_rows`) extracted from the existing parser. The backed builder appends CSC columns to `atlas.h5ad`; the fused aggregator reduces each row via `np.bincount` into pre-allocated `(n_groups, n_genes)` accumulators. Both paths share the same parser, same validation, and identical output schema.

**Tech Stack:** Python 3.10+, anndata ≥0.10 (`anndata.io.sparse_dataset`, `CSCDataset.append`), scipy sparse, numpy, pandas, Click CLI.

**Spec:** `docs_site/specs/2026-04-17-fused-ucsc-aggregation-design.md`

**Git workflow:** Project policy forbids direct commits to `dev`/`main`. Work continues on `feat/sc-aggregation-cli`. Open a PR into `dev` at the end. Commits require user approval per project policy.

---

## File map

| File | Action | Responsibility |
|---|---|---|
| `hvantk/tables/ucsc.py` | Modify | Extract `_iter_ucsc_rows`. Keep `create_anndata_from_ucsc_matrix` as in-memory path. Add `build_ucsc_atlas_backed` (CSC backed writes) and `summarize_ucsc_streaming` (fused aggregator). |
| `hvantk/tables/matrix_builders.py` | Modify | `build_ucsc_ad` gains `backed: bool \| None` kwarg with auto-select. |
| `hvantk/commands/make_matrix_cli.py` | Modify | `mkmatrix ucsc` gains `--backed / --no-backed` flag. |
| `hvantk/commands/summarize_expression_cli.py` | Modify | New `summarize-ucsc` subcommand. |
| `hvantk/tests/test_expression_builders_anndata.py` | Modify | Parametrize `TestBuildUcscAd` with `backed=[False, True]`. |
| `local/scripts/2026-04-17-ucsc-sandbox.py` | Create | Disposable sandbox validating pyarrow-free fused `np.bincount` and CSC backed-write patterns before touching hvantk. |

No new tests files. Existing `TestBuildUcscAd` and `TestSummarizeExpressionAd` are the regression harness.

---

## Task 1: Sandbox — validate CSC backed-write + `np.bincount` patterns on real data

**Files:**
- Create: `local/scripts/2026-04-17-ucsc-sandbox.py`

Per the prototyping-workflow memory, validate both load-bearing patterns against real data before touching hvantk. The sandbox proves three things on a small slice of `adult-ctx-meta-atlas`:

1. `anndata.io.sparse_dataset(group)` + `.append(csc_block)` builds a valid `.h5ad` that reopens cleanly with the correct shape
2. `np.bincount(group_idx, weights=row)` produces per-group sums that match the reference `sc.get.aggregate` output numerically
3. Peak RSS stays bounded (no accidental materialization)

- [ ] **Step 1: Write the sandbox script**

Create `local/scripts/2026-04-17-ucsc-sandbox.py`:

```python
"""Sandbox — validate CSC backed-write and np.bincount fused aggregation
on a slice of the real adult-ctx-meta-atlas file. Disposable; not imported
by hvantk."""
from __future__ import annotations
import gzip, itertools, os, resource, time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from anndata.io import sparse_dataset, write_elem
from scipy import sparse

EXPR = Path("local/data/UCSC/adult-ctx-meta-atlas/adult-ctx-meta-atlas/exprMatrix.tsv.gz")
META = Path("local/data/UCSC/adult-ctx-meta-atlas/adult-ctx-meta-atlas/meta.tsv")
OUT_BACKED = Path("local/scripts/sandbox-backed.h5ad")
N_ROWS_PROBE = 200   # small slice — fast feedback


def parse_rows(path: Path, max_rows: int):
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").rstrip("\r").split("\t")
        cell_ids = header[1:]
        for line in itertools.islice(fh, max_rows):
            line = line.rstrip("\n").rstrip("\r")
            tab = line.find("\t")
            gene = line[:tab].split("|", 1)[0]
            row = np.fromstring(line[tab + 1:], sep="\t", dtype=np.float32)
            if row.shape[0] != len(cell_ids):
                raise ValueError(f"bad row for {gene!r}")
            yield cell_ids, gene, row


def _check_api():
    """Pre-flight — fail fast if required anndata surface is missing."""
    from anndata.io import sparse_dataset, write_elem  # noqa: F401
    # Verify we can pass indptr_dtype via dataset_kwargs (load-bearing for
    # >2^31 nonzeros). The signature may vary across versions.
    import inspect
    sig = inspect.signature(write_elem)
    if "dataset_kwargs" not in sig.parameters:
        raise RuntimeError(
            f"anndata.io.write_elem lacks dataset_kwargs; signature={sig}. "
            f"Upgrade anndata (>=0.10) or adapt the backed builder."
        )
    print(f"[api] anndata.io.write_elem signature: {sig}")


def test_backed_write():
    print("=== Backed CSC write ===")
    if OUT_BACKED.exists():
        OUT_BACKED.unlink()
    t0 = time.time()

    rows_iter = parse_rows(EXPR, N_ROWS_PROBE)
    cell_ids, gene0, row0 = next(rows_iter)
    n_cells = len(cell_ids)
    gene_names = [gene0]

    # Initialize backed CSC X with an empty (n_cells, 0) seed. Pass
    # indptr_dtype=int64 so large atlases don't overflow the default int32.
    f = h5py.File(OUT_BACKED, "w")
    seed = sparse.csc_matrix((n_cells, 0), dtype=np.float32)
    write_elem(f, "X", seed, dataset_kwargs={"indptr_dtype": "int64"})
    X = sparse_dataset(f["X"])
    # Sanity: confirm the on-disk indptr dataset is int64.
    assert f["X/indptr"].dtype == np.int64, f["X/indptr"].dtype
    print(f"[check] indptr dtype on disk: {f['X/indptr'].dtype}")

    # Append first row as a 1-column CSC block.
    X.append(sparse.csc_matrix(row0.reshape(-1, 1)))

    for _, gene, row in rows_iter:
        X.append(sparse.csc_matrix(row.reshape(-1, 1)))
        gene_names.append(gene)

    # Write obs/var placeholders.
    obs = pd.DataFrame(index=pd.Index(cell_ids, name="cell_id"))
    var = pd.DataFrame(index=pd.Index(gene_names, name="gene"))
    write_elem(f, "obs", obs)
    write_elem(f, "var", var)
    f.close()

    rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024
    print(f"wrote {N_ROWS_PROBE} cols in {time.time()-t0:.1f}s, peak_rss~{rss_mb:.0f} MB")

    # Reopen and verify.
    adata = ad.read_h5ad(OUT_BACKED)
    assert adata.shape == (n_cells, N_ROWS_PROBE), adata.shape
    assert adata.var_names[:3].tolist() == gene_names[:3]
    assert sparse.issparse(adata.X)
    print(f"reopen OK: shape={adata.shape}, X.format={adata.X.format}")


def test_fused_aggregation():
    print("=== Fused np.bincount vs sc.get.aggregate ===")
    meta = pd.read_csv(META, sep="\t", index_col=0)
    meta.index.name = "cell_id"
    meta.columns = [c.replace(".", "_") for c in meta.columns]

    # Build a small reference AnnData from the first N rows, then aggregate
    # both ways and compare.
    rows_iter = parse_rows(EXPR, N_ROWS_PROBE)
    cell_ids, gene0, row0 = next(rows_iter)
    rest = list(rows_iter)
    all_rows = [row0] + [r for _, _, r in rest]
    all_genes = [gene0] + [g for _, g, _ in rest]

    # Reference: full dense AnnData of the slice, then sc.get.aggregate.
    X_dense = np.vstack(all_rows).T  # (n_cells, n_genes)
    meta_aligned = meta.reindex(cell_ids)
    group_col = "Class" if "Class" in meta_aligned.columns else meta_aligned.columns[0]
    mask = meta_aligned[group_col].notna().values
    ref_adata = ad.AnnData(
        X=sparse.csr_matrix(X_dense[mask]),
        obs=meta_aligned.loc[mask, [group_col]].copy(),
        var=pd.DataFrame(index=all_genes),
    )
    ref_agg = sc.get.aggregate(ref_adata, by=group_col, func=["mean"])
    ref_mean = np.asarray(ref_agg.layers["mean"])

    # Fused: np.bincount on the same slice.
    labels = meta_aligned.loc[mask, group_col].astype(str).values
    codes, uniques = pd.factorize(labels, sort=True)
    n_groups = len(uniques)
    n_cells_per_group = np.bincount(codes, minlength=n_groups).astype(np.int64)
    sum_matrix = np.zeros((n_groups, len(all_genes)), dtype=np.float64)
    for j, row in enumerate(all_rows):
        sum_matrix[:, j] = np.bincount(codes, weights=row[mask], minlength=n_groups)
    fused_mean = sum_matrix / n_cells_per_group[:, None]

    # Reorder ref to match uniques order (sc.get.aggregate sorts by category).
    ref_order = [list(ref_agg.obs_names).index(u) for u in uniques]
    ref_mean_reordered = ref_mean[ref_order]

    max_abs = np.abs(fused_mean - ref_mean_reordered).max()
    print(f"n_groups={n_groups}, n_genes={len(all_genes)}, max_abs_diff={max_abs:.3e}")
    assert max_abs < 1e-5, f"fused mean diverged from sc.get.aggregate by {max_abs}"
    print("fused path numerically agrees with sc.get.aggregate")


if __name__ == "__main__":
    _check_api()
    test_backed_write()
    test_fused_aggregation()
```

- [ ] **Step 2: Run the sandbox**

```bash
/Users/enrique/projects/github/pyvatk/.venv/bin/python local/scripts/2026-04-17-ucsc-sandbox.py
```

Expected: both blocks print success. `test_backed_write` reports `shape=(520014, 200), X.format=csc`. `test_fused_aggregation` prints `max_abs_diff < 1e-5` and `fused path numerically agrees with sc.get.aggregate`.

If either fails, stop and debug before touching hvantk. The sandbox is the cheapest place to find an API mismatch.

- [ ] **Step 3: Clean up the sandbox artifact**

```bash
rm -f local/scripts/sandbox-backed.h5ad
```

Keep the `.py` script committed as a reproduction receipt (it's in `local/` which is gitignored — so nothing to commit).

---

## Task 2: Extract `_iter_ucsc_rows` from `create_anndata_from_ucsc_matrix`

**Files:**
- Modify: `hvantk/tables/ucsc.py` (existing `create_anndata_from_ucsc_matrix`, lines ~75–217)
- Regression: `hvantk/tests/test_expression_builders_anndata.py::TestBuildUcscAd` (must still pass).

Pure refactor. Lift the parse loop into a module-private helper so the backed builder and fused aggregator can reuse it. The existing `create_anndata_from_ucsc_matrix` becomes a thin wrapper around it.

- [ ] **Step 1: Add `_iter_ucsc_rows` helper in `hvantk/tables/ucsc.py`**

Insert this function after `_open_text` and before `load_ucsc_metadata`:

```python
from typing import Iterator


def _iter_ucsc_rows(
    expression_matrix_path: str,
    delimiter: str = "\t",
    split_gene_field: bool = True,
) -> tuple[list[str], Iterator[tuple[str, np.ndarray]]]:
    """Open a UCSC expression TSV (plain or gzipped) and return
    ``(cell_ids, row_iterator)``.

    ``cell_ids`` is the header's cell-column list, extracted eagerly.
    ``row_iterator`` yields ``(gene_name, row_values_float32)`` per data
    line and closes the underlying handle when exhausted.

    Raises
    ------
    FileNotFoundError
        If ``expression_matrix_path`` does not exist.
    ValueError
        If the file is empty or a row has the wrong number of values.
    """
    if not os.path.exists(expression_matrix_path):
        raise FileNotFoundError(
            f"Expression matrix file not found: {expression_matrix_path}"
        )

    fh = _open_text(expression_matrix_path)
    try:
        header = fh.readline().rstrip("\n").rstrip("\r")
        if not header:
            fh.close()
            raise ValueError(
                f"Expression matrix {expression_matrix_path!r} is empty."
            )
        header_fields = header.split(delimiter)
        cell_ids = header_fields[1:]
        n_cells = len(cell_ids)
    except BaseException:
        fh.close()
        raise

    def _rows() -> Iterator[tuple[str, np.ndarray]]:
        try:
            for line in fh:
                line = line.rstrip("\n").rstrip("\r")
                if not line:
                    continue
                tab = line.find(delimiter)
                if tab < 0:
                    continue
                gene = line[:tab]
                if split_gene_field:
                    gene = gene.split("|", 1)[0]
                row = np.fromstring(
                    line[tab + 1:], sep=delimiter, dtype=np.float32
                )
                if row.shape[0] != n_cells:
                    raise ValueError(
                        f"Row for gene {gene!r} has {row.shape[0]} parseable "
                        f"float values; expected {n_cells}. This usually "
                        f"means a short row or a non-numeric token in the row."
                    )
                yield gene, row
        finally:
            fh.close()

    return cell_ids, _rows()
```

- [ ] **Step 2: Rewrite `create_anndata_from_ucsc_matrix` to consume the helper**

Replace the body of `create_anndata_from_ucsc_matrix` (the `if not os.path.exists(...)` block through the closing `return adata`) with:

```python
    logger.info(
        "Streaming expression matrix from %s (chunk_size=%d)",
        expression_matrix_path,
        chunk_size,
    )

    cell_ids, row_iter = _iter_ucsc_rows(
        expression_matrix_path,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
    )

    gene_name_chunks: list[np.ndarray] = []
    sparse_chunks: list[sparse.csr_matrix] = []
    buf_rows: list[np.ndarray] = []
    buf_genes: list[str] = []
    n_seen = 0
    n_chunks = 0

    def _flush() -> None:
        nonlocal n_chunks
        if not buf_rows:
            return
        sparse_chunks.append(sparse.csr_matrix(np.vstack(buf_rows)))
        gene_name_chunks.append(np.asarray(buf_genes))
        buf_rows.clear()
        buf_genes.clear()
        n_chunks += 1
        if n_chunks % 10 == 0:
            logger.info("  streamed %d genes", n_seen)

    for gene, row in row_iter:
        buf_rows.append(row)
        buf_genes.append(gene)
        n_seen += 1
        if len(buf_rows) >= chunk_size:
            _flush()
    _flush()

    if not sparse_chunks:
        raise ValueError(
            f"Expression matrix {expression_matrix_path!r} contained no data rows."
        )

    X = sparse.vstack(sparse_chunks, format="csr").T.tocsr()
    gene_names = np.concatenate(gene_name_chunks)

    var = pd.DataFrame(index=pd.Index(gene_names, name=gene_column))
    obs = pd.DataFrame(index=pd.Index(cell_ids, name=UCSC_CELL_ID_COLUMN))

    adata = ad.AnnData(X=X, obs=obs, var=var)

    if metadata_df is not None:
        common = adata.obs.index.intersection(metadata_df.index)
        if len(common) == 0:
            raise ValueError(
                f"No overlapping cell ids between expression matrix "
                f"{expression_matrix_path!r} and metadata."
            )
        meta_aligned = metadata_df.reindex(adata.obs.index)
        for col in meta_aligned.columns:
            adata.obs[col] = meta_aligned[col].values

    logger.info("Built AnnData: %d cells × %d genes", adata.shape[0], adata.shape[1])
    return adata
```

Note the **spec-approved behavior flip**: zero-overlap cell IDs now raise `ValueError` instead of `logger.warning`.

- [ ] **Step 3: Run the regression test**

Run: `/Users/enrique/projects/github/pyvatk/.venv/bin/python -m pytest hvantk/tests/test_expression_builders_anndata.py::TestBuildUcscAd -v`

Expected: all tests PASS. If any failure, inspect — the refactor diverged from the original contract and must be corrected before continuing.

- [ ] **Step 4: Run lint on the modified file**

Run: `/Users/enrique/projects/github/pyvatk/.venv/bin/python -m flake8 hvantk/tables/ucsc.py --count --select=E9,F63,F7,F82 --show-source --statistics`

Expected: `0`.

- [ ] **Step 5: Propose commit to user**

Ask the user to approve:

```bash
git add hvantk/tables/ucsc.py
git commit -m "$(cat <<'EOF'
refactor(ucsc): extract _iter_ucsc_rows helper from builder

Lift the gzip/np.fromstring parse loop into a private generator so that
the upcoming backed builder and fused aggregator can reuse it without
duplicating validation. create_anndata_from_ucsc_matrix becomes a thin
wrapper. Also flip the zero-overlap-cell-ids behavior from warn-only to
ValueError (approved in design spec).
EOF
)"
```

Do not commit without explicit user approval per project policy.

---

## Task 3: Implement `build_ucsc_atlas_backed` (CSC backed writer)

**Files:**
- Modify: `hvantk/tables/ucsc.py` — add `build_ucsc_atlas_backed` after `create_anndata_from_ucsc_matrix`.
- Regression: extend `TestBuildUcscAd` parametrically in Task 7.

- [ ] **Step 1: Add the function**

Add these imports at the top of `hvantk/tables/ucsc.py` (alongside existing ones):

```python
import h5py
from anndata.io import sparse_dataset, write_elem
```

Append this function after `create_anndata_from_ucsc_matrix`:

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
    uns: "dict | None" = None,
) -> str:
    """Stream-build an AnnData .h5ad file on disk, appending one batch of
    genes (CSC columns) at a time. Never materializes the full
    ``cells × genes`` matrix.

    Parameters
    ----------
    expression_matrix_path : str
        Path to the UCSC expression TSV (plain or gzipped).
    output_path : str
        Destination ``.h5ad`` path.
    metadata_df : pd.DataFrame, optional
        Cell metadata; reindexed to expression header ``cell_ids``.
    gene_column : str
        Name for the ``var`` index.
    delimiter : str
        Column delimiter in the expression TSV.
    split_gene_field : bool
        If True, split gene ids on ``|`` and keep the first element.
    column_batch : int
        Number of gene columns buffered before a CSC append to disk.
    overwrite : bool
        If False and ``output_path`` exists, raise ``FileExistsError``.
    uns : dict, optional
        Opaque dict written to the h5ad's ``uns`` group (e.g. provenance
        metadata). Callers remain responsible for the dict's structure.

    Returns
    -------
    str
        ``output_path``.
    """
    if os.path.exists(output_path) and not overwrite:
        raise FileExistsError(
            f"{output_path!r} already exists; pass overwrite=True to replace."
        )

    cell_ids, row_iter = _iter_ucsc_rows(
        expression_matrix_path,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
    )
    n_cells = len(cell_ids)
    logger.info(
        "Backed-write atlas → %s (n_cells=%d, column_batch=%d)",
        output_path, n_cells, column_batch,
    )

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    if os.path.exists(output_path):
        os.remove(output_path)

    f = h5py.File(output_path, "w")
    try:
        # Seed X as an empty (n_cells, 0) CSC; pass indptr dtype via
        # dataset_kwargs so large atlases don't overflow int32.
        seed = sparse.csc_matrix((n_cells, 0), dtype=np.float32)
        write_elem(f, "X", seed, dataset_kwargs={"indptr_dtype": "int64"})
        X = sparse_dataset(f["X"])

        gene_names: list[str] = []
        buf_cols: list[np.ndarray] = []
        buf_genes: list[str] = []
        n_appended = 0
        n_batches = 0

        def _flush() -> None:
            nonlocal n_appended, n_batches
            if not buf_cols:
                return
            # Stack column-vectors into (n_cells, batch_width) CSC block.
            dense_block = np.column_stack(buf_cols)  # (n_cells, batch_width)
            X.append(sparse.csc_matrix(dense_block))
            gene_names.extend(buf_genes)
            n_appended += len(buf_cols)
            n_batches += 1
            buf_cols.clear()
            buf_genes.clear()
            if n_batches % 20 == 0:
                logger.info("  appended %d genes", n_appended)

        for gene, row in row_iter:
            buf_cols.append(row)
            buf_genes.append(gene)
            if len(buf_cols) >= column_batch:
                _flush()
        _flush()

        if n_appended == 0:
            raise ValueError(
                f"Expression matrix {expression_matrix_path!r} contained no data rows."
            )

        # obs: cell_ids + optional metadata join (reject zero-overlap).
        obs = pd.DataFrame(index=pd.Index(cell_ids, name=UCSC_CELL_ID_COLUMN))
        if metadata_df is not None:
            common = obs.index.intersection(metadata_df.index)
            if len(common) == 0:
                raise ValueError(
                    f"No overlapping cell ids between expression matrix "
                    f"{expression_matrix_path!r} and metadata."
                )
            meta_aligned = metadata_df.reindex(obs.index)
            for col in meta_aligned.columns:
                obs[col] = meta_aligned[col].values

        var = pd.DataFrame(index=pd.Index(gene_names, name=gene_column))

        write_elem(f, "obs", obs)
        write_elem(f, "var", var)
        if uns is not None:
            write_elem(f, "uns", uns)
    finally:
        f.close()

    logger.info(
        "Backed atlas written: %d cells × %d genes → %s",
        n_cells, n_appended, output_path,
    )
    return output_path
```

- [ ] **Step 2: Smoke-test with a tiny synthetic input**

Run:

```bash
/Users/enrique/projects/github/pyvatk/.venv/bin/python - <<'PY'
import gzip, os, tempfile
import numpy as np
import pandas as pd
from hvantk.tables.ucsc import build_ucsc_atlas_backed
import anndata as ad

rng = np.random.default_rng(0)
n_genes, n_cells = 50, 300
vals = rng.random((n_genes, n_cells)).astype(np.float32)
vals[vals < 0.4] = 0.0
genes = [f"G{i}|ENST{i:05d}" for i in range(n_genes)]
cells = [f"C{j}" for j in range(n_cells)]

with tempfile.TemporaryDirectory() as td:
    src = os.path.join(td, "x.tsv.gz")
    with gzip.open(src, "wt") as fh:
        fh.write("gene\t" + "\t".join(cells) + "\n")
        for g, row in zip(genes, vals):
            fh.write(g + "\t" + "\t".join(f"{v:.6g}" for v in row) + "\n")

    meta = pd.DataFrame({"cell_type": ["A"]*150 + ["B"]*150},
                        index=pd.Index(cells, name="cell_id"))
    out = os.path.join(td, "atlas.h5ad")
    build_ucsc_atlas_backed(src, out, metadata_df=meta, column_batch=8)
    a = ad.read_h5ad(out)
    print(f"shape={a.shape}, X.format={a.X.format if hasattr(a.X,'format') else 'dense'}")
    print(f"obs.columns={list(a.obs.columns)}, var_names[:3]={list(a.var_names[:3])}")
    assert a.shape == (n_cells, n_genes), a.shape
    # Numerical equivalence vs. in-memory path
    from hvantk.tables.ucsc import create_anndata_from_ucsc_matrix
    b = create_anndata_from_ucsc_matrix(src, metadata_df=meta, chunk_size=10)
    X_a = a.X.toarray() if hasattr(a.X, "toarray") else a.X
    X_b = b.X.toarray() if hasattr(b.X, "toarray") else b.X
    diff = np.abs(X_a - X_b).max()
    print(f"max_abs_diff vs in-memory builder: {diff:.3e}")
    assert diff < 1e-6, diff
print("OK")
PY
```

Expected: `OK` printed, `shape=(300, 50)`, `max_abs_diff < 1e-6`.

- [ ] **Step 3: Run lint**

Run: `/Users/enrique/projects/github/pyvatk/.venv/bin/python -m flake8 hvantk/tables/ucsc.py --count --select=E9,F63,F7,F82 --show-source --statistics`

Expected: `0`.

- [ ] **Step 4: Propose commit to user**

```bash
git add hvantk/tables/ucsc.py
git commit -m "$(cat <<'EOF'
feat(ucsc): add build_ucsc_atlas_backed with CSC incremental writes

Streams UCSC gz → one CSC column block at a time into an h5ad on disk via
anndata.io.sparse_dataset and CSCDataset.append. Passes indptr_dtype=int64
so 500k+-cell atlases don't overflow the default int32. Never materializes
the full cells x genes matrix; peak RAM is column_batch x n_cells x 4B.
EOF
)"
```

---

## Task 4: Implement `summarize_ucsc_streaming` (fused aggregator)

**Files:**
- Modify: `hvantk/tables/ucsc.py` — add `summarize_ucsc_streaming` after `build_ucsc_atlas_backed`.

- [ ] **Step 1: Add the function**

Append to `hvantk/tables/ucsc.py`:

```python
def summarize_ucsc_streaming(
    expression_matrix_path: str,
    metadata_df: pd.DataFrame,
    group_by: "str | list[str]",
    filter_by: "dict[str, str | list[str]] | None" = None,
    min_cells_per_group: int = 10,
    gene_column: str = UCSC_GENE_COLUMN,
    delimiter: str = "\t",
    split_gene_field: bool = True,
) -> ad.AnnData:
    """Stream a UCSC expression matrix and aggregate per-group per-gene
    statistics in one pass.

    Returns an AnnData of shape ``(n_groups, n_genes)`` with ``X`` set to
    the per-group mean (float32) and ``layers`` containing ``sum``,
    ``count_nonzero``, ``fraction_expressed``. ``obs`` includes the
    original group-by columns plus ``n_cells``.

    Memory: one ``float32[n_cells]`` row buffer + accumulators of size
    ``n_groups × n_genes × 16 B`` (sum float64 + count_nz int64). Independent
    of ``n_cells`` past the single-row buffer.
    """
    by = [group_by] if isinstance(group_by, str) else list(group_by)

    # Apply filter on the DataFrame before factorizing groups.
    work = metadata_df.copy()
    if filter_by:
        for field, value in filter_by.items():
            if field not in work.columns:
                raise ValueError(
                    f"filter_by field {field!r} not in metadata columns: "
                    f"{sorted(work.columns)}"
                )
            if isinstance(value, (list, tuple, set)):
                work = work[work[field].isin(list(value))]
            else:
                work = work[work[field] == value]
        if work.empty:
            raise ValueError("filter_by produced zero cells.")

    for col in by:
        if col not in work.columns:
            raise ValueError(
                f"group_by column {col!r} not in metadata columns: "
                f"{sorted(work.columns)}"
            )

    cell_ids, row_iter = _iter_ucsc_rows(
        expression_matrix_path,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
    )
    n_cells = len(cell_ids)

    # Align metadata to expression header order, keep only cells present
    # in both, and drop NaN-in-group rows.
    aligned = work.reindex(cell_ids)
    keep_mask = aligned[by].notna().all(axis=1).to_numpy()
    if not keep_mask.any():
        raise ValueError(
            "No cells overlap between metadata (post-filter) and expression header."
        )
    n_dropped_nan = int((~keep_mask[aligned.index.isin(work.index)]).sum())
    if n_dropped_nan:
        example = cell_ids[int(np.where(~keep_mask)[0][0])]
        logger.info(
            "Dropped %d cells with NaN in group-by columns (example: %s).",
            n_dropped_nan, example,
        )

    labels = (
        aligned.loc[keep_mask, by].astype(str).agg("_".join, axis=1)
        if len(by) > 1
        else aligned.loc[keep_mask, by[0]].astype(str)
    )
    codes, uniques = pd.factorize(labels, sort=True)
    n_groups = len(uniques)

    # group_idx per cell in expression-header order, -1 for dropped cells.
    group_idx = np.full(n_cells, -1, dtype=np.int64)
    group_idx[keep_mask] = codes
    valid = group_idx >= 0
    n_cells_per_group = np.bincount(group_idx[valid], minlength=n_groups).astype(np.int64)

    sum_matrix = np.zeros((n_groups, 0), dtype=np.float64)
    count_nz_matrix = np.zeros((n_groups, 0), dtype=np.int64)
    gene_names: list[str] = []

    # Pre-size buffers in chunks to avoid per-gene resize; column-extend
    # accumulators lazily.
    BLOCK = 1024
    sum_block = np.zeros((n_groups, BLOCK), dtype=np.float64)
    count_block = np.zeros((n_groups, BLOCK), dtype=np.int64)
    block_fill = 0

    def _flush_block():
        nonlocal sum_matrix, count_nz_matrix, block_fill
        if block_fill == 0:
            return
        sum_matrix = np.concatenate([sum_matrix, sum_block[:, :block_fill]], axis=1)
        count_nz_matrix = np.concatenate(
            [count_nz_matrix, count_block[:, :block_fill]], axis=1
        )
        block_fill = 0

    for gene, row in row_iter:
        row_valid = row[valid]
        sum_block[:, block_fill] = np.bincount(
            group_idx[valid], weights=row_valid, minlength=n_groups
        )
        count_block[:, block_fill] = np.bincount(
            group_idx[valid],
            weights=(row_valid != 0).astype(np.float64),
            minlength=n_groups,
        ).astype(np.int64)
        gene_names.append(gene)
        block_fill += 1
        if block_fill == BLOCK:
            _flush_block()
    _flush_block()

    if len(gene_names) == 0:
        raise ValueError(
            f"Expression matrix {expression_matrix_path!r} contained no data rows."
        )

    n_cells_col = n_cells_per_group[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = np.where(n_cells_col > 0, sum_matrix / n_cells_col, 0.0).astype(np.float32)
        fraction_expressed = np.where(
            n_cells_col > 0, count_nz_matrix / n_cells_col, 0.0
        ).astype(np.float32)

    obs = pd.DataFrame({"n_cells": n_cells_per_group}, index=pd.Index(uniques, name="group"))
    # Unpack the composite label back into per-field columns.
    if len(by) == 1:
        obs[by[0]] = uniques
    else:
        parts = [u.split("_") for u in uniques]
        for i, col in enumerate(by):
            obs[col] = [p[i] if i < len(p) else "" for p in parts]

    var = pd.DataFrame(index=pd.Index(gene_names, name=gene_column))

    keep_groups = obs["n_cells"].to_numpy() >= min_cells_per_group
    if not keep_groups.any():
        # Build a short report for the error message.
        report = "; ".join(
            f"{u}={int(n)}" for u, n in zip(uniques, n_cells_per_group)
        )
        raise ValueError(
            f"All groups fell below min_cells_per_group={min_cells_per_group}. "
            f"Observed: {report}"
        )

    adata = ad.AnnData(
        X=mean[keep_groups],
        obs=obs[keep_groups].copy(),
        var=var,
        layers={
            "sum": sum_matrix[keep_groups].astype(np.float32),
            "count_nonzero": count_nz_matrix[keep_groups].astype(np.int64),
            "fraction_expressed": fraction_expressed[keep_groups],
            "mean": mean[keep_groups],
        },
    )
    logger.info(
        "Aggregated AnnData: %d groups × %d genes (dropped %d groups below min_cells)",
        adata.n_obs, adata.n_vars, int((~keep_groups).sum()),
    )
    return adata
```

- [ ] **Step 2: Smoke-test numerical equivalence against `summarize_expression_ad`**

Run:

```bash
/Users/enrique/projects/github/pyvatk/.venv/bin/python - <<'PY'
import gzip, os, tempfile
import numpy as np
import pandas as pd
import anndata as ad
from hvantk.tables.ucsc import summarize_ucsc_streaming, create_anndata_from_ucsc_matrix
from hvantk.utils.matrix_utils import summarize_expression_ad

rng = np.random.default_rng(42)
n_genes, n_cells = 40, 300
vals = rng.random((n_genes, n_cells)).astype(np.float32)
vals[vals < 0.5] = 0.0
genes = [f"G{i}" for i in range(n_genes)]
cells = [f"C{j}" for j in range(n_cells)]

with tempfile.TemporaryDirectory() as td:
    src = os.path.join(td, "x.tsv.gz")
    with gzip.open(src, "wt") as fh:
        fh.write("gene\t" + "\t".join(cells) + "\n")
        for g, row in zip(genes, vals):
            fh.write(g + "\t" + "\t".join(f"{v:.6g}" for v in row) + "\n")

    meta = pd.DataFrame(
        {"cell_type": ["A"] * 120 + ["B"] * 100 + ["C"] * 80},
        index=pd.Index(cells, name="cell_id"),
    )

    fused = summarize_ucsc_streaming(
        src, meta, group_by="cell_type", min_cells_per_group=5
    )

    ad_full = create_anndata_from_ucsc_matrix(src, metadata_df=meta, chunk_size=10)
    two_step = summarize_expression_ad(ad_full, group_by="cell_type", min_cells_per_group=5)

    # Align on group label order.
    fused = fused[sorted(fused.obs_names)].copy()
    two_step = two_step[sorted(two_step.obs_names)].copy()

    assert list(fused.obs_names) == list(two_step.obs_names), (
        list(fused.obs_names), list(two_step.obs_names)
    )
    assert fused.shape == two_step.shape, (fused.shape, two_step.shape)
    for layer in ("mean", "fraction_expressed"):
        a = np.asarray(fused.layers[layer])
        b = np.asarray(two_step.layers[layer])
        diff = np.abs(a - b).max()
        print(f"{layer}: max_abs_diff={diff:.3e}")
        assert diff < 1e-5, (layer, diff)
print("OK")
PY
```

Expected: `OK`. All layers (`mean`, `fraction_expressed`) agree within `1e-5`.

- [ ] **Step 3: Run lint**

Run: `/Users/enrique/projects/github/pyvatk/.venv/bin/python -m flake8 hvantk/tables/ucsc.py --count --select=E9,F63,F7,F82 --show-source --statistics`

Expected: `0`.

- [ ] **Step 4: Propose commit to user**

```bash
git add hvantk/tables/ucsc.py
git commit -m "$(cat <<'EOF'
feat(ucsc): add summarize_ucsc_streaming fused aggregator

Streams UCSC gz row-by-row, reducing each row into (n_groups, n_genes)
float64/int64 accumulators via np.bincount. Emits an AnnData whose
schema (obs[n_cells], layers {mean, sum, count_nonzero,
fraction_expressed}) matches summarize_expression_ad so downstream
consumers cannot tell which path produced a given summary.
EOF
)"
```

---

## Task 5: Plumb `backed` kwarg through `build_ucsc_ad`

**Files:**
- Modify: `hvantk/tables/matrix_builders.py` — function `build_ucsc_ad`.

- [ ] **Step 1: Add `backed` kwarg and dispatch logic**

In `hvantk/tables/matrix_builders.py`, find `build_ucsc_ad` (around line 24). Add a new constant at the top of the module near other constants (or add one if none exist yet):

```python
# Auto-select the backed builder for UCSC inputs larger than this threshold.
BACKED_BUILDER_THRESHOLD_BYTES = 1 * 1024 * 1024 * 1024  # 1 GiB
```

Update `build_ucsc_ad` to the signature:

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
    backed: "bool | None" = None,
    column_batch: int = 64,
) -> "ad.AnnData":
```

Replace the body with:

```python
    import anndata as ad

    from hvantk.tables.ucsc import (
        load_ucsc_metadata,
        create_anndata_from_ucsc_matrix,
        build_ucsc_atlas_backed,
    )
    from hvantk.core.anndata_utils import (
        build_anndata_metadata,
        annotate_column_summary_ad,
        save_anndata,
    )

    logger.info("Loading UCSC metadata from %s", metadata_path)
    metadata_df = load_ucsc_metadata(metadata_path, sep=delimiter)

    # Auto-select backed mode by input size unless caller forces.
    if backed is None:
        size = os.path.getsize(expression_matrix_path)
        backed = size > BACKED_BUILDER_THRESHOLD_BYTES
        logger.info(
            "build_ucsc_ad: auto-selected backed=%s (input size %.2f GiB, threshold %.2f GiB)",
            backed, size / (1024**3), BACKED_BUILDER_THRESHOLD_BYTES / (1024**3),
        )

    if backed:
        if output_path is None:
            raise ValueError("backed=True requires output_path (writes directly to disk).")
        provenance = {"hvantk_metadata": build_anndata_metadata("UCSC", expression_matrix_path)}
        build_ucsc_atlas_backed(
            expression_matrix_path=expression_matrix_path,
            output_path=output_path,
            metadata_df=metadata_df,
            gene_column=gene_column,
            delimiter=delimiter,
            split_gene_field=split_gene_field,
            column_batch=column_batch,
            overwrite=overwrite,
            uns=provenance,
        )
        # Return a backed-mode handle — shape + obs/var without materializing X.
        # annotate_column_summary_ad would need to scan the full X matrix and
        # is intentionally skipped for backed atlases (v1 trade-off).
        logger.info(
            "Backed atlas built at %s; skipping annotate_column_summary_ad "
            "(would materialize X). Returning a backed AnnData handle.",
            output_path,
        )
        return ad.read_h5ad(output_path, backed="r")

    logger.info("Creating AnnData from UCSC expression matrix (in-memory)")
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

Add `import os` at the top of the module if not already present.

- [ ] **Step 2: Run existing builder tests**

Run: `/Users/enrique/projects/github/pyvatk/.venv/bin/python -m pytest hvantk/tests/test_expression_builders_anndata.py::TestBuildUcscAd -v`

Expected: PASS.

- [ ] **Step 3: Run lint**

Run: `/Users/enrique/projects/github/pyvatk/.venv/bin/python -m flake8 hvantk/tables/matrix_builders.py --count --select=E9,F63,F7,F82 --show-source --statistics`

Expected: `0`.

- [ ] **Step 4: Propose commit to user**

```bash
git add hvantk/tables/matrix_builders.py
git commit -m "$(cat <<'EOF'
feat(ucsc): route build_ucsc_ad through backed builder for large inputs

Add backed: bool | None kwarg (default auto-select by input size >
BACKED_BUILDER_THRESHOLD_BYTES, 1 GiB). Dispatches to
build_ucsc_atlas_backed when backed, create_anndata_from_ucsc_matrix
otherwise. backed=True requires output_path.
EOF
)"
```

---

## Task 6: Expose `--backed / --no-backed` on `mkmatrix ucsc`

**Files:**
- Modify: `hvantk/commands/make_matrix_cli.py` — command `mkmatrix_ucsc`.

- [ ] **Step 1: Add the flag**

Find `mkmatrix_ucsc` in `hvantk/commands/make_matrix_cli.py`. Add a new click option after `--chunk-size`:

```python
@click.option(
    "--backed/--no-backed",
    default=None,
    help=(
        "Use the incremental backed-write builder (CSC columns appended to "
        ".h5ad on disk; no RAM materialization). Default: auto — backed if "
        "the input is larger than ~1 GiB."
    ),
)
@click.option(
    "--column-batch",
    type=int,
    default=64,
    show_default=True,
    help="Gene columns per backed-append block (backed mode only).",
)
```

Add `backed` and `column_batch` params to the `mkmatrix_ucsc` signature and pass them through to `_build_ucsc_ad`:

```python
def mkmatrix_ucsc(
    expression_matrix,
    metadata,
    output_path,
    gene_column,
    delimiter,
    split_gene_field,
    chunk_size,
    backed,
    column_batch,
    overwrite,
):
    ...
    adata = _build_ucsc_ad(
        expression_matrix_path=expression_matrix,
        metadata_path=metadata,
        output_path=output_path,
        gene_column=gene_column,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
        chunk_size=chunk_size,
        backed=backed,
        column_batch=column_batch,
        overwrite=overwrite,
    )
```

- [ ] **Step 2: Run existing CLI test**

Run: `/Users/enrique/projects/github/pyvatk/.venv/bin/python -m pytest hvantk/tests/test_expression_builders_anndata.py::TestMkmatrixCli::test_ucsc_produces_h5ad -v`

Expected: PASS.

- [ ] **Step 3: Sanity-check help output**

Run: `/Users/enrique/projects/github/pyvatk/.venv/bin/hvantk mkmatrix ucsc --help`

Expected: see `--backed/--no-backed` and `--column-batch` in the help text.

- [ ] **Step 4: Run lint**

Run: `/Users/enrique/projects/github/pyvatk/.venv/bin/python -m flake8 hvantk/commands/make_matrix_cli.py --count --select=E9,F63,F7,F82 --show-source --statistics`

Expected: `0`.

- [ ] **Step 5: Propose commit to user**

```bash
git add hvantk/commands/make_matrix_cli.py
git commit -m "feat(mkmatrix): add --backed/--no-backed + --column-batch to mkmatrix ucsc"
```

---

## Task 7: Parametrize `TestBuildUcscAd` with `backed=[False, True]`

**Files:**
- Modify: `hvantk/tests/test_expression_builders_anndata.py` — class `TestBuildUcscAd`.

- [ ] **Step 1: Inspect the current test class**

Run: `/Users/enrique/projects/github/pyvatk/.venv/bin/python -m pytest hvantk/tests/test_expression_builders_anndata.py::TestBuildUcscAd -v --collect-only`

Read the existing test methods and their fixtures. Note how `build_ucsc_ad` is invoked (positional vs kwargs) and which paths are expected.

- [ ] **Step 2: Parametrize**

At the top of `TestBuildUcscAd`, add a class-level parametrize:

```python
@pytest.mark.parametrize("backed", [False, True])
class TestBuildUcscAd:
    ...
```

In every test method, add `backed` as a parameter and pass it through to `build_ucsc_ad`:

```python
def test_...(self, tmp_path, backed):
    ...
    adata = build_ucsc_ad(
        ...,
        backed=backed,
        output_path=str(tmp_path / "atlas.h5ad"),  # required when backed=True
        overwrite=True,
    )
    ...
```

If any test currently passes `output_path=None`, force `output_path=str(tmp_path / "...h5ad")` so the backed variant can write. The in-memory variant tolerates a real path too.

- [ ] **Step 3: Run the parametrized tests**

Run: `/Users/enrique/projects/github/pyvatk/.venv/bin/python -m pytest hvantk/tests/test_expression_builders_anndata.py::TestBuildUcscAd -v`

Expected: test count doubles (each existing test now has `[False]` and `[True]` variants). All PASS.

- [ ] **Step 4: Run the full fast-test suite to check for collateral damage**

Run: `/Users/enrique/projects/github/pyvatk/.venv/bin/python -m pytest -q`

Expected: all PASS (excluding hail/network/slow markers as configured in `pytest.ini`).

- [ ] **Step 5: Propose commit to user**

```bash
git add hvantk/tests/test_expression_builders_anndata.py
git commit -m "test(ucsc): parametrize TestBuildUcscAd with backed=[False, True]"
```

---

## Task 8: Add `expression summarize-ucsc` CLI subcommand

**Files:**
- Modify: `hvantk/commands/summarize_expression_cli.py` — add `summarize_ucsc_cmd`.

- [ ] **Step 1: Add the command**

At the bottom of `hvantk/commands/summarize_expression_cli.py`, before any `__main__` stanza, add:

```python
@expression_group.command("summarize-ucsc")
@click.option(
    "-e",
    "--expression-matrix",
    "expression_matrix",
    type=click.Path(exists=True),
    required=True,
    help="Path to the UCSC expression TSV (plain or gzipped).",
)
@click.option(
    "-m",
    "--metadata",
    "metadata_path",
    type=click.Path(exists=True),
    required=True,
    help="Path to the UCSC cell metadata TSV.",
)
@click.option(
    "--group-by",
    multiple=True,
    required=True,
    help="Metadata column(s) to group by. Repeat for multi-field grouping.",
)
@click.option(
    "--filter-by",
    multiple=True,
    default=None,
    help=(
        "Pre-filter cells: FIELD=VALUE (repeatable). "
        "Example: --filter-by Region=Cortex --filter-by TimePoint=9wpc"
    ),
)
@click.option(
    "--min-cells",
    type=int,
    default=10,
    show_default=True,
    help="Drop groups with fewer cells than this threshold.",
)
@click.option(
    "--gene-column",
    default="gene",
    show_default=True,
)
@click.option(
    "--delimiter",
    default="\t",
    show_default=True,
)
@click.option(
    "--split-gene-field/--no-split-gene-field",
    default=True,
    show_default=True,
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
def summarize_ucsc_cmd(
    expression_matrix,
    metadata_path,
    group_by,
    filter_by,
    min_cells,
    gene_column,
    delimiter,
    split_gene_field,
    output,
    overwrite,
):
    """Fused stream-and-aggregate: UCSC expression + metadata → groups × genes .h5ad.

    Skips the intermediate cells × genes atlas.h5ad — streams the expression
    matrix row-by-row and accumulates per-group per-gene statistics in a
    single pass. Output schema matches `hvantk expression summarize`, so
    downstream consumers do not branch on which path produced the summary.

    \b
    Examples:

      # Single-field grouping
      hvantk expression summarize-ucsc \\
          -e exprMatrix.tsv.gz -m meta.tsv \\
          --group-by Class \\
          -o class_summary.h5ad

      # Multi-field grouping + pre-filter
      hvantk expression summarize-ucsc \\
          -e exprMatrix.tsv.gz -m meta.tsv \\
          --group-by Region --group-by TimePoint \\
          --filter-by Region=Cortex \\
          --min-cells 10 \\
          -o region_timepoint_summary.h5ad
    """
    from pathlib import Path

    from hvantk.tables.ucsc import load_ucsc_metadata, summarize_ucsc_streaming

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
                    f"Expected FIELD=VALUE format, got: {item!r}",
                    param_hint="--filter-by",
                )
            key, value = item.split("=", 1)
            filters[key.strip()] = value.strip()

    metadata_df = load_ucsc_metadata(metadata_path, sep=delimiter)

    missing = [col for col in group_by if col not in metadata_df.columns]
    if missing:
        raise click.BadParameter(
            f"--group-by column(s) not in metadata: {missing}. "
            f"Available: {sorted(metadata_df.columns)}",
            param_hint="--group-by",
        )

    summary = summarize_ucsc_streaming(
        expression_matrix_path=expression_matrix,
        metadata_df=metadata_df,
        group_by=list(group_by),
        filter_by=filters,
        min_cells_per_group=min_cells,
        gene_column=gene_column,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
    )

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

- [ ] **Step 2: Smoke-test the new command**

Run:

```bash
/Users/enrique/projects/github/pyvatk/.venv/bin/python - <<'PY'
import gzip, subprocess, tempfile, pathlib
import numpy as np, pandas as pd, anndata as ad

rng = np.random.default_rng(0)
n_genes, n_cells = 30, 100
vals = rng.random((n_genes, n_cells)).astype(np.float32); vals[vals < 0.5] = 0.0
genes = [f"G{i}" for i in range(n_genes)]
cells = [f"C{j}" for j in range(n_cells)]

with tempfile.TemporaryDirectory() as d:
    d = pathlib.Path(d)
    src = d / "x.tsv.gz"
    with gzip.open(src, "wt") as fh:
        fh.write("gene\t" + "\t".join(cells) + "\n")
        for g, row in zip(genes, vals):
            fh.write(g + "\t" + "\t".join(f"{v:.6g}" for v in row) + "\n")
    meta = pd.DataFrame(
        {"cell_type": ["A"]*40 + ["B"]*35 + ["C"]*25},
        index=pd.Index(cells, name="cell_id"),
    )
    meta_path = d / "meta.tsv"
    meta.to_csv(meta_path, sep="\t")

    out = d / "summary.h5ad"
    r = subprocess.run(
        ["hvantk", "expression", "summarize-ucsc",
         "-e", str(src), "-m", str(meta_path),
         "--group-by", "cell_type",
         "--min-cells", "5",
         "-o", str(out)],
        capture_output=True, text=True,
    )
    print("stdout:", r.stdout)
    print("stderr:", r.stderr)
    assert r.returncode == 0, r.returncode
    loaded = ad.read_h5ad(out)
    print("shape:", loaded.shape, "layers:", sorted(loaded.layers.keys()))
    assert set(loaded.layers.keys()) >= {"mean", "sum", "count_nonzero", "fraction_expressed"}
    assert loaded.shape == (3, n_genes)
print("OK")
PY
```

Expected: `OK`, `shape: (3, 30)`, all 4 layers present.

- [ ] **Step 3: Smoke-test error handling (bad --group-by)**

```bash
/Users/enrique/projects/github/pyvatk/.venv/bin/python - <<'PY'
import subprocess, tempfile, pathlib, gzip
import numpy as np, pandas as pd
rng = np.random.default_rng(0)
with tempfile.TemporaryDirectory() as d:
    d = pathlib.Path(d)
    src = d / "x.tsv.gz"
    with gzip.open(src, "wt") as fh:
        fh.write("gene\tC1\tC2\n")
        fh.write("G1\t1.0\t2.0\n")
    meta_path = d / "meta.tsv"
    pd.DataFrame({"cell_type": ["A", "B"]}, index=pd.Index(["C1", "C2"], name="cell_id")).to_csv(meta_path, sep="\t")
    r = subprocess.run(
        ["hvantk", "expression", "summarize-ucsc",
         "-e", str(src), "-m", str(meta_path),
         "--group-by", "NoSuchField",
         "-o", str(d / "out.h5ad")],
        capture_output=True, text=True,
    )
    print("exit:", r.returncode)
    print("stderr:", r.stderr[:300])
    assert r.returncode == 2, r.returncode
    assert "NoSuchField" in r.stderr
    assert "cell_type" in r.stderr
print("OK — bad --group-by fails fast with helpful message")
PY
```

Expected: exit 2, error message lists `NoSuchField` and `cell_type`.

- [ ] **Step 4: Run full fast-test suite**

Run: `/Users/enrique/projects/github/pyvatk/.venv/bin/python -m pytest -q`

Expected: all PASS.

- [ ] **Step 5: Run lint**

Run: `/Users/enrique/projects/github/pyvatk/.venv/bin/python -m flake8 hvantk/commands/summarize_expression_cli.py --count --select=E9,F63,F7,F82 --show-source --statistics`

Expected: `0`.

- [ ] **Step 6: Propose commit to user**

```bash
git add hvantk/commands/summarize_expression_cli.py
git commit -m "$(cat <<'EOF'
feat(expression): add summarize-ucsc CLI (fused stream-and-aggregate)

New hvantk expression summarize-ucsc subcommand: takes raw UCSC
expression + metadata + --group-by, streams and aggregates in one pass
via summarize_ucsc_streaming, emits a .h5ad whose schema matches
expression summarize. Skips the cells x genes atlas intermediate.
EOF
)"
```

---

## Task 9: Manual validation on the adult-cortex meta-atlas

No automated tests. Purpose: verify end-to-end behavior on real data before PR.

- [ ] **Step 1: Backed build of the atlas**

Run (in a separate terminal — it will take ~30–45 min):

```bash
/Users/enrique/projects/github/pyvatk/.venv/bin/hvantk -v mkmatrix ucsc \
    -e /Users/enrique/projects/github/pyvatk/local/data/UCSC/adult-ctx-meta-atlas/adult-ctx-meta-atlas/exprMatrix.tsv.gz \
    -m /Users/enrique/projects/github/pyvatk/local/data/UCSC/adult-ctx-meta-atlas/adult-ctx-meta-atlas/meta.tsv \
    -o /Users/enrique/projects/github/pyvatk/local/data/UCSC/adult-ctx-meta-atlas/adult-ctx-meta-atlas/adult_ctx_meta_atlas.h5ad \
    --backed \
    --column-batch 64 \
    --overwrite
```

Acceptance:
- Finishes without OOM.
- Peak RSS ≤ ~3 GB (spot-check via Activity Monitor or `ps`).
- Log lines show periodic `appended N genes` progress.
- Resulting file reopens cleanly: `adata = ad.read_h5ad(...)`; `adata.shape == (520014, 122697)`.

If peak RSS is above ~5 GB, investigate before continuing — the backed path is supposed to be bounded.

- [ ] **Step 2: Fused summarize of the same raw files**

```bash
/Users/enrique/projects/github/pyvatk/.venv/bin/hvantk -v expression summarize-ucsc \
    -e /Users/enrique/projects/github/pyvatk/local/data/UCSC/adult-ctx-meta-atlas/adult-ctx-meta-atlas/exprMatrix.tsv.gz \
    -m /Users/enrique/projects/github/pyvatk/local/data/UCSC/adult-ctx-meta-atlas/adult-ctx-meta-atlas/meta.tsv \
    --group-by Class \
    --min-cells 10 \
    -o /Users/enrique/projects/github/pyvatk/output/ptm_assembly/brain_metaatlas_gene_class_mean_expr.h5ad \
    --overwrite
```

Acceptance:
- Finishes in ≤ 45 min wall-clock.
- `summary.obs["n_cells"].sum()` equals the number of non-NaN `Class` cells in metadata.
- Report card shows ~10 groups (Jorstad Class labels), ~122k transcripts, layers `mean`, `sum`, `count_nonzero`, `fraction_expressed`.

- [ ] **Step 3: Numerical spot-check against the prior pickle (if present)**

Only if `output/ptm_assembly/brain_metaatlas_gene_class_mean_expr.pkl` still exists:

```bash
/Users/enrique/projects/github/pyvatk/.venv/bin/python - <<'PY'
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
print("OK — new fused pipeline agrees with old pickle within tolerance")
PY
```

Expected: `OK — new fused pipeline agrees with old pickle within tolerance`.

- [ ] **Step 4: Record the run as a validation receipt**

In `local/notebooks/` or `local/planning/`, create a short markdown note with: date, dataset, wall-clock per run, peak RSS, final shape, commit SHA at the time of the run. Commit it if it lives in a non-gitignored path; otherwise keep it local.

---

## Task 10: Open the PR

- [ ] **Step 1: Confirm with user that the branch is ready**

Ask: "Manual validation passed on `adult-ctx-meta-atlas`. Ready to push `feat/sc-aggregation-cli` and open a PR into `dev`?"

- [ ] **Step 2: Push branch (only after user says yes)**

```bash
git push -u origin feat/sc-aggregation-cli
```

- [ ] **Step 3: Create PR**

```bash
gh pr create --base dev --title "feat(expression): backed-write mkmatrix ucsc + fused summarize-ucsc" --body "$(cat <<'EOF'
## Summary

- `mkmatrix ucsc` can now build `atlas.h5ad` without materializing the full cells × genes matrix. Uses CSC incremental writes via `anndata.io.sparse_dataset`; `--backed` auto-enabled when the input is > 1 GiB.
- New `hvantk expression summarize-ucsc` subcommand: streams UCSC expression + metadata and aggregates in a single `np.bincount` pass, emitting a `groups × genes` `.h5ad` whose schema matches the existing `expression summarize`.
- Shared low-level row streamer (`_iter_ucsc_rows`) extracted from the existing builder so both new paths reuse the same validated parser.

## Test plan

- [x] `pytest hvantk/tests/test_expression_builders_anndata.py::TestBuildUcscAd -v` (parametrized on `backed=[False, True]`)
- [x] `pytest hvantk/tests/test_matrix_utils_anndata.py::TestSummarizeExpressionAd -v` (output schema parity)
- [x] `pytest -q` full fast suite
- [x] Manual backed build of `adult-ctx-meta-atlas` (520k × 122k) finished under 45 min; peak RSS ≤ 3 GB
- [x] Fused `summarize-ucsc` produced `groups × genes` summary agreeing with the prior pickle within 1e-3

## Deferred (tracked)

- `expression summarize` `--backed` mode for reading 40+ GB atlases from disk
- `describe` / `markers` backed-atlas support
- Same pattern applied to Expression Atlas and CPTAC builders

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

## Out of scope (NOT in this PR)

- `expression summarize` gains a `--backed` mode. Needed eventually (otherwise stage-2 against a 40 GB atlas OOMs), but separable and not blocking the fused path.
- `describe` / `markers` retrofits for backed atlases.
- Extending the same streaming pattern to Expression Atlas, CPTAC, or other wide-TSV sources.
- Promoting `var` / `std` / `log2_mean` to first-class aggregation layers.
