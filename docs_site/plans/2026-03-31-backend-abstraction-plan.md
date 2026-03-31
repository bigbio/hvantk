# Backend Abstraction Layer — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Decouple I/O from compute across the hvantk toolkit so algorithms can run on Hail, pandas, or DuckDB backends — with Hail Table parquet format as the universal persistence layer.

**Architecture:** A declarative `@algorithm` decorator marks each function's supported backends and expected I/O format. Backend-specific `DataReader` implementations read `.ht` parquet files and deliver data in the algorithm's preferred format. A `BackendRouter` selects the best backend at runtime based on data size heuristics. All outputs persist as Hail Tables.

**Tech Stack:** Python 3.10+, Hail, pandas, pyarrow (parquet I/O), DuckDB (optional), NumPy

**Design Spec:** `docs_site/specs/2026-03-31-backend-abstraction-design.md`

---

## File Structure

### New files

| File | Responsibility |
|------|----------------|
| `hvantk/core/backends.py` | `Backend` enum, `AlgorithmMeta` dataclass, `@algorithm` decorator, `get_algorithm_meta()` |
| `hvantk/core/readers.py` | `DataReader` protocol, `HailReader`, `PandasReader`, `DuckDBReader` implementations |
| `hvantk/core/writers.py` | `HailTableWriter` — converts DataFrame or Hail Table → persisted `.ht` |
| `hvantk/core/router.py` | `BackendRouter` (size-based backend selection), `ReaderFactory` |

### Modified files

| File | Changes |
|------|---------|
| `hvantk/core/__init__.py` | Export `Backend`, `algorithm`, `get_algorithm_meta`, `BackendRouter`, `ReaderFactory`, `HailTableWriter` |
| `hvantk/qtlcascade/coloc.py` | Add `@algorithm` decorator; extract Hail I/O from `run_coloc_per_gene` so it receives a DataFrame |
| `hvantk/qtlcascade/cascade.py` | Add `@algorithm` decorator to `build_cascade` |
| `hvantk/qtlcascade/gene_summary.py` | Add `@algorithm` decorator to `build_cascade_gene_summary` |
| `hvantk/qtlcascade/pipeline.py` | Integrate `BackendRouter` + readers into stage execution |
| `hvantk/qtlcascade/__init__.py` | Re-export new public symbols if needed |
| `pyproject.toml` | Add `duckdb` as optional dependency |

---

## Task 1: Backend Enum and Algorithm Decorator

**Files:**
- Create: `hvantk/core/backends.py`

This is the foundation — the decorator and metadata that everything else depends on.

- [ ] **Step 1: Create `hvantk/core/backends.py` with `Backend` enum and `AlgorithmMeta`**

```python
"""
Backend declarations and algorithm metadata.

The @algorithm decorator marks functions with the backends they support
and the data format they expect, enabling the BackendRouter to select
the best execution strategy at runtime.
"""

from dataclasses import dataclass, field
from enum import Enum
from functools import wraps
from typing import Callable, List, Optional


class Backend(Enum):
    """Supported compute/IO backends."""

    HAIL = "hail"
    PANDAS = "pandas"
    DUCKDB = "duckdb"


@dataclass
class AlgorithmMeta:
    """Metadata attached to algorithm functions by the @algorithm decorator.

    Attributes
    ----------
    backends : list[Backend]
        Backends the algorithm supports (hard constraint for the router).
    input_format : str
        Expected input type: ``"dataframe"`` (pandas) or ``"table"`` (Hail).
    output_format : str
        Return type: ``"dataframe"`` or ``"table"``.
    key_fields : list[str] or None
        Fields to key the output Hail Table by when persisting.
    """

    backends: List[Backend]
    input_format: str = "dataframe"
    output_format: str = "dataframe"
    key_fields: Optional[List[str]] = None


_IMPLICIT_HAIL_META = AlgorithmMeta(
    backends=[Backend.HAIL],
    input_format="table",
    output_format="table",
)


def algorithm(
    backends: List[Backend],
    input_format: str = "dataframe",
    output_format: str = "dataframe",
    key_fields: Optional[List[str]] = None,
) -> Callable:
    """Declare an algorithm's supported backends and I/O formats.

    This decorator attaches an :class:`AlgorithmMeta` instance to the
    function as ``_algorithm_meta``.  It does not alter the function's
    runtime behaviour.

    Parameters
    ----------
    backends : list[Backend]
        Which backends can execute this algorithm.
    input_format : str
        ``"dataframe"`` or ``"table"`` — what the function receives.
    output_format : str
        ``"dataframe"`` or ``"table"`` — what the function returns.
    key_fields : list[str], optional
        Hail Table key fields for persisting the output.
    """

    def decorator(fn: Callable) -> Callable:
        fn._algorithm_meta = AlgorithmMeta(
            backends=backends,
            input_format=input_format,
            output_format=output_format,
            key_fields=key_fields,
        )

        @wraps(fn)
        def wrapper(*args, **kwargs):
            return fn(*args, **kwargs)

        wrapper._algorithm_meta = fn._algorithm_meta
        return wrapper

    return decorator


def get_algorithm_meta(fn: Callable) -> AlgorithmMeta:
    """Retrieve backend metadata from a decorated function.

    Functions without ``@algorithm`` are treated as Hail-only.
    """
    return getattr(fn, "_algorithm_meta", _IMPLICIT_HAIL_META)
```

- [ ] **Step 2: Verify the module imports cleanly**

Run: `cd /Users/enrique/projects/github/pyvatk && python -c "from hvantk.core.backends import Backend, algorithm, get_algorithm_meta, AlgorithmMeta; print('OK')"`

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add hvantk/core/backends.py
git commit -m "feat(core): add Backend enum and @algorithm decorator"
```

---

## Task 2: DataReader Protocol and Implementations

**Files:**
- Create: `hvantk/core/readers.py`

Three reader implementations that read Hail Table parquet files and deliver data in the consumer's preferred format.

- [ ] **Step 1: Create `hvantk/core/readers.py`**

```python
"""
Backend-specific data readers for Hail Table parquet files.

Hail Tables on disk are directories containing partitioned parquet files
under ``rows/parts/``.  These readers allow reading that data via Hail
(requires Spark), pandas/pyarrow (local), or DuckDB (local).
"""

import logging
import os
from typing import Optional, Protocol, runtime_checkable

import pandas as pd

from hvantk.core.backends import Backend

logger = logging.getLogger(__name__)


@runtime_checkable
class DataReader(Protocol):
    """Protocol for reading Hail Table directories."""

    backend: Backend

    def read(self, path: str):
        """Read in the backend's native format."""
        ...

    def as_dataframe(self, path: str) -> pd.DataFrame:
        """Read and return as a pandas DataFrame."""
        ...

    def as_hail_table(self, path: str):
        """Read and return as a Hail Table."""
        ...

    def row_count_hint(self, path: str) -> Optional[int]:
        """Cheap row-count estimate from parquet metadata (no data loading)."""
        ...


def _parts_path(ht_path: str) -> str:
    """Return the ``rows/parts`` directory inside a ``.ht`` directory."""
    return os.path.join(ht_path, "rows", "parts")


def _parquet_row_count(ht_path: str) -> Optional[int]:
    """Estimate row count from parquet footer metadata."""
    try:
        import pyarrow.parquet as pq

        parts_dir = _parts_path(ht_path)
        dataset = pq.ParquetDataset(parts_dir)
        total = 0
        for fragment in dataset.fragments:
            total += fragment.metadata.num_rows
        return total
    except Exception:
        return None


class HailReader:
    """Read Hail Tables via Hail (requires a running Spark context)."""

    backend = Backend.HAIL

    def read(self, path: str):
        import hail as hl

        return hl.read_table(path)

    def as_dataframe(self, path: str) -> pd.DataFrame:
        return self.read(path).to_pandas()

    def as_hail_table(self, path: str):
        return self.read(path)

    def row_count_hint(self, path: str) -> Optional[int]:
        return _parquet_row_count(path)


class PandasReader:
    """Read Hail Table parquet files directly via pyarrow — no Spark needed."""

    backend = Backend.PANDAS

    def read(self, path: str) -> pd.DataFrame:
        parts_dir = _parts_path(path)
        return pd.read_parquet(parts_dir)

    def as_dataframe(self, path: str) -> pd.DataFrame:
        return self.read(path)

    def as_hail_table(self, path: str):
        import hail as hl

        return hl.Table.from_pandas(self.as_dataframe(path))

    def row_count_hint(self, path: str) -> Optional[int]:
        return _parquet_row_count(path)


class DuckDBReader:
    """Read Hail Table parquet files via DuckDB — no Spark needed."""

    backend = Backend.DUCKDB

    def __init__(self):
        import duckdb

        self._con = duckdb.connect()

    def read(self, path: str):
        parts_glob = os.path.join(_parts_path(path), "*.parquet")
        return self._con.sql(
            f"SELECT * FROM read_parquet('{parts_glob}')"
        )

    def as_dataframe(self, path: str) -> pd.DataFrame:
        return self.read(path).fetchdf()

    def as_hail_table(self, path: str):
        import hail as hl

        return hl.Table.from_pandas(self.as_dataframe(path))

    def row_count_hint(self, path: str) -> Optional[int]:
        try:
            parts_glob = os.path.join(_parts_path(path), "*.parquet")
            result = self._con.sql(
                f"SELECT count(*) AS n FROM read_parquet('{parts_glob}')"
            ).fetchone()
            return result[0] if result else None
        except Exception:
            return None
```

- [ ] **Step 2: Verify the module imports cleanly (pandas reader — no Spark)**

Run: `cd /Users/enrique/projects/github/pyvatk && python -c "from hvantk.core.readers import PandasReader, HailReader, DuckDBReader; print('OK')"`

Expected: `OK` (DuckDB import may fail if not installed — that's fine, the import of the module itself should succeed since duckdb is imported lazily inside `__init__`)

- [ ] **Step 3: Commit**

```bash
git add hvantk/core/readers.py
git commit -m "feat(core): add DataReader protocol with Hail, pandas, and DuckDB backends"
```

---

## Task 3: HailTableWriter

**Files:**
- Create: `hvantk/core/writers.py`

- [ ] **Step 1: Create `hvantk/core/writers.py`**

```python
"""
Hail Table writer — universal persistence layer.

All algorithm outputs are written as Hail Tables regardless of the
compute backend used, ensuring downstream stages can always consume
them via any reader.
"""

import logging
from typing import List, Optional, Union

import pandas as pd

logger = logging.getLogger(__name__)


class HailTableWriter:
    """Write data as a Hail Table.

    Accepts either a pandas DataFrame or a Hail Table.  DataFrames are
    converted via ``hl.Table.from_pandas()`` before writing.
    """

    def write(
        self,
        data: Union[pd.DataFrame, "hl.Table"],
        path: str,
        key: Optional[List[str]] = None,
        overwrite: bool = False,
    ) -> None:
        """Persist data as a Hail Table.

        Parameters
        ----------
        data : pd.DataFrame or hl.Table
            The data to write.
        path : str
            Output ``.ht`` directory path.
        key : list[str], optional
            Fields to key the Hail Table by.
        overwrite : bool
            Overwrite an existing table at *path*.
        """
        import hail as hl

        if isinstance(data, pd.DataFrame):
            logger.info("Converting DataFrame (%d rows) to Hail Table", len(data))
            ht = hl.Table.from_pandas(data)
        else:
            ht = data

        if key:
            ht = ht.key_by(*key)

        logger.info("Writing Hail Table to %s", path)
        ht.write(path, overwrite=overwrite)
```

- [ ] **Step 2: Verify the module imports cleanly**

Run: `cd /Users/enrique/projects/github/pyvatk && python -c "from hvantk.core.writers import HailTableWriter; print('OK')"`

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add hvantk/core/writers.py
git commit -m "feat(core): add HailTableWriter for universal persistence"
```

---

## Task 4: BackendRouter and ReaderFactory

**Files:**
- Create: `hvantk/core/router.py`

- [ ] **Step 1: Create `hvantk/core/router.py`**

```python
"""
Backend routing and reader factory.

The BackendRouter selects the best available backend for an algorithm
based on data size heuristics and the algorithm's declared capabilities.
"""

import logging
import os
from typing import List, Optional

from hvantk.core.backends import AlgorithmMeta, Backend
from hvantk.core.readers import (
    DataReader,
    HailReader,
    PandasReader,
    _parquet_row_count,
)

logger = logging.getLogger(__name__)

# Default row-count thresholds for backend preference
_SMALL_THRESHOLD = 500_000
_LARGE_THRESHOLD = 5_000_000


class BackendRouter:
    """Select the best backend for an algorithm given data size.

    Parameters
    ----------
    small_threshold : int
        Below this row count, prefer local backends (pandas > DuckDB).
    large_threshold : int
        Above this row count, prefer distributed backends (Hail > DuckDB).
    """

    def __init__(
        self,
        small_threshold: int = _SMALL_THRESHOLD,
        large_threshold: int = _LARGE_THRESHOLD,
    ):
        self.small_threshold = small_threshold
        self.large_threshold = large_threshold

    def resolve(
        self,
        meta: AlgorithmMeta,
        data_paths: Optional[List[str]] = None,
    ) -> Backend:
        """Pick the best backend from *meta.backends*.

        Parameters
        ----------
        meta : AlgorithmMeta
            Algorithm's declared backend support.
        data_paths : list[str], optional
            Paths to input ``.ht`` directories for size estimation.

        Returns
        -------
        Backend
            The selected backend.

        Raises
        ------
        RuntimeError
            If none of the algorithm's declared backends are available.
        """
        available = self._filter_available(meta.backends)
        if not available:
            raise RuntimeError(
                f"No available backends among {meta.backends}. "
                "Install missing dependencies (duckdb, hail)."
            )

        row_hint = self._estimate_size(data_paths) if data_paths else None

        if row_hint is not None and row_hint < self.small_threshold:
            preference = [Backend.PANDAS, Backend.DUCKDB, Backend.HAIL]
        elif row_hint is not None and row_hint < self.large_threshold:
            preference = [Backend.DUCKDB, Backend.PANDAS, Backend.HAIL]
        else:
            # Large data or unknown size → prefer Hail
            preference = [Backend.HAIL, Backend.DUCKDB, Backend.PANDAS]

        selected = next(b for b in preference if b in available)
        logger.info(
            "BackendRouter: selected %s (row_hint=%s, available=%s)",
            selected.value,
            row_hint,
            [b.value for b in available],
        )
        return selected

    @staticmethod
    def _filter_available(backends: List[Backend]) -> List[Backend]:
        """Return backends whose dependencies are importable."""
        available = []
        for b in backends:
            if b == Backend.PANDAS:
                available.append(b)  # pandas is always available
            elif b == Backend.HAIL:
                try:
                    import hail  # noqa: F401

                    available.append(b)
                except ImportError:
                    pass
            elif b == Backend.DUCKDB:
                try:
                    import duckdb  # noqa: F401

                    available.append(b)
                except ImportError:
                    pass
        return available

    @staticmethod
    def _estimate_size(data_paths: List[str]) -> Optional[int]:
        """Estimate total row count from parquet metadata."""
        total = 0
        for path in data_paths:
            hint = _parquet_row_count(path)
            if hint is None:
                return None  # Can't estimate → unknown
            total += hint
        return total


class ReaderFactory:
    """Create a DataReader for a given backend."""

    @staticmethod
    def create(backend: Backend) -> DataReader:
        if backend == Backend.HAIL:
            return HailReader()
        elif backend == Backend.PANDAS:
            return PandasReader()
        elif backend == Backend.DUCKDB:
            from hvantk.core.readers import DuckDBReader

            return DuckDBReader()
        raise ValueError(f"Unknown backend: {backend}")
```

- [ ] **Step 2: Verify the module imports cleanly**

Run: `cd /Users/enrique/projects/github/pyvatk && python -c "from hvantk.core.router import BackendRouter, ReaderFactory; print('OK')"`

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add hvantk/core/router.py
git commit -m "feat(core): add BackendRouter and ReaderFactory"
```

---

## Task 5: Export Public API from `hvantk/core`

**Files:**
- Modify: `hvantk/core/__init__.py`

- [ ] **Step 1: Update `hvantk/core/__init__.py` to export the new public API**

Replace the entire file content with:

```python
# hvantk.core package
from hvantk.core.backends import (
    AlgorithmMeta,
    Backend,
    algorithm,
    get_algorithm_meta,
)
from hvantk.core.readers import (
    DataReader,
    DuckDBReader,
    HailReader,
    PandasReader,
)
from hvantk.core.router import BackendRouter, ReaderFactory
from hvantk.core.writers import HailTableWriter

__all__ = [
    "AlgorithmMeta",
    "Backend",
    "BackendRouter",
    "DataReader",
    "DuckDBReader",
    "HailReader",
    "HailTableWriter",
    "PandasReader",
    "ReaderFactory",
    "algorithm",
    "get_algorithm_meta",
]
```

- [ ] **Step 2: Verify the public API imports**

Run: `cd /Users/enrique/projects/github/pyvatk && python -c "from hvantk.core import Backend, algorithm, BackendRouter, ReaderFactory, HailTableWriter; print('OK')"`

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add hvantk/core/__init__.py
git commit -m "feat(core): export backend abstraction public API"
```

---

## Task 6: Add DuckDB as Optional Dependency

**Files:**
- Modify: `pyproject.toml`

- [ ] **Step 1: Add `duckdb` to optional dependencies in `pyproject.toml`**

In the `[tool.poetry.dependencies]` section, add after the `plotly` line:

```toml
duckdb = { version = ">=0.9.0", optional = true }
```

In the `[tool.poetry.extras]` section, add:

```toml
duckdb = ["duckdb"]
```

- [ ] **Step 2: Verify pyproject.toml is valid**

Run: `cd /Users/enrique/projects/github/pyvatk && python -c "import tomllib; tomllib.load(open('pyproject.toml','rb')); print('OK')"`

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add pyproject.toml
git commit -m "build: add duckdb as optional dependency"
```

---

## Task 7: Decorate `coloc.py` Algorithms and Refactor `run_coloc_per_gene`

**Files:**
- Modify: `hvantk/qtlcascade/coloc.py`

This is the most impactful refactor. The goal is to:
1. Add `@algorithm` to `compute_log_abf` and `coloc_abf` (pure NumPy, pandas backend)
2. Split `run_coloc_per_gene` into two parts:
   - A Hail-based data preparation function (joins allpairs tables, exports to DataFrame)
   - A pure pandas/numpy coloc driver that receives a DataFrame

- [ ] **Step 1: Add the `@algorithm` import to `coloc.py`**

At `hvantk/qtlcascade/coloc.py`, add to the imports (after line 29, after `import pandas as pd`):

```python
from hvantk.core.backends import Backend, algorithm
```

- [ ] **Step 2: Decorate `compute_log_abf` (line 47)**

Add the decorator before the function definition at line 47:

```python
@algorithm(
    backends=[Backend.PANDAS],
    input_format="dataframe",
    output_format="dataframe",
)
def compute_log_abf(
```

- [ ] **Step 3: Decorate `coloc_abf` (line 87)**

Add the decorator before the function definition at line 87:

```python
@algorithm(
    backends=[Backend.PANDAS],
    input_format="dataframe",
    output_format="dataframe",
)
def coloc_abf(
```

- [ ] **Step 4: Create `prepare_coloc_data` — the Hail I/O extraction function**

Add a new function before `run_coloc_per_gene` (before line 166). This extracts the Hail-specific data loading logic that currently lives inside `run_coloc_per_gene` (lines 202–267):

```python
@algorithm(
    backends=[Backend.HAIL],
    input_format="table",
    output_format="dataframe",
)
def prepare_coloc_data(
    eqtl_allpairs_ht_path: str,
    pqtl_allpairs_ht_path: str,
    cascade_genes: list,
    tissue: Optional[str] = None,
) -> pd.DataFrame:
    """Load and join allpairs eQTL/pQTL data for coloc.

    Reads both allpairs Hail Tables, filters to cascade genes,
    inner-joins on ``(locus, alleles, gene_id)``, and flattens to a
    pandas DataFrame.

    Parameters
    ----------
    eqtl_allpairs_ht_path : str
        Path to allpairs eQTL Hail Table.
    pqtl_allpairs_ht_path : str
        Path to allpairs pQTL Hail Table.
    cascade_genes : list[str]
        Gene IDs with both eQTL and pQTL evidence.
    tissue : str, optional
        Filter allpairs tables to this tissue.

    Returns
    -------
    pd.DataFrame
        Columns: gene_id, pos, eqtl_beta, eqtl_se, eqtl_p, pqtl_beta,
        pqtl_se.
    """
    import hail as hl

    eqtl_ht = hl.read_table(eqtl_allpairs_ht_path)
    pqtl_ht = hl.read_table(pqtl_allpairs_ht_path)

    # Filter to cascade genes (single Spark filter)
    gene_set = hl.literal(set(cascade_genes))
    eqtl_ht = eqtl_ht.filter(gene_set.contains(eqtl_ht.gene_id))
    pqtl_ht = pqtl_ht.filter(gene_set.contains(pqtl_ht.gene_id))

    eqtl_has_tissue = "tissue" in list(eqtl_ht.row)
    pqtl_has_tissue = "tissue" in list(pqtl_ht.row)

    if tissue:
        if eqtl_has_tissue:
            eqtl_ht = eqtl_ht.filter(eqtl_ht.tissue == tissue)
        if pqtl_has_tissue:
            pqtl_ht = pqtl_ht.filter(pqtl_ht.tissue == tissue)
    elif eqtl_has_tissue or pqtl_has_tissue:
        logger.warning(
            "No tissue filter provided but allpairs table(s) contain a "
            "'tissue' field. Inner join on (locus, alleles, gene_id) may "
            "cross-multiply rows from different tissues."
        )

    # Select fields for join; annotate position for windowing
    eqtl_sel = eqtl_ht.select(
        eqtl_beta=eqtl_ht.beta,
        eqtl_se=eqtl_ht.se,
        eqtl_p=eqtl_ht.p_value,
        position=eqtl_ht.locus.position,
    )
    pqtl_sel = pqtl_ht.select(
        pqtl_beta=pqtl_ht.beta,
        pqtl_se=pqtl_ht.se,
    )

    # Inner join on (locus, alleles, gene_id) — single Spark job
    joined = eqtl_sel.join(pqtl_sel, how="inner")

    # Flatten for pandas export
    joined = joined.key_by()
    joined = joined.annotate(
        contig=joined.locus.contig,
        pos=joined.locus.position,
    )
    joined = joined.select(
        "gene_id",
        "pos",
        "eqtl_beta",
        "eqtl_se",
        "eqtl_p",
        "pqtl_beta",
        "pqtl_se",
    )

    logger.info(
        "Exporting joined allpairs for coloc (%d cascade genes)", len(cascade_genes)
    )
    return joined.to_pandas()
```

- [ ] **Step 5: Refactor `run_coloc_per_gene` to use `prepare_coloc_data`**

Replace the entire `run_coloc_per_gene` function (lines 166–305) with:

```python
@algorithm(
    backends=[Backend.PANDAS],
    input_format="dataframe",
    output_format="dataframe",
    key_fields=["gene_id"],
)
def run_coloc_per_gene(
    eqtl_allpairs_ht_path: str,
    pqtl_allpairs_ht_path: str,
    cascade_genes: list,
    tissue: Optional[str] = None,
    window_kb: int = DEFAULT_COLOC_WINDOW_KB,
    p1: float = DEFAULT_COLOC_P1,
    p2: float = DEFAULT_COLOC_P2,
    p12: float = DEFAULT_COLOC_P12,
    W: float = DEFAULT_COLOC_W,
) -> pd.DataFrame:
    """Run coloc for all cascade genes.

    Uses ``prepare_coloc_data`` for bulk data extraction and NumPy for
    per-gene ABF computation.

    Parameters
    ----------
    eqtl_allpairs_ht_path : str
        Path to allpairs eQTL Hail Table.
    pqtl_allpairs_ht_path : str
        Path to allpairs pQTL Hail Table.
    cascade_genes : list[str]
        Gene IDs with both eQTL and pQTL evidence.
    tissue : str, optional
        Filter allpairs tables to this tissue.
    window_kb : int
        Window (±kb) around lead variant for regional extraction.
    p1, p2, p12, W : float
        Coloc prior parameters.

    Returns
    -------
    pd.DataFrame
        Columns: gene_id, tissue, H0–H4, n_variants.
    """
    result_cols = ["gene_id", "tissue", "H0", "H1", "H2", "H3", "H4", "n_variants"]
    empty = pd.DataFrame(columns=result_cols)

    if not cascade_genes:
        return empty

    df = prepare_coloc_data(
        eqtl_allpairs_ht_path=eqtl_allpairs_ht_path,
        pqtl_allpairs_ht_path=pqtl_allpairs_ht_path,
        cascade_genes=cascade_genes,
        tissue=tissue,
    )

    if df.empty:
        logger.warning("No overlapping variants found between allpairs tables")
        return empty

    # Per-gene coloc with regional windowing (pure Python)
    window_bp = window_kb * 1000
    results = []

    for gene_id, group in df.groupby("gene_id"):
        # Window around lead eQTL variant
        lead_pos = group.loc[group["eqtl_p"].idxmin(), "pos"]
        region = group[np.abs(group["pos"] - lead_pos) <= window_bp]

        if len(region) < 2:
            continue

        row = coloc_abf(
            eqtl_beta=region["eqtl_beta"].values,
            eqtl_se=region["eqtl_se"].values,
            pqtl_beta=region["pqtl_beta"].values,
            pqtl_se=region["pqtl_se"].values,
            p1=p1,
            p2=p2,
            p12=p12,
            W=W,
        )
        row["gene_id"] = gene_id
        row["tissue"] = tissue or "unknown"
        results.append(row)

    if not results:
        logger.warning("Coloc produced no results (check variant overlap)")
        return empty

    logger.info("Coloc completed for %d / %d genes", len(results), len(cascade_genes))
    return pd.DataFrame(results)[result_cols]
```

- [ ] **Step 6: Verify `coloc.py` imports cleanly**

Run: `cd /Users/enrique/projects/github/pyvatk && python -c "from hvantk.qtlcascade.coloc import run_coloc_per_gene, prepare_coloc_data, compute_log_abf, coloc_abf; print('OK')"`

Expected: `OK`

- [ ] **Step 7: Commit**

```bash
git add hvantk/qtlcascade/coloc.py
git commit -m "refactor(qtlcascade): decouple coloc I/O from computation with @algorithm decorators"
```

---

## Task 8: Decorate `cascade.py` and `gene_summary.py`

**Files:**
- Modify: `hvantk/qtlcascade/cascade.py`
- Modify: `hvantk/qtlcascade/gene_summary.py`

- [ ] **Step 1: Add `@algorithm` decorator to `build_cascade` in `cascade.py`**

Add the import at the top of `hvantk/qtlcascade/cascade.py` (after line 27, after the constants import):

```python
from hvantk.core.backends import Backend, algorithm
```

Add the decorator before the function definition at line 32:

```python
@algorithm(
    backends=[Backend.HAIL, Backend.DUCKDB],
    input_format="table",
    output_format="table",
    key_fields=["locus", "alleles", "gene_id"],
)
def build_cascade(
```

- [ ] **Step 2: Add `@algorithm` decorator to `build_cascade_gene_summary` in `gene_summary.py`**

Add the import at the top of `hvantk/qtlcascade/gene_summary.py` (after line 14, after the constants import):

```python
from hvantk.core.backends import Backend, algorithm
```

Add the decorator before the function definition at line 19:

```python
@algorithm(
    backends=[Backend.HAIL, Backend.DUCKDB],
    input_format="table",
    output_format="table",
    key_fields=["gene_id"],
)
def build_cascade_gene_summary(
```

- [ ] **Step 3: Verify both modules import cleanly**

Run: `cd /Users/enrique/projects/github/pyvatk && python -c "from hvantk.qtlcascade.cascade import build_cascade; from hvantk.qtlcascade.gene_summary import build_cascade_gene_summary; print('OK')"`

Expected: `OK`

- [ ] **Step 4: Commit**

```bash
git add hvantk/qtlcascade/cascade.py hvantk/qtlcascade/gene_summary.py
git commit -m "feat(qtlcascade): add @algorithm decorators to cascade and gene_summary"
```

---

## Task 9: Update `qtlcascade/__init__.py` Exports

**Files:**
- Modify: `hvantk/qtlcascade/__init__.py`

- [ ] **Step 1: Add `prepare_coloc_data` to the exports**

In `hvantk/qtlcascade/__init__.py`, update the coloc import line (around line 39) to include the new function:

Change:
```python
from hvantk.qtlcascade.coloc import coloc_abf, compute_log_abf, run_coloc_per_gene
```

To:
```python
from hvantk.qtlcascade.coloc import (
    coloc_abf,
    compute_log_abf,
    prepare_coloc_data,
    run_coloc_per_gene,
)
```

Also add `"prepare_coloc_data"` to the `__all__` list (in alphabetical position within the list).

- [ ] **Step 2: Verify the package exports**

Run: `cd /Users/enrique/projects/github/pyvatk && python -c "from hvantk.qtlcascade import prepare_coloc_data; print('OK')"`

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add hvantk/qtlcascade/__init__.py
git commit -m "feat(qtlcascade): export prepare_coloc_data"
```

---

## Task 10: Integrate BackendRouter into CascadePipeline

**Files:**
- Modify: `hvantk/qtlcascade/pipeline.py`

This wires the router into the pipeline. The pipeline currently calls algorithm functions directly with Hail Table paths. We add the router so it can resolve backends and deliver data in the right format. For this first integration, we keep the existing direct-call pattern for stages that are still Hail-only (`build_cascade`, `gene_summary`), and use the router for the coloc stage which now has explicit I/O separation.

- [ ] **Step 1: Add router imports to `pipeline.py`**

At the top of `hvantk/qtlcascade/pipeline.py`, add after the existing imports (after line 21, after `import pandas as pd`):

```python
from hvantk.core.backends import get_algorithm_meta
from hvantk.core.router import BackendRouter, ReaderFactory
```

- [ ] **Step 2: Add router to `CascadePipeline.__init__`**

In `hvantk/qtlcascade/pipeline.py`, at the end of `__init__` (after line 153), add:

```python
        self._router = BackendRouter()
        self._reader_factory = ReaderFactory()
```

- [ ] **Step 3: Update `_stage_coloc` to log the resolved backend**

In the `_stage_coloc` method (line 335), after the import of `run_coloc_per_gene` (line 337), add backend resolution logging:

```python
        meta = get_algorithm_meta(run_coloc_per_gene)
        backend = self._router.resolve(
            meta,
            [self.config.eqtl_allpairs_ht, self.config.pqtl_allpairs_ht],
        )
        logger.info("Coloc backend: %s", backend.value)
```

This doesn't change behavior yet — `run_coloc_per_gene` still calls `prepare_coloc_data` internally. But it establishes the routing pattern and logs which backend _would_ be selected, paving the way for full router-driven execution in a follow-up.

- [ ] **Step 4: Verify `pipeline.py` imports cleanly**

Run: `cd /Users/enrique/projects/github/pyvatk && python -c "from hvantk.qtlcascade.pipeline import CascadePipeline, CascadeConfig; print('OK')"`

Expected: `OK`

- [ ] **Step 5: Commit**

```bash
git add hvantk/qtlcascade/pipeline.py
git commit -m "feat(qtlcascade): integrate BackendRouter into CascadePipeline"
```

---

## Task 11: Run Existing Tests and Verify No Regressions

**Files:** None (verification only)

- [ ] **Step 1: Run the fast test suite**

Run: `cd /Users/enrique/projects/github/pyvatk && pytest -q`

Expected: All existing tests pass. The decorators don't change runtime behavior, so nothing should break.

- [ ] **Step 2: Verify decorator metadata is accessible at runtime**

Run:
```bash
cd /Users/enrique/projects/github/pyvatk && python -c "
from hvantk.core.backends import get_algorithm_meta, Backend
from hvantk.qtlcascade.coloc import run_coloc_per_gene, coloc_abf, prepare_coloc_data
from hvantk.qtlcascade.cascade import build_cascade
from hvantk.qtlcascade.gene_summary import build_cascade_gene_summary

# Check coloc functions are pandas-only
meta = get_algorithm_meta(run_coloc_per_gene)
assert meta.backends == [Backend.PANDAS], f'Expected [PANDAS], got {meta.backends}'
assert meta.input_format == 'dataframe'

meta = get_algorithm_meta(coloc_abf)
assert meta.backends == [Backend.PANDAS]

# Check prepare_coloc_data is hail-only
meta = get_algorithm_meta(prepare_coloc_data)
assert meta.backends == [Backend.HAIL]

# Check cascade has both backends
meta = get_algorithm_meta(build_cascade)
assert Backend.HAIL in meta.backends
assert Backend.DUCKDB in meta.backends

# Check gene_summary has both backends
meta = get_algorithm_meta(build_cascade_gene_summary)
assert Backend.HAIL in meta.backends
assert Backend.DUCKDB in meta.backends

print('All metadata checks passed')
"
```

Expected: `All metadata checks passed`

- [ ] **Step 3: Commit (only if tests revealed issues that needed fixing)**

If any fixes were needed, commit them:
```bash
git add -u
git commit -m "fix: address test regressions from backend abstraction"
```

---

## Summary

| Task | What it delivers |
|------|-----------------|
| 1 | `Backend` enum, `@algorithm` decorator, `AlgorithmMeta` |
| 2 | `DataReader` protocol + Hail, Pandas, DuckDB readers |
| 3 | `HailTableWriter` — universal persistence |
| 4 | `BackendRouter` + `ReaderFactory` — automatic backend selection |
| 5 | Public API exports from `hvantk/core` |
| 6 | DuckDB optional dependency in `pyproject.toml` |
| 7 | Coloc refactor — I/O decoupled from ABF computation |
| 8 | `@algorithm` on `build_cascade` and `build_cascade_gene_summary` |
| 9 | Updated `qtlcascade/__init__.py` exports |
| 10 | Router integrated into `CascadePipeline` |
| 11 | Regression verification |
