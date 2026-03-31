# Backend Abstraction Layer — Design Spec

**Date:** 2026-03-31
**Status:** Draft
**Scope:** Toolkit-wide (hvantk/core), proving ground in qtlcascade

## Problem

The hvantk toolkit persists data as Hail Tables (parquet on disk) and uses Hail/Spark as the default compute engine. However, many algorithms — colocalization ABF, gene-level aggregations, plotting — don't need Spark. They currently work around this with ad-hoc `.to_pandas()` calls scattered inside algorithm functions, coupling I/O concerns with compute logic.

This creates three problems:

1. **Unnecessary Spark dependency** — algorithms like coloc ABF are pure NumPy/pandas but require a running Spark context because they read data through Hail.
2. **No backend flexibility** — there's no way to run a pipeline stage on DuckDB or pandas even when the data fits in memory and Spark adds overhead.
3. **Implicit boundaries** — the transition between Hail and local compute is buried inside functions rather than explicit at the architecture level.

## Goals

- Decouple I/O (reading/writing Hail Tables) from compute (algorithm execution).
- Make local algorithms (pandas, DuckDB, NumPy) first-class citizens, not escape hatches.
- Enable reading `.ht` files without starting Spark, using pyarrow or DuckDB directly.
- Automatically select the best backend based on data size and algorithm support.
- Maintain full backward compatibility with existing Builder/Streamer/Downloader code.
- Persist all outputs as Hail Tables (the common format for downstream stages).

## Non-Goals

- Replacing Hail for large-scale distributed computation.
- Supporting non-parquet persistence formats.
- User-facing backend selection (the system decides).
- Migrating all existing builders/streamers in one pass.

## Architecture

### Core Abstractions

Four new components in `hvantk/core/`:

#### 1. Backend Enum (`backends.py`)

```python
from enum import Enum

class Backend(Enum):
    HAIL = "hail"
    PANDAS = "pandas"
    DUCKDB = "duckdb"
```

#### 2. Algorithm Decorator (`backends.py`)

Declares which backends a function supports and what data format it expects/produces.

```python
from dataclasses import dataclass
from typing import List, Callable

@dataclass
class AlgorithmMeta:
    backends: List[Backend]
    input_format: str      # "dataframe" or "table"
    output_format: str     # "dataframe" or "table"
    key_fields: List[str] = None  # for writer to key the output Hail Table

def algorithm(
    backends: List[Backend],
    input_format: str = "dataframe",
    output_format: str = "dataframe",
    key_fields: List[str] = None,
) -> Callable:
    """Decorator that attaches backend metadata to a function."""
    def decorator(fn):
        fn._algorithm_meta = AlgorithmMeta(
            backends=backends,
            input_format=input_format,
            output_format=output_format,
            key_fields=key_fields,
        )
        return fn
    return decorator

def get_algorithm_meta(fn) -> AlgorithmMeta:
    """Retrieve metadata from a decorated function.

    Functions without @algorithm are implicitly Hail-only.
    """
    return getattr(fn, "_algorithm_meta", AlgorithmMeta(
        backends=[Backend.HAIL],
        input_format="table",
        output_format="table",
    ))
```

#### 3. DataReader Protocol (`readers.py`)

One implementation per backend. All read Hail Table directories (parquet on disk) and deliver data in the consumer's preferred format.

```python
from typing import Protocol, Optional, Any
import pandas as pd

class DataReader(Protocol):
    backend: Backend

    def read(self, path: str) -> Any:
        """Read in backend's native format."""
        ...

    def as_dataframe(self, path: str) -> pd.DataFrame:
        """Read and return as pandas DataFrame."""
        ...

    def as_hail_table(self, path: str):
        """Read and return as Hail Table."""
        ...

    def row_count_hint(self, path: str) -> Optional[int]:
        """Cheap row count estimate from parquet metadata."""
        ...
```

**HailReader:**
- `read()` → `hl.read_table(path)`
- `as_dataframe()` → `hl.read_table(path).to_pandas()`
- `as_hail_table()` → `hl.read_table(path)`
- `row_count_hint()` → reads parquet metadata via pyarrow (no Spark)

**PandasReader:**
- `read()` → `pd.read_parquet(os.path.join(path, "rows", "parts"))`
- `as_dataframe()` → same as `read()`
- `as_hail_table()` → `hl.Table.from_pandas(self.as_dataframe(path))`
- `row_count_hint()` → parquet footer metadata

**DuckDBReader:**
- `read()` → DuckDB relation from `read_parquet()` glob
- `as_dataframe()` → `.fetchdf()` on the relation
- `as_hail_table()` → via DataFrame intermediate
- `row_count_hint()` → `SELECT count(*) FROM read_parquet(...)` (fast on parquet)

Lazy imports: DuckDB and pyarrow are imported inside class methods, making them optional dependencies. The router skips backends whose imports fail.

#### 4. HailTableWriter (`writers.py`)

Always persists output as Hail Table format.

```python
class HailTableWriter:
    def write(
        self,
        data,
        path: str,
        key: List[str] = None,
        overwrite: bool = False,
    ) -> None:
        """Write data as a Hail Table.

        Args:
            data: pd.DataFrame or hl.Table
            path: output .ht path
            key: fields to key the table by
            overwrite: overwrite existing table
        """
        if isinstance(data, pd.DataFrame):
            import hail as hl
            ht = hl.Table.from_pandas(data)
            if key:
                ht = ht.key_by(*key)
            ht.write(path, overwrite=overwrite)
        else:
            # Already a Hail Table
            if key:
                data = data.key_by(*key)
            data.write(path, overwrite=overwrite)
```

Future optimization: write parquet directly via pyarrow in Hail Table directory layout, avoiding Spark entirely for the write path. Not in initial scope.

### Backend Router (`router.py`)

Selects the best backend at runtime based on data size and algorithm capabilities.

```python
class BackendRouter:
    DEFAULT_THRESHOLDS = {
        "small": 500_000,     # rows — prefer pandas
        "large": 5_000_000,   # rows — prefer Hail
    }

    def __init__(self, thresholds: dict = None):
        self.thresholds = thresholds or self.DEFAULT_THRESHOLDS

    def resolve(
        self,
        algorithm_meta: AlgorithmMeta,
        data_paths: List[str],
    ) -> Backend:
        available = self._filter_available(algorithm_meta.backends)
        row_hint = self._estimate_size(data_paths)

        if row_hint is not None and row_hint < self.thresholds["small"]:
            preference = [Backend.PANDAS, Backend.DUCKDB, Backend.HAIL]
        elif row_hint is not None and row_hint < self.thresholds["large"]:
            preference = [Backend.DUCKDB, Backend.PANDAS, Backend.HAIL]
        else:
            preference = [Backend.HAIL, Backend.DUCKDB, Backend.PANDAS]

        return next(b for b in preference if b in available)

    def _filter_available(self, backends: List[Backend]) -> List[Backend]:
        """Filter to backends whose dependencies are importable."""
        available = []
        for b in backends:
            if b == Backend.HAIL:
                try:
                    import hail
                    available.append(b)
                except ImportError:
                    pass
            elif b == Backend.DUCKDB:
                try:
                    import duckdb
                    available.append(b)
                except ImportError:
                    pass
            elif b == Backend.PANDAS:
                available.append(b)  # pandas is always available
        return available

    def _estimate_size(self, data_paths: List[str]) -> Optional[int]:
        """Estimate total row count from parquet metadata."""
        try:
            import pyarrow.parquet as pq
            total = 0
            for path in data_paths:
                parts_dir = os.path.join(path, "rows", "parts")
                dataset = pq.ParquetDataset(parts_dir)
                for fragment in dataset.fragments:
                    total += fragment.metadata.num_rows
            return total
        except Exception:
            return None  # Unknown size → fall back to Hail preference
```

### Reader Factory (`router.py`)

Creates the appropriate reader for a resolved backend:

```python
class ReaderFactory:
    def create(self, backend: Backend) -> DataReader:
        if backend == Backend.HAIL:
            return HailReader()
        elif backend == Backend.PANDAS:
            return PandasReader()
        elif backend == Backend.DUCKDB:
            return DuckDBReader()
        raise ValueError(f"Unknown backend: {backend}")
```

### Pipeline Integration

The pipeline orchestrator wires reader → algorithm → writer:

```python
class PipelineStep:
    def __init__(self, router: BackendRouter, reader_factory: ReaderFactory,
                 writer: HailTableWriter):
        self.router = router
        self.reader_factory = reader_factory
        self.writer = writer

    def run(self, step_fn, input_paths: List[str], output_path: str,
            key: List[str] = None, **kwargs):
        meta = get_algorithm_meta(step_fn)
        backend = self.router.resolve(meta, input_paths)
        reader = self.reader_factory.create(backend)

        # Read in the format the algorithm expects
        if meta.input_format == "dataframe":
            inputs = [reader.as_dataframe(p) for p in input_paths]
        else:
            inputs = [reader.as_hail_table(p) for p in input_paths]

        result = step_fn(*inputs, **kwargs)

        # Persist as Hail Table
        output_key = key or meta.key_fields
        self.writer.write(result, output_path, key=output_key)
```

## Backward Compatibility

- **Existing Builder/Streamer/Downloader code** — unchanged. Functions without `@algorithm` are implicitly `Backend.HAIL` only.
- **Existing CLI commands** — unchanged. They call pipeline stages, unaware of backend.
- **Existing test fixtures** — `hail_session` fixture continues to work for Hail-backend tests.
- **Migration is incremental** — add `@algorithm` to functions one at a time, starting with qtlcascade.

## Reading Hail Tables Without Hail

Hail Table on-disk layout:
```
my_table.ht/
├── rows/
│   └── parts/
│       ├── part-0-*.parquet
│       ├── part-1-*.parquet
│       └── ...
├── metadata.json.gz
└── _SUCCESS
```

The parquet files under `rows/parts/` are standard parquet and can be read by pyarrow or DuckDB directly.

**Hail-specific type handling:**
- `hl.Locus` → parquet struct `{contig: string, position: int32}` → flattened to columns or kept as struct depending on algorithm needs
- `hl.tarray(hl.tstr)` (alleles) → parquet list of strings → pandas list column
- Nested Hail structs → pyarrow nested columns → may need flattening per-algorithm

Type handling logic lives in the readers, not in algorithms.

## Application to qtlcascade

| Stage | Function | Current Engine | Declared Backends | Input Format |
|-------|----------|---------------|-------------------|-------------|
| Build cascade | `build_cascade()` | Hail | HAIL, DUCKDB | table |
| Coloc ABF | `run_coloc_per_gene()` | Hail→pandas→NumPy | PANDAS | dataframe |
| Gene summary | `build_cascade_gene_summary()` | Hail + pandas overlay | HAIL, DUCKDB | table |
| Plot | `plot_*()` | matplotlib | PANDAS | dataframe |
| Report | `generate_report()` | pandas/jinja | PANDAS | dataframe |

**Key refactor — coloc.py:**
- Current: `run_coloc_per_gene(eqtl_path, pqtl_path, ...)` reads Hail Tables internally, calls `.to_pandas()`, runs ABF
- New: `run_coloc_per_gene(joined_df: pd.DataFrame, ...)` receives a DataFrame from the reader. The Hail→pandas boundary moves to the reader layer.
- `compute_log_abf()` and `coloc_abf()` — pure NumPy, no changes needed.

## File Layout

New files:
```
hvantk/core/
├── backends.py    # Backend enum, AlgorithmMeta, @algorithm decorator, get_algorithm_meta()
├── readers.py     # DataReader protocol, HailReader, PandasReader, DuckDBReader
├── writers.py     # HailTableWriter
├── router.py      # BackendRouter, ReaderFactory
```

Modified files:
```
hvantk/core/__init__.py           # Export new public API
hvantk/qtlcascade/coloc.py        # Refactor to receive DataFrames, add @algorithm
hvantk/qtlcascade/cascade.py      # Add @algorithm decorator
hvantk/qtlcascade/gene_summary.py # Add @algorithm decorator
hvantk/qtlcascade/pipeline.py     # Use BackendRouter + readers
```

## Dependencies

- **pyarrow** — already a transitive dependency via Hail. Used by PandasReader for parquet I/O and by all readers for `row_count_hint()`.
- **duckdb** — new optional dependency. DuckDBReader and DuckDB backend are skipped if not installed.
- **pandas, numpy** — already required dependencies.

## Testing Strategy

- **Unit tests per reader** — verify each reader can read a small `.ht` directory and return correct data
- **Unit tests for router** — verify backend selection logic with mocked row counts
- **Integration tests for qtlcascade** — run coloc with pandas backend on test data, compare results to current Hail-based run
- **Backward compatibility tests** — existing builders without `@algorithm` continue to work unchanged
- **Hail type handling tests** — verify Locus/alleles structs are correctly read by pyarrow and DuckDB

## Future Extensions

- **Richer size heuristics** — metadata-driven backend hints (per-algorithm preferred size ranges) instead of global thresholds.
- **Direct parquet write** — HailTableWriter writes Hail Table directory layout via pyarrow, bypassing Spark entirely for local outputs.
- **Polars backend** — new `Backend.POLARS` entry, PolarReader implementation. The abstraction supports it without changes to algorithms.
- **Per-algorithm backend overrides** — allow algorithms to influence preference order beyond just declaring support.
