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

from hvantk.core.models.backends import Backend

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

    @staticmethod
    def _safe_glob(path: str) -> str:
        """Build a parquet glob path with escaped single quotes."""
        parts_glob = os.path.join(_parts_path(path), "*.parquet")
        return parts_glob.replace("'", "''")

    def read(self, path: str):
        safe = self._safe_glob(path)
        return self._con.sql(f"SELECT * FROM read_parquet('{safe}')")

    def as_dataframe(self, path: str) -> pd.DataFrame:
        return self.read(path).fetchdf()

    def as_hail_table(self, path: str):
        import hail as hl

        return hl.Table.from_pandas(self.as_dataframe(path))

    def row_count_hint(self, path: str) -> Optional[int]:
        try:
            safe = self._safe_glob(path)
            result = self._con.sql(
                f"SELECT count(*) AS n FROM read_parquet('{safe}')"
            ).fetchone()
            return result[0] if result else None
        except Exception:
            return None
