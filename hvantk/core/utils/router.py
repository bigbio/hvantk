"""
Backend routing and reader factory.

The BackendRouter selects the best available backend for an algorithm
based on data size heuristics and the algorithm's declared capabilities.
"""

import logging
import os
from typing import List, Optional

from hvantk.core.models.backends import AlgorithmMeta, Backend
from hvantk.core.utils.readers import (
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
            from hvantk.core.utils.readers import DuckDBReader

            return DuckDBReader()
        raise ValueError(f"Unknown backend: {backend}")
