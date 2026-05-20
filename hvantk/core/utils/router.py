"""
Backend routing and reader factory.

The BackendRouter selects the best available backend for an algorithm
based on data size heuristics and the algorithm's declared capabilities.

Note: backend values are compared as strings ("hail", "pandas", "duckdb") so
this module stays free of hvantk.core.models imports (intra-core direction rule).
Callers pass Backend enum instances whose .value attribute gives the string.
"""

import logging
import os
from typing import Any, List, Optional

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

# Backend string constants (mirror Backend enum .value)
_BACKEND_HAIL = "hail"
_BACKEND_PANDAS = "pandas"
_BACKEND_DUCKDB = "duckdb"


def _backend_value(b: Any) -> str:
    """Return the string value of a Backend enum or plain string."""
    return b.value if hasattr(b, "value") else str(b)


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
        meta: Any,
        data_paths: Optional[List[str]] = None,
    ) -> Any:
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
            preference = [_BACKEND_PANDAS, _BACKEND_DUCKDB, _BACKEND_HAIL]
        elif row_hint is not None and row_hint < self.large_threshold:
            preference = [_BACKEND_DUCKDB, _BACKEND_PANDAS, _BACKEND_HAIL]
        else:
            # Large data or unknown size → prefer Hail
            preference = [_BACKEND_HAIL, _BACKEND_DUCKDB, _BACKEND_PANDAS]

        # Pick from `available` the backend that ranks earliest in `preference`.
        # Backends not in `preference` are sorted to the end as a tiebreaker.
        def _rank(b):
            bv = _backend_value(b)
            try:
                return preference.index(bv)
            except ValueError:
                return len(preference)

        selected = min(available, key=_rank)
        logger.info(
            "BackendRouter: selected %s (row_hint=%s, available=%s)",
            _backend_value(selected),
            row_hint,
            [_backend_value(b) for b in available],
        )
        return selected

    @staticmethod
    def _filter_available(backends: List[Any]) -> List[Any]:
        """Return backends whose dependencies are importable."""
        available = []
        for b in backends:
            bv = _backend_value(b)
            if bv == _BACKEND_PANDAS:
                available.append(b)  # pandas is always available
            elif bv == _BACKEND_HAIL:
                try:
                    import hail  # noqa: F401

                    available.append(b)
                except ImportError:
                    pass
            elif bv == _BACKEND_DUCKDB:
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
    def create(backend: Any) -> DataReader:
        bv = _backend_value(backend)
        if bv == _BACKEND_HAIL:
            return HailReader()
        elif bv == _BACKEND_PANDAS:
            return PandasReader()
        elif bv == _BACKEND_DUCKDB:
            from hvantk.core.utils.readers import DuckDBReader

            return DuckDBReader()
        raise ValueError(f"Unknown backend: {backend}")
