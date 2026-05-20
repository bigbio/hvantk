"""AnnotationTable: a keyed annotation table, backend-agnostic.

Algorithms consume AnnotationTable instances and use the portable Expr API
to filter / select / join. Backend (hail vs pandas) is an implementation
detail — algorithms must not branch on it.

This file delivers construction, conversion, and identity. The portable
query API (filter, select, etc.) lands in subsequent tasks.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

import pandas as pd

from hvantk.core.models.provenance import Provenance
from hvantk.core.models._expr import Expr  # noqa: F401 (used in type hints only)


_BACKENDS = ("hail", "pandas")


def _normalize_dtype(dtype: Any) -> str:
    """Stringify a pandas / hail dtype to a stable schema label."""
    s = str(dtype).lower()
    if "int" in s:
        return "int"
    if "float" in s or "double" in s:
        return "float"
    if "bool" in s:
        return "bool"
    return "str"


@dataclass
class AnnotationTable:
    backend: Literal["hail", "pandas"]
    provenance: Provenance
    schema: dict[str, str]
    _table: Any = field(repr=False)

    def __post_init__(self) -> None:
        if self.backend not in _BACKENDS:
            raise ValueError(f"backend must be one of {_BACKENDS}; got {self.backend!r}")

    # --- constructors ---

    @classmethod
    def from_pandas(cls, df: pd.DataFrame, *, provenance: Provenance) -> "AnnotationTable":
        schema = {col: _normalize_dtype(dtype) for col, dtype in df.dtypes.items()}
        return cls(backend="pandas", provenance=provenance, schema=schema, _table=df)

    @classmethod
    def from_hail(cls, ht: Any, *, provenance: Provenance) -> "AnnotationTable":
        schema = {name: _normalize_dtype(t) for name, t in ht.row.dtype.items()}
        return cls(backend="hail", provenance=provenance, schema=schema, _table=ht)

    # --- escape hatches ---

    def to_pandas(self) -> pd.DataFrame:
        if self.backend == "pandas":
            return self._table.copy()
        return self._table.to_pandas()

    def to_hail(self) -> Any:
        import hail as hl

        if self.backend == "hail":
            return self._table
        return hl.Table.from_pandas(self._table)

    # --- query operations (immutable, return new AnnotationTable) ---

    def filter(self, predicate: "Expr") -> "AnnotationTable":
        from hvantk.core.models._compile import (
            compile_to_hail,
            compile_to_pandas,
        )

        if self.backend == "pandas":
            mask = compile_to_pandas(predicate, self._table)
            return AnnotationTable.from_pandas(
                self._table[mask].reset_index(drop=True),
                provenance=self.provenance,
            )
        # hail
        hail_expr = compile_to_hail(predicate, self._table)
        return AnnotationTable.from_hail(
            self._table.filter(hail_expr), provenance=self.provenance
        )

    # --- persistence (stubbed; lands in Task 13) ---

    def save(self, path: str | Path) -> None:
        raise NotImplementedError("save lands in Task 13 alongside core/io")
