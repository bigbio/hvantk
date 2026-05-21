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
from hvantk.core.models._expr import AggOp, Expr  # noqa: F401


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

    def select(self, *columns: str) -> "AnnotationTable":
        if self.backend == "pandas":
            return AnnotationTable.from_pandas(
                self._table[list(columns)].copy(), provenance=self.provenance
            )
        return AnnotationTable.from_hail(
            self._table.select(*columns), provenance=self.provenance
        )

    def with_columns(self, **assignments: "Expr") -> "AnnotationTable":
        from hvantk.core.models._compile import (
            compile_to_hail,
            compile_to_pandas,
        )

        if self.backend == "pandas":
            df = self._table.copy()
            for name, expr in assignments.items():
                df[name] = compile_to_pandas(expr, df)
            return AnnotationTable.from_pandas(df, provenance=self.provenance)
        # hail
        kwargs = {
            name: compile_to_hail(expr, self._table)
            for name, expr in assignments.items()
        }
        return AnnotationTable.from_hail(
            self._table.annotate(**kwargs), provenance=self.provenance
        )

    def rename(self, **mapping: str) -> "AnnotationTable":
        if self.backend == "pandas":
            return AnnotationTable.from_pandas(
                self._table.rename(columns=mapping), provenance=self.provenance
            )
        return AnnotationTable.from_hail(
            self._table.rename(mapping), provenance=self.provenance
        )

    def join(
        self, other: "AnnotationTable", on: str | list[str], how: str = "inner",
        *, suffixes: tuple[str, str] | None = None,
    ) -> "AnnotationTable":
        if self.backend != other.backend:
            raise ValueError(
                f"cannot join AnnotationTables on different backends "
                f"({self.backend} vs {other.backend}); convert one first"
            )

        # Detect overlapping non-key column names — pandas would auto-suffix,
        # Hail would raise. Make the failure mode uniform: raise eagerly unless
        # the caller explicitly passed `suffixes` to opt in to renaming.
        keys = {on} if isinstance(on, str) else set(on)
        self_cols = set(self.schema.keys()) - keys
        other_cols = set(other.schema.keys()) - keys
        overlap = self_cols & other_cols
        if overlap and suffixes is None:
            raise ValueError(
                f"join: non-key columns overlap on both sides: {sorted(overlap)}. "
                f"Either rename them before joining, or pass suffixes=('_left', '_right') "
                f"to opt in to pandas-style suffixing (Hail backend does not support suffixes; "
                f"rename is the portable fix)."
            )

        if self.backend == "pandas":
            merge_kwargs: dict[str, Any] = {
                "on": on if isinstance(on, str) else list(on),
                "how": how,
            }
            if suffixes is not None:
                merge_kwargs["suffixes"] = suffixes
            merged = self._table.merge(other._table, **merge_kwargs)
            return AnnotationTable.from_pandas(merged, provenance=self.provenance)
        # hail: rekey both sides defensively (key_by is idempotent)
        keys_list = [on] if isinstance(on, str) else list(on)
        self_keyed = self._table.key_by(*keys_list)
        other_keyed = other._table.key_by(*keys_list)
        joined = self_keyed.join(other_keyed, how=how)
        return AnnotationTable.from_hail(joined, provenance=self.provenance)

    def distinct(self, subset: list[str] | None = None) -> "AnnotationTable":
        if self.backend == "pandas":
            df = self._table.drop_duplicates(subset=subset).reset_index(drop=True)
            return AnnotationTable.from_pandas(df, provenance=self.provenance)
        # hail
        if subset:
            return AnnotationTable.from_hail(
                self._table.key_by(*subset).distinct(), provenance=self.provenance
            )
        return AnnotationTable.from_hail(
            self._table.distinct(), provenance=self.provenance
        )

    def head(self, n: int = 5) -> "AnnotationTable":
        if self.backend == "pandas":
            return AnnotationTable.from_pandas(
                self._table.head(n).reset_index(drop=True), provenance=self.provenance
            )
        return AnnotationTable.from_hail(
            self._table.head(n), provenance=self.provenance
        )

    def group_by(self, *columns: str) -> "_GroupedAnnotationTable":
        return _GroupedAnnotationTable(self, list(columns))

    # --- terminal operations ---

    def collect(self) -> list[dict]:
        if self.backend == "pandas":
            return self._table.to_dict(orient="records")
        return [dict(r) for r in self._table.collect()]

    def count(self) -> int:
        if self.backend == "pandas":
            return len(self._table)
        return self._table.count()

    # --- persistence (stubbed; lands in Task 13) ---

    def save(self, path: str | Path) -> None:
        from hvantk.core import io as core_io
        core_io.save(self, path)


class _GroupedAnnotationTable:
    """Intermediate handle between group_by() and agg(). Not constructed directly."""

    def __init__(self, parent: AnnotationTable, columns: list[str]) -> None:
        self._parent = parent
        self._columns = columns

    def agg(self, **aggregations: "AggOp") -> AnnotationTable:
        from hvantk.core.models._compile import (
            compile_agg_to_hail,
            compile_agg_to_pandas,
        )

        if self._parent.backend == "pandas":
            df = self._parent._table
            rows = []
            for keys, group_df in df.groupby(self._columns):
                # groupby(list) always yields a tuple of keys
                if not isinstance(keys, tuple):
                    keys = (keys,)
                out = dict(zip(self._columns, keys))
                for name, agg in aggregations.items():
                    out[name] = compile_agg_to_pandas(agg, group_df)
                rows.append(out)
            result = pd.DataFrame(rows)
            return AnnotationTable.from_pandas(
                result, provenance=self._parent.provenance
            )
        # hail
        ht = self._parent._table
        kwargs = {
            name: compile_agg_to_hail(agg, ht) for name, agg in aggregations.items()
        }
        grouped_ht = ht.group_by(*self._columns).aggregate(**kwargs)
        return AnnotationTable.from_hail(grouped_ht, provenance=self._parent.provenance)
