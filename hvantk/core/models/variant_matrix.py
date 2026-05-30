"""VariantMatrix: multi-sample variant cohort backed by a Hail MatrixTable.

The natural artifact for genotype data — variants × samples × multi-field
entries (GT, AD, DP, GQ, PL, ...). Hail-only by design; AnnData is the wrong
shape for cohort genotype data at TB scale.

`samples` and `variants` accessors return AnnotationTable views over
`mt.cols()` / `mt.rows()`. Entry-field access is deliberately not exposed —
callers go through `to_hail_mt()` for any operation that touches entries.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from hvantk.core.models._expr import Expr
from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.provenance import Provenance


@dataclass
class VariantMatrix:
    provenance: Provenance
    _mt: Any = field(repr=False)
    _n_samples_cached: int | None = field(default=None, repr=False)
    _n_variants_cached: int | None = field(default=None, repr=False)

    # --- constructors ---

    @classmethod
    def from_hail_mt(cls, mt: Any, *, provenance: Provenance) -> "VariantMatrix":
        return cls(provenance=provenance, _mt=mt)

    # --- counts ---

    @property
    def n_samples(self) -> int:
        if self._n_samples_cached is None:
            self._n_samples_cached = self._mt.count_cols()
        return self._n_samples_cached

    @property
    def n_variants(self) -> int:
        if self._n_variants_cached is None:
            self._n_variants_cached = self._mt.count_rows()
        return self._n_variants_cached

    # --- metadata accessors ---

    @property
    def samples(self) -> AnnotationTable:
        return AnnotationTable.from_hail(self._mt.cols(), provenance=self.provenance)

    @property
    def variants(self) -> AnnotationTable:
        return AnnotationTable.from_hail(self._mt.rows(), provenance=self.provenance)

    # --- subsetting ---

    def subset_samples(self, predicate: Expr) -> "VariantMatrix":
        from hvantk.core.models._compile import compile_to_hail_mt_col
        hail_expr = compile_to_hail_mt_col(predicate, self._mt)
        return VariantMatrix.from_hail_mt(
            self._mt.filter_cols(hail_expr), provenance=self.provenance,
        )

    def subset_variants(self, predicate: Expr) -> "VariantMatrix":
        from hvantk.core.models._compile import compile_to_hail_mt_row
        hail_expr = compile_to_hail_mt_row(predicate, self._mt)
        return VariantMatrix.from_hail_mt(
            self._mt.filter_rows(hail_expr), provenance=self.provenance,
        )

    # --- escape hatch ---

    def to_hail_mt(self) -> Any:
        """Return the underlying Hail MatrixTable. No copy is made."""
        return self._mt

    # --- persistence ---

    def save(self, path: str | Path) -> None:
        from hvantk.core import io as core_io
        core_io.save(self, path)

    @classmethod
    def load(cls, path: str | Path) -> "VariantMatrix":
        from hvantk.core import io as core_io
        result = core_io.load(path)
        if not isinstance(result, cls):
            raise TypeError(
                f"{path} contains a {type(result).__name__}, not an {cls.__name__}"
            )
        return result
