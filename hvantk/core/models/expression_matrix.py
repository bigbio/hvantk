"""ExpressionMatrix: samples × features with sample/feature metadata.

Backend choices:
  - ``anndata``: backed by an :class:`anndata.AnnData` object (in-memory).
  - ``hail-mt``: backed by a Hail :class:`~hail.MatrixTable` (distributed).

obs and var are exposed as AnnotationTable instances so the same Expr DSL
works across both backends.

Phase J notes
-------------
- ``from_hail_mt`` / ``to_hail_mt`` are fully implemented.
- ``to_anndata()`` from a hail-mt backend is a lossy dense materialization
  intended for small fixtures only.
- ``to_hail_mt()`` from an anndata backend delegates to
  :func:`hvantk.core.utils.converters.anndata_to_hail_mt`; metadata
  (obs/var columns) is propagated but obs/var key fields are named
  ``sample_id`` / ``gene_id`` by default.
- ``X()`` on a hail-mt backend triggers an eager Hail ``entries()``
  collect — use ``to_hail_mt()`` for large-scale, distributed operations.
- ``subset_obs`` / ``subset_var`` compile the predicate against the
  MatrixTable's cols/rows Table (not the MatrixTable directly) to avoid
  Hail's key-based ``mt[field]`` semantics, then feed the resulting
  boolean expression to ``filter_cols`` / ``filter_rows``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

import anndata as ad
import pandas as pd

from hvantk.core.models._compile import compile_to_pandas
from hvantk.core.models._expr import Expr
from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.provenance import Provenance


_BACKENDS = ("anndata", "hail-mt")


def _materialize_metadata_df(df: pd.DataFrame, axis_id_name: str) -> pd.DataFrame:
    """Reset a metadata DataFrame's index to an explicit column named ``axis_id_name``.

    AnnData obs/var DataFrames may have either a named or unnamed index.
    This helper produces a stable schema with ``axis_id_name`` regardless.
    """
    out = df.copy()
    out.index = out.index.rename(None)
    return out.reset_index().rename(columns={"index": axis_id_name})


def _entry_field_name(mt: Any) -> str:
    """Return the canonical entry field name from a Hail MatrixTable.

    Prefers ``value``, then ``x``, then the first available field.
    Raises ``ValueError`` if there are no entry fields.
    """
    entry_fields = list(mt.entry)
    if not entry_fields:
        raise ValueError("MatrixTable has no entry fields")
    for preferred in ("value", "x", "X"):
        if preferred in entry_fields:
            return preferred
    return entry_fields[0]


@dataclass
class ExpressionMatrix:
    backend: Literal["anndata", "hail-mt"]
    provenance: Provenance
    n_obs: int
    n_vars: int
    _matrix: Any = field(repr=False)

    def __post_init__(self) -> None:
        if self.backend not in _BACKENDS:
            raise ValueError(f"backend must be one of {_BACKENDS}; got {self.backend!r}")

    # --- constructors ---

    @classmethod
    def from_anndata(
        cls, adata: ad.AnnData, *, provenance: Provenance
    ) -> "ExpressionMatrix":
        return cls(
            backend="anndata",
            provenance=provenance,
            n_obs=adata.n_obs,
            n_vars=adata.n_vars,
            _matrix=adata,
        )

    @classmethod
    def from_hail_mt(
        cls, mt: Any, *, provenance: Provenance
    ) -> "ExpressionMatrix":
        """Construct an ExpressionMatrix from a Hail MatrixTable.

        The MatrixTable's row fields become var metadata; column fields become
        obs metadata.  ``n_obs`` and ``n_vars`` are resolved eagerly via
        ``count_cols()`` / ``count_rows()`` — Hail triggers a computation here.

        Parameters
        ----------
        mt:
            A ``hail.MatrixTable`` instance.
        provenance:
            Provenance record for this artifact.
        """
        return cls(
            backend="hail-mt",
            provenance=provenance,
            n_obs=mt.count_cols(),
            n_vars=mt.count_rows(),
            _matrix=mt,
        )

    # --- metadata accessors ---

    @property
    def obs(self) -> AnnotationTable:
        if self.backend == "anndata":
            return AnnotationTable.from_pandas(
                _materialize_metadata_df(self._matrix.obs, "obs_id"),
                provenance=self.provenance,
            )
        # hail-mt: cols() returns a Hail Table of column metadata
        cols_ht = self._matrix.cols()
        return AnnotationTable.from_hail(cols_ht, provenance=self.provenance)

    @property
    def var(self) -> AnnotationTable:
        if self.backend == "anndata":
            return AnnotationTable.from_pandas(
                _materialize_metadata_df(self._matrix.var, "var_id"),
                provenance=self.provenance,
            )
        # hail-mt: rows() returns a Hail Table of row metadata
        rows_ht = self._matrix.rows()
        return AnnotationTable.from_hail(rows_ht, provenance=self.provenance)

    # --- subsetting ---

    def subset_obs(self, predicate: Expr) -> "ExpressionMatrix":
        if self.backend == "anndata":
            obs_df = _materialize_metadata_df(self._matrix.obs, "obs_id")
            mask = compile_to_pandas(predicate, obs_df)
            return ExpressionMatrix.from_anndata(
                self._matrix[mask.values, :].copy(), provenance=self.provenance
            )
        # hail-mt: compile against the MatrixTable's col scope so the resulting
        # expression is bound to the MT (not a separate Table), then filter_cols.
        # compile_to_hail_mt_col uses mt.col["field"] which returns an expression
        # in the correct scope for filter_cols().
        from hvantk.core.models._compile import compile_to_hail_mt_col
        hail_expr = compile_to_hail_mt_col(predicate, self._matrix)
        mt_filtered = self._matrix.filter_cols(hail_expr)
        return ExpressionMatrix.from_hail_mt(mt_filtered, provenance=self.provenance)

    def subset_var(self, predicate: Expr) -> "ExpressionMatrix":
        if self.backend == "anndata":
            var_df = _materialize_metadata_df(self._matrix.var, "var_id")
            mask = compile_to_pandas(predicate, var_df)
            return ExpressionMatrix.from_anndata(
                self._matrix[:, mask.values].copy(), provenance=self.provenance
            )
        # hail-mt: compile against the MatrixTable's row scope so the resulting
        # expression is bound to the MT (not a separate Table), then filter_rows.
        # compile_to_hail_mt_row uses mt.row["field"] which returns an expression
        # in the correct scope for filter_rows().
        from hvantk.core.models._compile import compile_to_hail_mt_row
        hail_expr = compile_to_hail_mt_row(predicate, self._matrix)
        mt_filtered = self._matrix.filter_rows(hail_expr)
        return ExpressionMatrix.from_hail_mt(mt_filtered, provenance=self.provenance)

    # --- value access ---

    def X(self, layer: str | None = None) -> Any:
        """Return the expression matrix.

        For the ``hail-mt`` backend this materializes the full matrix into
        a NumPy array via ``entries().to_pandas()`` — intended for small
        fixtures only.  Large-scale consumers should use ``to_hail_mt()``
        and Hail-native distributed operations.

        Returns an array shaped ``(n_obs, n_vars)``.
        """
        if self.backend == "anndata":
            return self._matrix.X if layer is None else self._matrix.layers[layer]
        # hail-mt
        import numpy as np
        field_name = layer if layer is not None else _entry_field_name(self._matrix)
        # Collect the full entries table; row_key and col_key fields are included.
        df = self._matrix.entries().select(field_name).to_pandas()
        # entries() returns one row per (row_key, col_key) pair in key order.
        # Pivot so that the result is (n_obs, n_vars) = (n_cols, n_rows).
        row_key_fields = list(self._matrix.row_key)
        col_key_fields = list(self._matrix.col_key)
        row_key = row_key_fields[0]
        col_key = col_key_fields[0]
        pivot = df.pivot(index=col_key, columns=row_key, values=field_name)
        return pivot.values.astype(np.float64)

    def layers(self) -> dict[str, Any]:
        if self.backend == "anndata":
            return dict(self._matrix.layers)
        # hail-mt: every entry field is a layer
        return {name: self.X(layer=name) for name in self._matrix.entry}

    # --- aggregation ---

    def aggregate_obs(self, by: str, func: str = "mean") -> "ExpressionMatrix":
        if self.backend == "anndata":
            grouped_df = self._matrix.obs.groupby(by, observed=True).indices
            import numpy as np

            rows = []
            idx = []
            for grp_value, row_idx in grouped_df.items():
                X = self._matrix.X[row_idx, :]
                if hasattr(X, "toarray"):  # sparse
                    X = X.toarray()
                rows.append(getattr(np, func)(X, axis=0))
                idx.append(grp_value)
            X_out = np.vstack(rows)
            new_obs = pd.DataFrame({by: idx}).set_index(pd.Index(idx, name=by))
            new = ad.AnnData(X=X_out, obs=new_obs, var=self._matrix.var.copy())
            return ExpressionMatrix.from_anndata(new, provenance=self.provenance)
        # hail-mt: group_cols_by + aggregate_entries
        import hail as hl
        field_name = _entry_field_name(self._matrix)
        grouped = self._matrix.group_cols_by(self._matrix[by]).aggregate(
            **{field_name: getattr(hl.agg, func)(self._matrix[field_name])}
        )
        return ExpressionMatrix.from_hail_mt(grouped, provenance=self.provenance)

    # --- backend conversions ---

    def to_anndata(self) -> ad.AnnData:
        """Return an AnnData representation.

        For the ``hail-mt`` backend this is a dense materialization; all
        entry values, row, and col metadata are collected into memory.
        Use only for small fixtures.
        """
        if self.backend == "anndata":
            return self._matrix.copy()
        # hail-mt → anndata via the shared converter utility
        from hvantk.core.utils.converters import hail_mt_to_anndata
        field_name = _entry_field_name(self._matrix)
        return hail_mt_to_anndata(self._matrix, entry_field=field_name)

    def to_hail_mt(self) -> Any:
        """Return a Hail MatrixTable representation.

        For the ``anndata`` backend this converts using
        :func:`hvantk.core.utils.converters.anndata_to_hail_mt`.
        The conversion is intentionally simplified: obs/var index values are
        propagated as ``sample_id`` / ``gene_id`` row/col keys, and obs/var
        column metadata is annotated as col/row fields.  Entry values are
        taken from ``adata.X`` and stored as the ``x`` entry field.

        For the ``hail-mt`` backend this returns the underlying MatrixTable
        directly (no copy).
        """
        if self.backend == "hail-mt":
            return self._matrix
        # anndata → hail-mt via shared converter
        from hvantk.core.utils.converters import anndata_to_hail_mt
        return anndata_to_hail_mt(self._matrix)

    # --- persistence ---

    def save(self, path: str | Path) -> None:
        from hvantk.core import io as core_io
        core_io.save(self, path)
