"""ExpressionMatrix: samples × features expression data, backed by AnnData.

Originally dual-backend (anndata + hail-mt), but in practice the hail-mt
backend was always playing the role of a multi-sample variant cohort —
that role is now held by VariantMatrix. ExpressionMatrix is the
appropriate artifact for expression data; VariantMatrix is the
appropriate artifact for genotype cohorts.

obs and var are exposed as AnnotationTable instances so the same Expr DSL
works against them as against any other AnnotationTable.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from hvantk.core.models._compile import compile_to_pandas
from hvantk.core.models._expr import Expr
from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.provenance import Provenance


def _materialize_metadata_df(df: pd.DataFrame, axis_id_name: str) -> pd.DataFrame:
    """Reset a metadata DataFrame's index to an explicit column named ``axis_id_name``."""
    out = df.copy()
    out.index = out.index.rename(None)
    return out.reset_index().rename(columns={"index": axis_id_name})


@dataclass
class ExpressionMatrix:
    provenance: Provenance
    _matrix: Any = field(repr=False)
    _n_obs_cached: int | None = field(default=None, repr=False)
    _n_vars_cached: int | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        if self._n_obs_cached is None or self._n_vars_cached is None:
            raise ValueError(
                "ExpressionMatrix requires _n_obs_cached and _n_vars_cached; "
                "use ExpressionMatrix.from_anndata() to construct."
            )

    @property
    def n_obs(self) -> int:
        return self._n_obs_cached

    @property
    def n_vars(self) -> int:
        return self._n_vars_cached

    @classmethod
    def from_anndata(cls, adata: ad.AnnData, *, provenance: Provenance) -> "ExpressionMatrix":
        return cls(
            provenance=provenance,
            _matrix=adata,
            _n_obs_cached=adata.n_obs,
            _n_vars_cached=adata.n_vars,
        )

    @property
    def obs(self) -> AnnotationTable:
        return AnnotationTable.from_pandas(
            _materialize_metadata_df(self._matrix.obs, "obs_id"),
            provenance=self.provenance,
        )

    @property
    def var(self) -> AnnotationTable:
        return AnnotationTable.from_pandas(
            _materialize_metadata_df(self._matrix.var, "var_id"),
            provenance=self.provenance,
        )

    def subset_obs(self, predicate: Expr) -> "ExpressionMatrix":
        obs_df = _materialize_metadata_df(self._matrix.obs, "obs_id")
        mask = compile_to_pandas(predicate, obs_df)
        return ExpressionMatrix.from_anndata(
            self._matrix[mask.values, :].copy(), provenance=self.provenance,
        )

    def subset_var(self, predicate: Expr) -> "ExpressionMatrix":
        var_df = _materialize_metadata_df(self._matrix.var, "var_id")
        mask = compile_to_pandas(predicate, var_df)
        return ExpressionMatrix.from_anndata(
            self._matrix[:, mask.values].copy(), provenance=self.provenance,
        )

    def X(self, layer: str | None = None) -> Any:
        return self._matrix.X if layer is None else self._matrix.layers[layer]

    def layers(self) -> dict[str, Any]:
        return dict(self._matrix.layers)

    def aggregate_obs(self, by: str, func: str = "mean") -> "ExpressionMatrix":
        grouped_df = self._matrix.obs.groupby(by, observed=True).indices
        rows = []
        idx = []
        for grp_value, row_idx in grouped_df.items():
            X = self._matrix.X[row_idx, :]
            if hasattr(X, "toarray"):
                X = X.toarray()
            rows.append(getattr(np, func)(X, axis=0))
            idx.append(grp_value)
        X_out = np.vstack(rows)
        new_obs = pd.DataFrame({by: idx}).set_index(pd.Index(idx, name=by))
        new = ad.AnnData(X=X_out, obs=new_obs, var=self._matrix.var.copy())
        return ExpressionMatrix.from_anndata(new, provenance=self.provenance)

    def to_anndata(self) -> ad.AnnData:
        """Return a copy of the underlying AnnData. Mutations to the result do not affect this artifact."""
        return self._matrix.copy()

    def save(self, path: str | Path) -> None:
        from hvantk.core import io as core_io
        core_io.save(self, path)

    @classmethod
    def load(cls, path: "str | Path") -> "ExpressionMatrix":
        from hvantk.core import io as core_io
        result = core_io.load(path)
        if not isinstance(result, cls):
            raise TypeError(
                f"{path} contains a {type(result).__name__}, not an {cls.__name__}"
            )
        return result
