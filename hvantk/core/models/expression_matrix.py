"""ExpressionMatrix: samples × features with sample/feature metadata.

Backend choices: anndata (the only supported backend in Phase A) or
hail-mt (raises NotImplementedError until Phase J). obs and var are
exposed as AnnotationTable instances so the same Expr DSL works.
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

    @property
    def obs(self) -> AnnotationTable:
        if self.backend == "anndata":
            return AnnotationTable.from_pandas(
                self._matrix.obs.reset_index().rename(columns={"index": "obs_id"}),
                provenance=self.provenance,
            )
        raise NotImplementedError("hail-mt backend lands in Phase J")

    @property
    def var(self) -> AnnotationTable:
        if self.backend == "anndata":
            return AnnotationTable.from_pandas(
                self._matrix.var.reset_index().rename(columns={"index": "var_id"}),
                provenance=self.provenance,
            )
        raise NotImplementedError("hail-mt backend lands in Phase J")

    def subset_obs(self, predicate: Expr) -> "ExpressionMatrix":
        if self.backend == "anndata":
            obs_df = self._matrix.obs.reset_index().rename(columns={"index": "obs_id"})
            mask = compile_to_pandas(predicate, obs_df)
            return ExpressionMatrix.from_anndata(
                self._matrix[mask.values, :].copy(), provenance=self.provenance
            )
        raise NotImplementedError("hail-mt backend lands in Phase J")

    def subset_var(self, predicate: Expr) -> "ExpressionMatrix":
        if self.backend == "anndata":
            var_df = self._matrix.var.reset_index().rename(columns={"index": "var_id"})
            mask = compile_to_pandas(predicate, var_df)
            return ExpressionMatrix.from_anndata(
                self._matrix[:, mask.values].copy(), provenance=self.provenance
            )
        raise NotImplementedError("hail-mt backend lands in Phase J")

    def X(self, layer: str | None = None) -> Any:
        if self.backend == "anndata":
            return self._matrix.X if layer is None else self._matrix.layers[layer]
        raise NotImplementedError("hail-mt backend lands in Phase J")

    def layers(self) -> dict[str, Any]:
        if self.backend == "anndata":
            return dict(self._matrix.layers)
        raise NotImplementedError("hail-mt backend lands in Phase J")

    def aggregate_obs(self, by: str, func: str = "mean") -> "ExpressionMatrix":
        if self.backend != "anndata":
            raise NotImplementedError("hail-mt backend lands in Phase J")
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

    def to_anndata(self) -> ad.AnnData:
        if self.backend == "anndata":
            return self._matrix.copy()
        raise NotImplementedError("hail-mt backend lands in Phase J")

    def to_hail_mt(self) -> Any:
        raise NotImplementedError("hail-mt backend lands in Phase J")

    def save(self, path: str | Path) -> None:
        from hvantk.core import io as core_io
        core_io.save(self, path)
