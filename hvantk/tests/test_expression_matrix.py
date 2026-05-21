"""Tests for ExpressionMatrix: construction, obs/var access, subsetting."""
from __future__ import annotations

from datetime import datetime, timezone

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from hvantk.core.models._expr import col
from hvantk.core.models.expression_matrix import ExpressionMatrix
from hvantk.core.models.provenance import Provenance


def _prov():
    return Provenance(
        plugin="t",
        dataset="t:expr",
        plugin_version="0.0",
        source_fingerprint="sha256:abc",
        schema_id="t-expr-v1",
        build_timestamp=datetime(2026, 5, 20, tzinfo=timezone.utc),
        builder_commit=None,
    )


@pytest.fixture()
def adata():
    obs = pd.DataFrame(
        {"tissue": ["liver", "brain", "liver"]},
        index=["s1", "s2", "s3"],
    )
    var = pd.DataFrame(
        {"gene": ["BRCA1", "TP53"]},
        index=["g1", "g2"],
    )
    X = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    return ad.AnnData(X=X, obs=obs, var=var)


def test_from_anndata(adata):
    em = ExpressionMatrix.from_anndata(adata, provenance=_prov())
    assert em.backend == "anndata"
    assert em.n_obs == 3
    assert em.n_vars == 2
    assert em.provenance == _prov()


def test_obs_and_var_are_annotation_tables(adata):
    em = ExpressionMatrix.from_anndata(adata, provenance=_prov())
    assert em.obs.collect()[0]["tissue"] == "liver"
    assert em.var.collect()[0]["gene"] == "BRCA1"


def test_subset_obs(adata):
    em = ExpressionMatrix.from_anndata(adata, provenance=_prov())
    sub = em.subset_obs(col("tissue") == "liver")
    assert sub.n_obs == 2
    assert sub.n_vars == 2


def test_subset_var(adata):
    em = ExpressionMatrix.from_anndata(adata, provenance=_prov())
    sub = em.subset_var(col("gene") == "BRCA1")
    assert sub.n_obs == 3
    assert sub.n_vars == 1


def test_X_returns_array(adata):
    em = ExpressionMatrix.from_anndata(adata, provenance=_prov())
    X = em.X()
    assert X.shape == (3, 2)
    assert X[0, 0] == 1.0


@pytest.mark.hail
def test_to_hail_mt_returns_matrix_table(adata):
    """Phase J: anndata backend converts to a Hail MatrixTable without error."""
    em = ExpressionMatrix.from_anndata(adata, provenance=_prov())
    mt = em.to_hail_mt()
    import hail as hl
    assert isinstance(mt, hl.MatrixTable)


def test_subset_obs_using_obs_id(adata):
    """Predicates referencing the synthetic obs_id column (exposed via .obs)
    must also work in subset_obs. Otherwise the two APIs disagree."""
    em = ExpressionMatrix.from_anndata(adata, provenance=_prov())
    sub = em.subset_obs(col("obs_id") == "s1")
    assert sub.n_obs == 1


def test_subset_var_using_var_id(adata):
    em = ExpressionMatrix.from_anndata(adata, provenance=_prov())
    sub = em.subset_var(col("var_id") == "g1")
    assert sub.n_vars == 1


def test_subset_obs_with_named_index(tmp_path):
    """When AnnData index has a name (e.g. 'sample'), obs_id is still the synthetic column."""
    obs = pd.DataFrame(
        {"tissue": ["liver", "brain"]},
        index=pd.Index(["s1", "s2"], name="sample"),  # NAMED index
    )
    var = pd.DataFrame({"gene": ["BRCA1"]}, index=["g1"])
    adata = ad.AnnData(X=np.array([[1.0], [2.0]]), obs=obs, var=var)
    em = ExpressionMatrix.from_anndata(adata, provenance=_prov())

    # obs column should be obs_id regardless of the original index name
    rows = em.obs.collect()
    assert "obs_id" in rows[0]
    assert rows[0]["obs_id"] == "s1"

    sub = em.subset_obs(col("obs_id") == "s1")
    assert sub.n_obs == 1


def test_subset_var_with_named_index(tmp_path):
    obs = pd.DataFrame({"tissue": ["liver"]}, index=["s1"])
    var = pd.DataFrame(
        {"gene": ["BRCA1", "TP53"]},
        index=pd.Index(["g1", "g2"], name="ensembl_id"),  # NAMED index
    )
    adata = ad.AnnData(X=np.array([[1.0, 2.0]]), obs=obs, var=var)
    em = ExpressionMatrix.from_anndata(adata, provenance=_prov())

    sub = em.subset_var(col("var_id") == "g1")
    assert sub.n_vars == 1
