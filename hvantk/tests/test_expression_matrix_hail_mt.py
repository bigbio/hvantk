"""Phase J: tests for the hail-mt backend of ExpressionMatrix.

All tests in this file are automatically marked ``@pytest.mark.hail``
because the filename ends in ``_hail_mt.py`` — the conftest auto-applies
the ``hail`` marker to filenames ending in ``_hail.py``.  For filenames
that don't match the suffix, we mark them explicitly.

Actually: the conftest suffix-to-mark table maps ``_hail.py`` suffix, but
this file ends in ``_hail_mt.py`` which also ends with ``_hail.py``'s
pattern? No — ``_hail_mt.py`` does NOT end with ``_hail.py``.  We
therefore mark every test here explicitly with ``@pytest.mark.hail``.
"""
from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pytest

from hvantk.core.models._expr import col
from hvantk.core.models.expression_matrix import ExpressionMatrix
from hvantk.core.models.provenance import Provenance


def _prov():
    return Provenance(
        plugin="t",
        dataset="t:mt",
        plugin_version="0.0",
        source_fingerprint="sha256:abc",
        schema_id="t-mt-v1",
        build_timestamp=datetime(2026, 5, 21, tzinfo=timezone.utc),
        builder_commit=None,
    )


@pytest.fixture
def small_mt(hail_session):
    """A 3×2 Hail MatrixTable (3 variants / rows × 2 samples / cols).

    Entry value layout (row_idx × col_idx):
        col 0  col 1
    row 0:  0.0    1.0
    row 1: 10.0   11.0
    row 2: 20.0   21.0

    After ExpressionMatrix wrapping (n_obs=2, n_vars=3), X() should return
    a (2, 3) array:
        var0  var1  var2
    obs0:  0.0  10.0  20.0
    obs1:  1.0  11.0  21.0
    """
    import hail as hl
    rows = []
    for r in range(3):   # rows = variants / vars
        for c in range(2):   # cols = samples / obs
            rows.append({"row_idx": r, "col_idx": c, "value": float(r * 10 + c)})
    ht = hl.Table.parallelize(
        rows,
        hl.tstruct(row_idx=hl.tint32, col_idx=hl.tint32, value=hl.tfloat64),
    )
    return ht.to_matrix_table(
        row_key=["row_idx"],
        col_key=["col_idx"],
    )


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

@pytest.mark.hail
def test_from_hail_mt_attributes(small_mt):
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    assert em.backend == "hail-mt"
    assert em.n_obs == 2
    assert em.n_vars == 3
    assert em.provenance == _prov()


@pytest.mark.hail
def test_from_hail_mt_stores_matrix(small_mt):
    import hail as hl
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    assert isinstance(em._matrix, hl.MatrixTable)


# ---------------------------------------------------------------------------
# obs / var accessors
# ---------------------------------------------------------------------------

@pytest.mark.hail
def test_obs_is_annotation_table_hail(small_mt):
    from hvantk.core.models.annotation_table import AnnotationTable
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    obs = em.obs
    assert isinstance(obs, AnnotationTable)
    assert obs.backend == "hail"


@pytest.mark.hail
def test_var_is_annotation_table_hail(small_mt):
    from hvantk.core.models.annotation_table import AnnotationTable
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    var = em.var
    assert isinstance(var, AnnotationTable)
    assert var.backend == "hail"


@pytest.mark.hail
def test_obs_count_matches_n_obs(small_mt):
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    assert em.obs.count() == em.n_obs


@pytest.mark.hail
def test_var_count_matches_n_vars(small_mt):
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    assert em.var.count() == em.n_vars


# ---------------------------------------------------------------------------
# X() materialization
# ---------------------------------------------------------------------------

@pytest.mark.hail
def test_X_shape(small_mt):
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    X = em.X()
    assert X.shape == (2, 3)  # (n_obs, n_vars)


@pytest.mark.hail
def test_X_values(small_mt):
    """Check concrete values: obs0 = [0, 10, 20], obs1 = [1, 11, 21]."""
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    X = em.X()
    expected = np.array([[0.0, 10.0, 20.0], [1.0, 11.0, 21.0]])
    np.testing.assert_array_almost_equal(X, expected)


@pytest.mark.hail
def test_layers_returns_dict(small_mt):
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    layers = em.layers()
    assert isinstance(layers, dict)
    assert "value" in layers
    assert layers["value"].shape == (2, 3)


# ---------------------------------------------------------------------------
# to_anndata conversion
# ---------------------------------------------------------------------------

@pytest.mark.hail
def test_to_anndata_shape(small_mt):
    import anndata as ad
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    adata = em.to_anndata()
    assert isinstance(adata, ad.AnnData)
    assert adata.n_obs == 2
    assert adata.n_vars == 3


@pytest.mark.hail
def test_anndata_hail_mt_X_parity(small_mt):
    """X values are identical between hail-mt and anndata conversions."""
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    X_mt = em.X()
    adata = em.to_anndata()
    X_anndata = adata.X if not hasattr(adata.X, "toarray") else adata.X.toarray()
    assert X_mt.shape == X_anndata.shape
    np.testing.assert_array_almost_equal(
        X_mt.astype(np.float32), X_anndata.astype(np.float32)
    )


# ---------------------------------------------------------------------------
# to_hail_mt (anndata → hail-mt round-trip)
# ---------------------------------------------------------------------------

@pytest.mark.hail
def test_to_hail_mt_from_anndata_returns_matrix_table(hail_session):
    """anndata backend → to_hail_mt() returns a Hail MatrixTable."""
    import anndata as ad
    import hail as hl
    import pandas as pd

    obs = pd.DataFrame({"tissue": ["liver", "brain"]}, index=["s1", "s2"])
    var = pd.DataFrame({"gene": ["BRCA1", "TP53"]}, index=["g1", "g2"])
    X = np.array([[1.0, 2.0], [3.0, 4.0]])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    em = ExpressionMatrix.from_anndata(adata, provenance=_prov())
    mt = em.to_hail_mt()
    assert isinstance(mt, hl.MatrixTable)
    n_rows, n_cols = mt.count()
    assert n_rows == 2
    assert n_cols == 2


# ---------------------------------------------------------------------------
# subset_obs / subset_var
# ---------------------------------------------------------------------------

@pytest.mark.hail
def test_subset_obs_hail_mt(small_mt):
    """Filter cols (obs) by col_idx == 0 should leave 1 obs."""
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    sub = em.subset_obs(col("col_idx") == 0)
    assert sub.n_obs == 1
    assert sub.n_vars == 3


@pytest.mark.hail
def test_subset_var_hail_mt(small_mt):
    """Filter rows (var) by row_idx == 0 should leave 1 var."""
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    sub = em.subset_var(col("row_idx") == 0)
    assert sub.n_obs == 2
    assert sub.n_vars == 1


@pytest.mark.hail
def test_subset_obs_gt_predicate(small_mt):
    """col_idx > 0 selects obs with col_idx=1 only."""
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    sub = em.subset_obs(col("col_idx") > 0)
    assert sub.n_obs == 1


@pytest.mark.hail
def test_subset_var_gt_predicate(small_mt):
    """row_idx > 0 selects vars with row_idx=1 and row_idx=2."""
    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    sub = em.subset_var(col("row_idx") > 0)
    assert sub.n_vars == 2


# ---------------------------------------------------------------------------
# core/io round-trip: save → load
# ---------------------------------------------------------------------------

@pytest.mark.hail
def test_save_load_mt_round_trip(tmp_path, small_mt):
    """Save an ExpressionMatrix as .mt/ via core/io; load it back."""
    from hvantk.core import io as core_io

    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    out = tmp_path / "expr.mt"
    core_io.save(em, out)

    # The .mt directory must exist and a sidecar provenance file alongside it
    assert out.exists()
    assert (out.with_name("expr.mt.provenance.json")).exists()

    loaded = core_io.load(out)
    assert isinstance(loaded, ExpressionMatrix)
    assert loaded.backend == "hail-mt"
    assert loaded.provenance == em.provenance
    assert loaded.n_obs == 2
    assert loaded.n_vars == 3


@pytest.mark.hail
def test_save_load_mt_X_parity(tmp_path, small_mt):
    """Entry values survive a save/load round-trip."""
    from hvantk.core import io as core_io

    em = ExpressionMatrix.from_hail_mt(small_mt, provenance=_prov())
    out = tmp_path / "expr2.mt"
    core_io.save(em, out)

    loaded = core_io.load(out)
    np.testing.assert_array_almost_equal(em.X(), loaded.X())


@pytest.mark.hail
def test_save_load_anndata_backend_via_mt(tmp_path, hail_session):
    """anndata-backend ExpressionMatrix saved as .mt/ round-trips correctly."""
    import anndata as ad
    import pandas as pd
    from hvantk.core import io as core_io

    obs = pd.DataFrame({"tissue": ["liver"]}, index=["s1"])
    var = pd.DataFrame({"gene": ["BRCA1"]}, index=["g1"])
    adata = ad.AnnData(X=np.array([[1.0]]), obs=obs, var=var)
    em = ExpressionMatrix.from_anndata(adata, provenance=_prov())
    out = tmp_path / "anndata_as_mt.mt"
    core_io.save(em, out)

    loaded = core_io.load(out)
    assert isinstance(loaded, ExpressionMatrix)
    assert loaded.backend == "hail-mt"
    assert loaded.n_obs == 1
    assert loaded.n_vars == 1


# ---------------------------------------------------------------------------
# Lazy count behavior
# ---------------------------------------------------------------------------

@pytest.mark.hail
def test_from_hail_mt_does_not_eagerly_count(monkeypatch, hail_session):
    """from_hail_mt should not trigger count_cols / count_rows at construction."""
    import hail as hl

    rows = []
    for r in range(3):
        for c in range(2):
            rows.append({"row_idx": r, "col_idx": c, "value": float(r * 10 + c)})
    ht = hl.Table.parallelize(
        rows,
        hl.tstruct(row_idx=hl.tint32, col_idx=hl.tint32, value=hl.tfloat64),
    )
    mt = ht.to_matrix_table(row_key=["row_idx"], col_key=["col_idx"])

    # Track count calls via monkeypatching the type
    counts = {"cols": 0, "rows": 0}
    orig_count_cols = type(mt).count_cols
    orig_count_rows = type(mt).count_rows

    def spy_count_cols(self):
        counts["cols"] += 1
        return orig_count_cols(self)

    def spy_count_rows(self):
        counts["rows"] += 1
        return orig_count_rows(self)

    monkeypatch.setattr(type(mt), "count_cols", spy_count_cols)
    monkeypatch.setattr(type(mt), "count_rows", spy_count_rows)

    em = ExpressionMatrix.from_hail_mt(mt, provenance=_prov())

    # Construction must NOT have triggered any counts
    assert counts == {"cols": 0, "rows": 0}, (
        f"from_hail_mt triggered eager counts: {counts}"
    )

    # Accessing n_obs triggers exactly one count_cols call, then caches
    _ = em.n_obs
    _ = em.n_obs  # second access — should NOT re-count
    assert counts["cols"] == 1, "n_obs should trigger count_cols exactly once"
    assert counts["rows"] == 0, "n_obs must not trigger count_rows"

    # Accessing n_vars triggers exactly one count_rows call, then caches
    _ = em.n_vars
    _ = em.n_vars  # second access — should NOT re-count
    assert counts["rows"] == 1, "n_vars should trigger count_rows exactly once"
    assert counts["cols"] == 1, "n_vars must not re-trigger count_cols"
