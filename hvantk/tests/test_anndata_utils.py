"""Tests for hvantk.core.models.anndata_utils."""

import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest


def _make_test_adata() -> ad.AnnData:
    """Create a small AnnData for testing (20 obs, 10 var)."""
    np.random.seed(42)
    X = np.random.rand(20, 10)
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(
                np.random.choice(["A", "B", "C"], size=20)
            ),
            "n_counts": np.random.rand(20) * 1000,
        },
        index=[f"cell_{i}" for i in range(20)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(10)])
    return ad.AnnData(X=X, obs=obs, var=var)


class TestAnnotateColumnSummary:
    def test_adds_column_summary_to_uns(self):
        from hvantk.core.models.anndata_utils import annotate_column_summary_ad

        adata = _make_test_adata()
        annotate_column_summary_ad(adata)
        assert "column_summary" in adata.uns

    def test_numeric_column(self):
        from hvantk.core.models.anndata_utils import annotate_column_summary_ad

        adata = _make_test_adata()
        annotate_column_summary_ad(adata)
        summary = adata.uns["column_summary"]
        n_counts = summary["n_counts"]
        assert n_counts["dtype"] == "numeric"
        assert "min" in n_counts
        assert "max" in n_counts
        assert "mean" in n_counts
        assert "n_missing" in n_counts

    def test_categorical_column(self):
        from hvantk.core.models.anndata_utils import annotate_column_summary_ad

        adata = _make_test_adata()
        annotate_column_summary_ad(adata)
        summary = adata.uns["column_summary"]
        ct = summary["cell_type"]
        assert ct["dtype"] == "categorical"
        assert ct["n_unique"] == 3
        assert sorted(ct["levels"]) == ["A", "B", "C"]
        assert "n_missing" in ct

    def test_high_cardinality_categorical(self):
        from hvantk.core.models.anndata_utils import annotate_column_summary_ad

        adata = _make_test_adata()
        # Add a high-cardinality column
        adata.obs["cell_id"] = [f"id_{i}" for i in range(20)]
        annotate_column_summary_ad(adata, max_levels=5, top_n_levels=3)
        summary = adata.uns["column_summary"]
        cid = summary["cell_id"]
        assert cid["dtype"] == "categorical"
        assert cid["n_unique"] == 20
        assert "top_levels" in cid
        assert len(cid["top_levels"]) <= 3
        assert "levels" not in cid


class TestSaveLoadAnndata:
    def test_roundtrip(self, tmp_path):
        from hvantk.core.io.anndata_io import load_anndata, save_anndata

        adata = _make_test_adata()
        path = str(tmp_path / "test.h5ad")
        save_anndata(adata, path)
        loaded = load_anndata(path)
        assert loaded.shape == adata.shape
        np.testing.assert_array_almost_equal(loaded.X, adata.X)

    def test_overwrite_false_raises(self, tmp_path):
        from hvantk.core.io.anndata_io import save_anndata

        adata = _make_test_adata()
        path = str(tmp_path / "test.h5ad")
        save_anndata(adata, path)
        with pytest.raises(FileExistsError):
            save_anndata(adata, path, overwrite=False)

    def test_overwrite_true_succeeds(self, tmp_path):
        from hvantk.core.io.anndata_io import save_anndata

        adata = _make_test_adata()
        path = str(tmp_path / "test.h5ad")
        save_anndata(adata, path)
        save_anndata(adata, path, overwrite=True)
        assert os.path.exists(path)
