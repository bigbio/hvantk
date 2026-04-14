"""Tests for hvantk.core.anndata_utils."""

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


class TestBuildAnndataMetadata:
    def test_returns_required_keys(self):
        from hvantk.core.anndata_utils import build_anndata_metadata

        meta = build_anndata_metadata("ClinVar", "/data/clinvar.vcf")
        required = {
            "hvantk_version",
            "source_name",
            "source_description",
            "raw_input_path",
            "build_date",
        }
        assert required.issubset(meta.keys())

    def test_known_source_ucsc(self):
        from hvantk.core.anndata_utils import build_anndata_metadata

        meta = build_anndata_metadata("UCSC", "/data/ucsc.h5ad")
        assert meta["source_name"] == "UCSC"
        assert "UCSC Cell Browser" in meta["source_description"]

    def test_known_source_expressionatlas(self):
        from hvantk.core.anndata_utils import build_anndata_metadata

        meta = build_anndata_metadata("Expression-Atlas", "/data/ea.h5ad")
        assert "Expression Atlas" in meta["source_description"]

    def test_known_source_cptac(self):
        from hvantk.core.anndata_utils import build_anndata_metadata

        meta = build_anndata_metadata("CPTAC", "/data/cptac.h5ad")
        assert "Proteomics" in meta["source_description"]

    def test_unknown_source(self):
        from hvantk.core.anndata_utils import build_anndata_metadata

        meta = build_anndata_metadata("UnknownDB", "/data/unk.tsv")
        assert meta["source_description"] == ""

    def test_input_path_stored(self):
        from hvantk.core.anndata_utils import build_anndata_metadata

        meta = build_anndata_metadata("test", "/my/path.tsv")
        assert meta["raw_input_path"] == "/my/path.tsv"


class TestAnnotateColumnSummary:
    def test_adds_column_summary_to_uns(self):
        from hvantk.core.anndata_utils import annotate_column_summary_ad

        adata = _make_test_adata()
        annotate_column_summary_ad(adata)
        assert "column_summary" in adata.uns

    def test_numeric_column(self):
        from hvantk.core.anndata_utils import annotate_column_summary_ad

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
        from hvantk.core.anndata_utils import annotate_column_summary_ad

        adata = _make_test_adata()
        annotate_column_summary_ad(adata)
        summary = adata.uns["column_summary"]
        ct = summary["cell_type"]
        assert ct["dtype"] == "categorical"
        assert ct["n_unique"] == 3
        assert sorted(ct["levels"]) == ["A", "B", "C"]
        assert "n_missing" in ct

    def test_high_cardinality_categorical(self):
        from hvantk.core.anndata_utils import annotate_column_summary_ad

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
        from hvantk.core.anndata_utils import load_anndata, save_anndata

        adata = _make_test_adata()
        path = str(tmp_path / "test.h5ad")
        save_anndata(adata, path)
        loaded = load_anndata(path)
        assert loaded.shape == adata.shape
        np.testing.assert_array_almost_equal(loaded.X, adata.X)

    def test_overwrite_false_raises(self, tmp_path):
        from hvantk.core.anndata_utils import save_anndata

        adata = _make_test_adata()
        path = str(tmp_path / "test.h5ad")
        save_anndata(adata, path)
        with pytest.raises(FileExistsError):
            save_anndata(adata, path, overwrite=False)

    def test_overwrite_true_succeeds(self, tmp_path):
        from hvantk.core.anndata_utils import save_anndata

        adata = _make_test_adata()
        path = str(tmp_path / "test.h5ad")
        save_anndata(adata, path)
        save_anndata(adata, path, overwrite=True)
        assert os.path.exists(path)
