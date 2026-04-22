"""Tests for AnnData-based expression analysis functions in matrix_utils."""

import numpy as np
import pandas as pd
import pytest

import anndata as ad

from hvantk.core.anndata_utils import annotate_column_summary_ad
from hvantk.utils.matrix_utils import (
    describe_expression_ad,
    filter_by_metadata_ad,
    summarize_expression_ad,
)


@pytest.fixture
def test_adata():
    """Create a test AnnData object with 100 obs, 50 var."""
    np.random.seed(42)
    X = np.random.rand(100, 50).astype(np.float32)
    # Sprinkle some zeros for fraction_expressed testing
    X[X < 0.2] = 0.0

    obs = pd.DataFrame(
        {
            "cell_type": np.random.choice(
                ["neuron", "astrocyte", "microglia"], size=100
            ),
            "tissue": np.random.choice(["brain", "spinal_cord"], size=100),
        },
        index=[f"cell_{i}" for i in range(100)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(50)])

    adata = ad.AnnData(X=X, obs=obs, var=var)
    annotate_column_summary_ad(adata)
    return adata


class TestDescribeExpressionAd:
    def test_returns_summary_dict(self, test_adata):
        result = describe_expression_ad(test_adata)
        assert isinstance(result, dict)
        assert result["n_obs"] == 100
        assert result["n_vars"] == 50
        assert any(f["name"] == "cell_type" for f in result["fields"])


class TestFilterByMetadataAd:
    def test_filters_by_single_value(self, test_adata):
        filtered = filter_by_metadata_ad(test_adata, {"cell_type": "neuron"})
        assert all(filtered.obs["cell_type"] == "neuron")
        assert filtered.n_obs > 0

    def test_filters_by_list(self, test_adata):
        filtered = filter_by_metadata_ad(
            test_adata, {"cell_type": ["neuron", "astrocyte"]}
        )
        assert set(filtered.obs["cell_type"].unique()).issubset(
            {"neuron", "astrocyte"}
        )
        assert filtered.n_obs > 0


class TestSummarizeExpressionAd:
    def test_returns_anndata(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        assert isinstance(result, ad.AnnData)

    def test_shape_is_groups_by_genes(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        expected_groups = len(set(test_adata.obs["cell_type"]))
        assert result.shape == (expected_groups, test_adata.n_vars)

    def test_has_required_layers(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        expected_layers = {"mean", "sum", "count_nonzero", "fraction_expressed"}
        assert expected_layers.issubset(set(result.layers.keys()))

    def test_obs_has_n_cells(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        assert "n_cells" in result.obs.columns
        assert (result.obs["n_cells"] > 0).all()
        assert int(result.obs["n_cells"].sum()) == test_adata.n_obs

    def test_groups_match_unique_cell_types(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        assert set(result.obs_names) == {"neuron", "astrocyte", "microglia"}

    def test_fraction_expressed_between_0_and_1(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        frac = np.asarray(result.layers["fraction_expressed"])
        assert frac.min() >= 0.0
        assert frac.max() <= 1.0

    def test_min_cells_filters_small_groups(self, test_adata):
        # all three cell types should be well above n=5 in the fixture (100 cells / 3),
        # but at min_cells=10_000 we should drop every group.
        result = summarize_expression_ad(
            test_adata, group_by="cell_type", min_cells_per_group=10_000
        )
        assert result.n_obs == 0
