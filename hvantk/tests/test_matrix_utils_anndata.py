"""Tests for AnnData-based expression analysis functions in matrix_utils."""

import numpy as np
import pandas as pd
import pytest

import anndata as ad

from hvantk.core.anndata_utils import annotate_column_summary_ad
import importlib.util
import sys
from pathlib import Path

# Load matrix_utils directly to sidestep hvantk.utils.__init__ which imports
# modules using Python 3.10+ syntax (e.g. ``int | None`` in correction.py).
_mu_path = str(Path(__file__).resolve().parent.parent / "utils" / "matrix_utils.py")
_spec = importlib.util.spec_from_file_location("_matrix_utils_direct", _mu_path)
_mod = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _mod
_spec.loader.exec_module(_mod)

describe_expression_ad = _mod.describe_expression_ad
filter_by_metadata_ad = _mod.filter_by_metadata_ad
summarize_expression_ad = _mod.summarize_expression_ad


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
    def test_returns_dataframe(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        assert isinstance(result, pd.DataFrame)

    def test_has_expected_columns(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        expected_cols = {"gene_id", "group", "mean", "fraction_expressed", "n_cells"}
        assert expected_cols == set(result.columns)

    def test_groups_match_unique_cell_types(self, test_adata):
        result = summarize_expression_ad(test_adata, group_by="cell_type")
        assert set(result["group"].unique()) == {"neuron", "astrocyte", "microglia"}
