"""Tests for Phase 4 visualization functions."""
import pandas as pd
import pytest
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for testing

from hvantk.enrichex.plot import (
    plot_celltype_burden_heatmap,
    plot_burden_volcano,
    plot_celltype_forest,
)


def _make_burden_results():
    """Create a realistic combined burden results DataFrame."""
    rows = []
    for collection in ["heart", "brain"]:
        for vc in ["lof", "missense", "synonymous"]:
            for gs in [f"celltype_{i}" for i in range(5)]:
                rows.append({
                    "gene_set_name": gs,
                    "variant_class": vc,
                    "collection": collection,
                    "p_value": 0.001 + len(gs) * 0.01,
                    "p_adjusted": 0.005 + len(gs) * 0.02,
                    "odds_ratio": 1.5 + len(gs) * 0.1,
                    "beta": 0.3 + len(gs) * 0.05,
                    "ci_lower": 1.1,
                    "ci_upper": 2.2,
                    "n_carriers": 15,
                })
    return pd.DataFrame(rows)


class TestCelltypeBurdenHeatmap:
    def test_basic(self, tmp_path):
        df = _make_burden_results()
        fig = plot_celltype_burden_heatmap(df, str(tmp_path / "heatmap.png"))
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_empty_df(self, tmp_path):
        fig = plot_celltype_burden_heatmap(pd.DataFrame(), str(tmp_path / "heatmap.png"))
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_custom_variant_classes(self, tmp_path):
        df = _make_burden_results()
        fig = plot_celltype_burden_heatmap(
            df, str(tmp_path / "heatmap.png"),
            variant_classes=["lof", "missense"]
        )
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)


class TestBurdenVolcano:
    def test_basic(self, tmp_path):
        df = _make_burden_results()
        fig = plot_burden_volcano(df, str(tmp_path / "volcano.png"))
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_empty_df(self, tmp_path):
        fig = plot_burden_volcano(pd.DataFrame(), str(tmp_path / "volcano.png"))
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_color_by_collection(self, tmp_path):
        df = _make_burden_results()
        fig = plot_burden_volcano(
            df, str(tmp_path / "volcano.png"),
            color_by="collection"
        )
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)


class TestCelltypeForest:
    def test_basic(self, tmp_path):
        df = _make_burden_results()
        fig = plot_celltype_forest(
            df, str(tmp_path / "forest.png"),
            cell_type="celltype_1"
        )
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_empty_df(self, tmp_path):
        fig = plot_celltype_forest(
            pd.DataFrame(), str(tmp_path / "forest.png"),
            cell_type="foo"
        )
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_missing_celltype(self, tmp_path):
        df = _make_burden_results()
        fig = plot_celltype_forest(
            df, str(tmp_path / "forest.png"),
            cell_type="nonexistent"
        )
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)
