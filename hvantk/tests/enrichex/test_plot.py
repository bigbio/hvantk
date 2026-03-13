"""Tests for enrichex visualization functions (core + phase 4)."""

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest

matplotlib.use("Agg")

from hvantk.enrichex.plot import (
    encode_figure_to_base64,
    plot_burden_forest,
    plot_burden_volcano,
    plot_celltype_burden_heatmap,
    plot_celltype_forest,
    plot_enrichment_barplot,
    plot_enrichment_dotplot,
)


def _mock_enrichment_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_set_name": ["Microglia", "Astrocytes", "Neurons", "Endothelium"],
            "n_overlap": [12, 8, 6, 5],
            "odds_ratio": [2.5, 1.8, 1.2, 0.9],
            "p_value": [1e-6, 5e-4, 0.02, 0.2],
            "p_adjusted": [5e-6, 1e-3, 0.05, 0.5],
            "significant": [True, True, False, False],
            "overlap_genes": [
                "APOE,TREM2,CD33",
                "GFAP,S100B",
                "NRGN,GRIN2A",
                "PECAM1,CLDN5",
            ],
        }
    )


def _mock_burden_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_set_name": ["Microglia", "Astrocytes", "Neurons"],
            "odds_ratio": [1.6, 1.2, 0.95],
            "ci_lower": [1.2, 0.9, 0.7],
            "ci_upper": [2.1, 1.5, 1.3],
            "p_value": [0.001, 0.08, 0.6],
            "p_adjusted": [0.002, 0.12, 0.9],
            "significant": [True, False, False],
        }
    )


def _make_burden_results():
    """Create a realistic combined burden results DataFrame for phase 4 tests."""
    rows = []
    for collection in ["heart", "brain"]:
        for vc in ["lof", "missense", "synonymous"]:
            for gs in [f"celltype_{i}" for i in range(5)]:
                rows.append(
                    {
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
                    }
                )
    return pd.DataFrame(rows)


# --- Core plot tests ---


def test_plot_enrichment_dotplot(tmp_path):
    fig = plot_enrichment_dotplot(
        _mock_enrichment_df(), output_path=str(tmp_path / "dotplot.png"), top_n=3, label_top_n=2
    )
    assert (tmp_path / "dotplot.png").exists()
    assert fig.get_axes()
    plt.close(fig)


def test_plot_burden_forest(tmp_path):
    fig = plot_burden_forest(
        _mock_burden_df(),
        output_path=str(tmp_path / "forest.png"),
        phenotype_type="binary",
        top_n=2,
        color_by="significant",
        show_values=True,
    )
    assert (tmp_path / "forest.png").exists()
    assert fig.get_axes()
    plt.close(fig)


def test_plot_enrichment_barplot(tmp_path):
    fig = plot_enrichment_barplot(
        _mock_enrichment_df(),
        output_path=str(tmp_path / "bar.png"),
        value="-log10_p",
        top_n=2,
        orientation="vertical",
    )
    assert (tmp_path / "bar.png").exists()
    assert fig.get_axes()
    plt.close(fig)


def test_encode_figure_to_base64_returns_string():
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    encoded = encode_figure_to_base64(fig)
    assert isinstance(encoded, str)
    assert encoded.startswith(("iVBOR", "/9j/"))
    plt.close(fig)


# --- Phase 4 plot tests ---


def test_celltype_burden_heatmap(tmp_path):
    fig = plot_celltype_burden_heatmap(_make_burden_results(), str(tmp_path / "heatmap.png"))
    assert fig is not None
    plt.close(fig)


def test_celltype_burden_heatmap_custom_variant_classes(tmp_path):
    fig = plot_celltype_burden_heatmap(
        _make_burden_results(),
        str(tmp_path / "heatmap.png"),
        variant_classes=["lof", "missense"],
    )
    assert fig is not None
    plt.close(fig)


def test_burden_volcano(tmp_path):
    fig = plot_burden_volcano(_make_burden_results(), str(tmp_path / "volcano.png"))
    assert fig is not None
    plt.close(fig)


def test_burden_volcano_color_by_collection(tmp_path):
    fig = plot_burden_volcano(
        _make_burden_results(), str(tmp_path / "volcano.png"), color_by="collection"
    )
    assert fig is not None
    plt.close(fig)


def test_celltype_forest(tmp_path):
    fig = plot_celltype_forest(
        _make_burden_results(), str(tmp_path / "forest.png"), cell_type="celltype_1"
    )
    assert fig is not None
    plt.close(fig)


def test_celltype_forest_missing_celltype(tmp_path):
    fig = plot_celltype_forest(
        _make_burden_results(), str(tmp_path / "forest.png"), cell_type="nonexistent"
    )
    assert fig is not None
    plt.close(fig)


# --- Empty DataFrame placeholder tests (all plot functions) ---


@pytest.mark.parametrize(
    "plot_fn,kwargs",
    [
        (plot_enrichment_dotplot, {}),
        (plot_enrichment_barplot, {}),
        (plot_burden_forest, {}),
        (plot_celltype_burden_heatmap, {}),
        (plot_burden_volcano, {}),
        (plot_celltype_forest, {"cell_type": "foo"}),
    ],
    ids=[
        "dotplot",
        "barplot",
        "forest",
        "heatmap",
        "volcano",
        "celltype_forest",
    ],
)
def test_empty_df_returns_placeholder(tmp_path, plot_fn, kwargs):
    """All plot functions handle empty DataFrames gracefully."""
    output_path = tmp_path / "empty.png"
    # burden_forest needs specific columns to avoid KeyError
    if plot_fn == plot_burden_forest:
        df = pd.DataFrame(
            columns=["gene_set_name", "odds_ratio", "ci_lower", "ci_upper", "p_value"]
        )
    else:
        df = pd.DataFrame()
    fig = plot_fn(df, str(output_path), **kwargs)
    assert fig is not None
    plt.close(fig)
