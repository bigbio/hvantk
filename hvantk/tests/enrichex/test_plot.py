import matplotlib.pyplot as plt
import pandas as pd
import pytest

from hvantk.enrichex.plot import (
    encode_figure_to_base64,
    plot_burden_forest,
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


def test_plot_enrichment_dotplot(tmp_path):
    df = _mock_enrichment_df()
    output_path = tmp_path / "dotplot.png"
    fig = plot_enrichment_dotplot(
        df,
        output_path=str(output_path),
        top_n=3,
        label_top_n=2,
    )
    assert output_path.exists()
    assert fig.get_axes()
    plt.close(fig)


def test_plot_burden_forest(tmp_path):
    df = _mock_burden_df()
    output_path = tmp_path / "forest.png"
    fig = plot_burden_forest(
        df,
        output_path=str(output_path),
        phenotype_type="binary",
        top_n=2,
        color_by="significant",
        show_values=True,
    )
    assert output_path.exists()
    assert fig.get_axes()
    plt.close(fig)


def test_plot_enrichment_barplot(tmp_path):
    df = _mock_enrichment_df()
    output_path = tmp_path / "bar.png"
    fig = plot_enrichment_barplot(
        df,
        output_path=str(output_path),
        value="-log10_p",
        top_n=2,
        orientation="vertical",
    )
    assert output_path.exists()
    assert fig.get_axes()
    plt.close(fig)


def test_encode_figure_to_base64_returns_string():
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    encoded = encode_figure_to_base64(fig)
    assert isinstance(encoded, str)
    assert encoded.startswith(("iVBOR", "/9j/"))
    plt.close(fig)


def test_plot_enrichment_dotplot_empty_df_raises():
    with pytest.raises(ValueError):
        plot_enrichment_dotplot(pd.DataFrame(), output_path="out.png")
