"""Tests for the placeholder ("no data") branches of qtlcascade visualizations.

#308: hvantk.algorithms.qtlcascade.plot rendered its own ad hoc "No data" box in two
places instead of the shared hvantk.algorithms.visualization.base.empty_figure style
that every other plotting domain (enrichex, ptm) uses. These tests pin the shared
style so both placeholder branches stay in sync with the rest of the codebase.
"""

import pandas as pd

from hvantk.algorithms.qtlcascade.plot import plot_attenuation, plot_cascade_classes


def test_empty_cascade_plot_uses_the_shared_placeholder_style():
    """#308: qtlcascade rendered its own 'No data' box while every other domain uses
    visualization.base.empty_figure, so the same empty state looked like two products."""
    # plot_cascade_classes reads a class_counts dict (gene_set_name -> count), not a
    # DataFrame; an empty dict exercises the same "nothing to plot" branch.
    fig = plot_cascade_classes({})
    try:
        ax = fig.axes[0]
        texts = [t.get_text() for t in ax.texts]
        assert "No data available" in texts, texts
        assert not any(spine.get_visible() for spine in ax.spines.values())
    finally:
        import matplotlib.pyplot as plt

        plt.close(fig)


def test_empty_attenuation_plot_uses_the_shared_placeholder_style():
    """#308: plot_attenuation's empty branch inlined the same placeholder idiom as
    plot_cascade_classes; both must render through the shared empty_figure style."""
    df = pd.DataFrame({"cascade_class": []})
    fig = plot_attenuation(df)
    try:
        ax = fig.axes[0]
        texts = [t.get_text() for t in ax.texts]
        assert "No data available" in texts, texts
        assert not any(spine.get_visible() for spine in ax.spines.values())
    finally:
        import matplotlib.pyplot as plt

        plt.close(fig)
