"""Symbolic <*>/<NON_REF> rows are reference blocks, not variants.

A gVCF reference block densifies into a row like ("T", "<*>") -- the "any other allele"
placeholder DeepVariant and GATK emit. It is not a variant, carries AC=0 on the alt, and
has no business in a variant-QC population.

Measured on the 1005-sample chr20 dense MatrixTable: 5,139,209 of 11,396,989 rows
(45.09%) were <*>. A QC report over that is describing reference blocks as much as
variants, which is why its allele-frequency spectrum looked nothing like a real one.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from hvantk.algorithms.hgc.qc import SYMBOLIC_ALT_ALLELES
from hvantk.algorithms.visualization.qc_plots import (
    _stats_box,
    plot_allele_frequency_spectrum,
)


def test_symbolic_allele_set_covers_the_gvcf_placeholders():
    """Both spellings occur in the wild: DeepVariant emits <*>, GATK <NON_REF>."""
    assert "<*>" in SYMBOLIC_ALT_ALLELES
    assert "<NON_REF>" in SYMBOLIC_ALT_ALLELES
    assert "*" in SYMBOLIC_ALT_ALLELES, "the VCF spanning-deletion allele too"
    # Real alleles must never be caught by it.
    for real in ("A", "C", "G", "T", "AT", "TTTC"):
        assert real not in SYMBOLIC_ALT_ALLELES


@pytest.mark.parametrize(
    "legend_loc,expected_x,expected_ha",
    [("upper right", 0.02, "left"), ("upper left", 0.98, "right")],
)
def test_stats_box_avoids_the_legend_corner(legend_loc, expected_x, expected_ha):
    """The box must sit opposite the legend.

    plot_allele_frequency_spectrum pinned the stats box at "upper right" while also
    calling ax.legend(loc="upper right"), so "Mean: 0.2071" rendered on top of
    "5% MAF". plot_variant_qc_overview had the same collision whenever log_transform
    was False.
    """
    fig, ax = plt.subplots()
    _stats_box(ax, "Mean: 0.5", legend_loc=legend_loc)

    texts = [t for t in ax.texts if "Mean" in t.get_text()]
    assert len(texts) == 1
    t = texts[0]
    assert t.get_position()[0] == pytest.approx(expected_x)
    assert t.get_horizontalalignment() == expected_ha
    plt.close(fig)


def test_af_spectrum_stats_box_and_legend_do_not_overlap():
    """End-to-end on the real plotting function, in figure coordinates."""
    rng = np.random.default_rng(0)
    df = pd.DataFrame({"AF": np.clip(rng.beta(0.4, 3.0, 5000), 0, 1)})

    fig = plot_allele_frequency_spectrum(df, figsize=(11, 5.6))
    ax = fig.axes[0]
    fig.canvas.draw()  # positions are only real after a draw

    stats = [t for t in ax.texts if "Mean:" in t.get_text()]
    assert stats, "stats box missing"
    stats_bb = stats[0].get_window_extent(fig.canvas.get_renderer())

    legend = ax.get_legend()
    assert legend is not None, "legend missing"
    legend_bb = legend.get_window_extent(fig.canvas.get_renderer())

    assert not stats_bb.overlaps(legend_bb), (
        f"stats box {stats_bb} overlaps legend {legend_bb} -- they are back in the "
        f"same corner"
    )
    plt.close(fig)
