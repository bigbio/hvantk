"""
QTL Cascade: Molecular QTL Cascade Analysis
============================================

Traces variant effects across omics layers:

    Variant → eQTL (mRNA) → pQTL (protein) → Disease

Main entry points:

* :func:`build_cascade` — outer-join eQTL + pQTL tables, classify pairs
* :func:`build_cascade_gene_summary` — gene-level aggregation + overlays
* :func:`coloc_abf` — colocalization ABF for a single region
* :func:`run_coloc_per_gene` — coloc across cascade genes (Hail + NumPy)
* :class:`CascadePipeline` — stage-based pipeline orchestration
* :func:`generate_report` — static HTML report

CLI commands::

    hvantk qtlcascade cascade   — build the cascade join
    hvantk qtlcascade coloc     — run colocalization
    hvantk qtlcascade run       — full pipeline
    hvantk qtlcascade report    — generate HTML report
"""

from __future__ import annotations

import logging

from hvantk.core.utils.lazy_exports import install_lazy_exports

logger = logging.getLogger(__name__)

# Constants
from hvantk.algorithms.qtlcascade.constants import (
    CASCADE_CLASSES,
    CASCADE_CLASS_COLORS,
    CASCADE_CLASS_LABELS,
    DEFAULT_COLOC_H4_THRESHOLD,
    DEFAULT_COLOC_P1,
    DEFAULT_COLOC_P12,
    DEFAULT_COLOC_P2,
    DEFAULT_COLOC_W,
    DEFAULT_COLOC_WINDOW_KB,
    DEFAULT_EQTL_P_THRESHOLD,
    DEFAULT_PQTL_P_THRESHOLD,
)


def _missing_hint(exc: ModuleNotFoundError) -> "str | None":
    if (exc.name or "").split(".")[0] in {"matplotlib", "seaborn"}:
        return (
            "matplotlib is required for qtlcascade plotting and reporting, and is not "
            "part of the base install. Install the 'viz' extra:\n"
            "    pip install 'hvantk[viz]'\n"
            "    poetry install --extras viz"
        )
    return None


# Everything below is resolved on ATTRIBUTE ACCESS, not at package import (PEP 562).
#
# The qtlcascade CLI needs constants (e.g. CASCADE_CLASSES) at decorator time, so it
# imports hvantk.algorithms.qtlcascade.constants, which runs this module -- and this
# module eagerly imported cascade, gene_summary (pandas), coloc, pipeline, plot
# (matplotlib) and report, whether or not the command being run needed any of them.
# `hvantk qtlcascade --help` paid the full cost just to print its own help text (#306).
#
# The names below stay importable and stay in __all__, so this is not an API change; the
# import simply happens on first use, via install_lazy_exports (hvantk/core/utils/lazy_exports.py).
install_lazy_exports(
    globals(),
    {
        "build_cascade": "hvantk.algorithms.qtlcascade.cascade",
        "build_cascade_gene_summary": "hvantk.algorithms.qtlcascade.gene_summary",
        "coloc_abf": "hvantk.algorithms.qtlcascade.coloc",
        "compute_log_abf": "hvantk.algorithms.qtlcascade.coloc",
        "prepare_coloc_data": "hvantk.algorithms.qtlcascade.coloc",
        "run_coloc_per_gene": "hvantk.algorithms.qtlcascade.coloc",
        "CascadeConfig": "hvantk.algorithms.qtlcascade.pipeline",
        "CascadePipeline": "hvantk.algorithms.qtlcascade.pipeline",
        "CascadeResult": "hvantk.algorithms.qtlcascade.pipeline",
        "encode_figure_to_base64": "hvantk.algorithms.qtlcascade.plot",
        "plot_attenuation": "hvantk.algorithms.qtlcascade.plot",
        "plot_cascade_classes": "hvantk.algorithms.qtlcascade.plot",
        "plot_coloc_posteriors": "hvantk.algorithms.qtlcascade.plot",
        "plot_cross_tissue_heatmap": "hvantk.algorithms.qtlcascade.plot",
        "plot_loeuf_by_cascade_class": "hvantk.algorithms.qtlcascade.plot",
        "generate_report": "hvantk.algorithms.qtlcascade.report",
    },
    missing_hint=_missing_hint,
)

__all__ = [
    # Cascade core
    "build_cascade",
    "build_cascade_gene_summary",
    # Coloc
    "compute_log_abf",
    "coloc_abf",
    "prepare_coloc_data",
    "run_coloc_per_gene",
    # Pipeline
    "CascadeConfig",
    "CascadePipeline",
    "CascadeResult",
    # Plots
    "plot_cascade_classes",
    "plot_attenuation",
    "plot_coloc_posteriors",
    "plot_cross_tissue_heatmap",
    "plot_loeuf_by_cascade_class",
    "encode_figure_to_base64",
    # Report
    "generate_report",
    # Constants
    "CASCADE_CLASSES",
    "CASCADE_CLASS_COLORS",
    "CASCADE_CLASS_LABELS",
    "DEFAULT_EQTL_P_THRESHOLD",
    "DEFAULT_PQTL_P_THRESHOLD",
    "DEFAULT_COLOC_P1",
    "DEFAULT_COLOC_P2",
    "DEFAULT_COLOC_P12",
    "DEFAULT_COLOC_W",
    "DEFAULT_COLOC_H4_THRESHOLD",
    "DEFAULT_COLOC_WINDOW_KB",
]

__version__ = "0.1.0"
