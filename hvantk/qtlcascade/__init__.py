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

logger = logging.getLogger(__name__)

# Core cascade
from hvantk.qtlcascade.cascade import build_cascade

# Gene-level summary
from hvantk.qtlcascade.gene_summary import build_cascade_gene_summary

# Colocalization (pure NumPy — always available)
from hvantk.qtlcascade.coloc import coloc_abf, compute_log_abf, run_coloc_per_gene

# Constants
from hvantk.qtlcascade.constants import (
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
    EQTL_SOURCES,
    FANG_TISSUES,
    FANG_TISSUE_EQTL_MAPPING,
    PQTL_SOURCES,
)

# Pipeline
from hvantk.qtlcascade.pipeline import CascadeConfig, CascadePipeline, CascadeResult

# Plotting
from hvantk.qtlcascade.plot import (
    encode_figure_to_base64,
    plot_attenuation,
    plot_cascade_classes,
    plot_coloc_posteriors,
    plot_cross_tissue_heatmap,
    plot_loeuf_by_cascade_class,
)

# Report
from hvantk.qtlcascade.report import generate_report

__all__ = [
    # Cascade core
    "build_cascade",
    "build_cascade_gene_summary",
    # Coloc
    "compute_log_abf",
    "coloc_abf",
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
    "FANG_TISSUES",
    "FANG_TISSUE_EQTL_MAPPING",
    "EQTL_SOURCES",
    "PQTL_SOURCES",
]

__version__ = "0.1.0"
