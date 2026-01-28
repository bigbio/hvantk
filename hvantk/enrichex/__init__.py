"""
EnrichEx: Gene Set Enrichment Analysis

This module provides tools for gene set enrichment analysis, including:
- Gene set data structures and I/O
- Overlap enrichment testing (Fisher's exact test)
- Burden testing with Hail-native regression
- Multiple testing correction

Main entry points:
- GeneSet: Named set of genes
- GeneSetCollection: Collection of gene sets with background
- compute_overlap_enrichment: Fisher's exact test for overlap
- run_burden_analysis: Case-control burden testing

CLI commands:
- hvantk enrichex overlap: Test gene list enrichment
- hvantk enrichex burden: Case-control burden analysis
"""

from __future__ import annotations

import logging
from typing import Any

logger = logging.getLogger(__name__)

try:
    from hvantk.enrichex.burden import (
        VariantFilter,
        compute_geneset_burden_mt,
        linear_burden_test,
        logistic_burden_test,
        run_burden_analysis,
    )
except (
    ModuleNotFoundError
) as exc:  # pragma: no cover - depends on optional Hail install
    _exc = exc
    logger.warning(
        "Hail-dependent EnrichEx burden methods unavailable: %s. "
        "Install hvantk with the required extras to enable them.",
        exc,
    )

    def _missing_dependency(*_: Any, **__: Any) -> None:
        raise ImportError(
            "Burden analysis requires the optional Hail dependency. "
            "Install hvantk with the 'hail' requirements."
        ) from _exc

    VariantFilter = None  # type: ignore[assignment]
    compute_geneset_burden_mt = _missing_dependency  # type: ignore[assignment]
    linear_burden_test = _missing_dependency  # type: ignore[assignment]
    logistic_burden_test = _missing_dependency  # type: ignore[assignment]
    run_burden_analysis = _missing_dependency  # type: ignore[assignment]
from hvantk.enrichex.constants import (
    CORRECTION_METHODS,
    DEFAULT_ALPHA,
    DEFAULT_CORRECTION_METHOD,
    DEFAULT_MAX_AF,
    DEFAULT_MIN_CADD,
    GENOTYPE_AGGREGATION_METHODS,
    PHENOTYPE_TYPES,
)
from hvantk.enrichex.correction import apply_correction, fdr_threshold
from hvantk.enrichex.gene_sets import (
    GeneSet,
    GeneSetCollection,
    load_gene_sets_from_dict,
    load_marker_genes,
)
from hvantk.enrichex.plot import (
    encode_figure_to_base64,
    plot_burden_forest,
    plot_enrichment_barplot,
    plot_enrichment_dotplot,
)
from hvantk.enrichex.overlap import (
    OverlapResult,
    compute_overlap_enrichment,
    compute_overlap_enrichment_pandas,
)
from hvantk.enrichex.report import generate_report

__all__ = [
    # Core data structures
    "GeneSet",
    "GeneSetCollection",
    # Loading functions
    "load_gene_sets_from_dict",
    "load_marker_genes",
    # Overlap enrichment
    "OverlapResult",
    "compute_overlap_enrichment",
    "compute_overlap_enrichment_pandas",
    # Burden testing
    "VariantFilter",
    "compute_geneset_burden_mt",
    "logistic_burden_test",
    "linear_burden_test",
    "run_burden_analysis",
    # Multiple testing correction
    "apply_correction",
    "fdr_threshold",
    # Plotting
    "plot_enrichment_dotplot",
    "plot_enrichment_barplot",
    "plot_burden_forest",
    "encode_figure_to_base64",
    # Reporting
    "generate_report",
    # Constants
    "DEFAULT_ALPHA",
    "DEFAULT_CORRECTION_METHOD",
    "DEFAULT_MAX_AF",
    "DEFAULT_MIN_CADD",
    "CORRECTION_METHODS",
    "GENOTYPE_AGGREGATION_METHODS",
    "PHENOTYPE_TYPES",
]

__version__ = "0.1.0"
