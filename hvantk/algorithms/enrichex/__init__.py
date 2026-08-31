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

logger = logging.getLogger(__name__)

try:
    from hvantk.algorithms.enrichex.burden import (
        VariantFilter,
        build_variant_classes_from_presets,
        permutation_burden_test,
        compute_geneset_burden_mt,
        compute_per_gene_burden_mt,
        linear_burden_test,
        logistic_burden_test,
        run_burden_analysis,
        run_stratified_burden_analysis,
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

    class VariantFilter:
        def __init__(self, *a, **k):
            raise ImportError(
                "Burden analysis requires the optional Hail dependency. "
                "Install hvantk with the 'hail' requirements."
            ) from _exc

    build_variant_classes_from_presets = VariantFilter  # type: ignore[assignment]
    permutation_burden_test = VariantFilter  # type: ignore[assignment]
    compute_geneset_burden_mt = VariantFilter  # type: ignore[assignment]
    compute_per_gene_burden_mt = VariantFilter  # type: ignore[assignment]
    linear_burden_test = VariantFilter  # type: ignore[assignment]
    logistic_burden_test = VariantFilter  # type: ignore[assignment]
    run_burden_analysis = VariantFilter  # type: ignore[assignment]
    run_stratified_burden_analysis = VariantFilter  # type: ignore[assignment]
from hvantk.algorithms.enrichex.constants import (
    CORRECTION_METHODS,
    DEFAULT_ALPHA,
    DEFAULT_CORRECTION_METHOD,
    DEFAULT_MAX_AF,
    DEFAULT_MIN_SCORE,
    DEFAULT_N_PERMUTATIONS,
    GENOTYPE_AGGREGATION_METHODS,
    PHENOTYPE_TYPES,
    VARIANT_CLASS_PRESETS,
)
from hvantk.algorithms.statistics.correction import apply_correction, fdr_threshold
from hvantk.core.utils.gene_sets import (
    GeneSet,
    GeneSetCollection,
    load_gene_sets_from_dict,
    load_marker_genes,
)
from hvantk.algorithms.enrichex.overlap import (
    OverlapResult,
    compute_overlap_enrichment,
    compute_overlap_enrichment_pandas,
)
from hvantk.algorithms.enrichex.pipeline import (
    BurdenConfig,
    BurdenPipeline,
    BurdenRunResult,
)
from hvantk.algorithms.enrichex.simulation import check_type_i_error

# Plotting and reporting are resolved on ATTRIBUTE ACCESS, not at package import (PEP 562).
#
# They used to be plain `from ... import` lines here, and that broke `hvantk` outright: the
# CLI imports hvantk.algorithms.enrichex.constants, which runs this module, which pulled in
# enrichex.plot / enrichex.report / visualization.base -- all three import matplotlib at
# module scope. matplotlib is optional, so on a base install `pip install hvantk` produced a
# console script where even `hvantk --help` raised ModuleNotFoundError. Every command paid
# for plotting whether or not it plotted.
#
# The names below stay importable and stay in __all__, so this is not an API change; the
# import simply happens on first use. The CLI never trips it -- overlap_cli and burden_cli
# already import `generate_report` inside the functions that need it -- so a base install
# now runs, and the `enrichex` extra (which carries matplotlib) is what a plotting caller
# installs.
_LAZY = {
    "encode_figure_to_base64": "hvantk.algorithms.visualization.base",
    "generate_report": "hvantk.algorithms.enrichex.report",
    "plot_burden_forest": "hvantk.algorithms.enrichex.plot",
    "plot_burden_volcano": "hvantk.algorithms.enrichex.plot",
    "plot_celltype_burden_heatmap": "hvantk.algorithms.enrichex.plot",
    "plot_celltype_forest": "hvantk.algorithms.enrichex.plot",
    "plot_enrichment_barplot": "hvantk.algorithms.enrichex.plot",
    "plot_enrichment_dotplot": "hvantk.algorithms.enrichex.plot",
}


_PLOTTING_HINT = (
    "matplotlib is required for EnrichEx plotting and reporting, and is not part of the "
    "base install. Install the 'enrichex' extra:\n"
    "    pip install 'hvantk[enrichex]'\n"
    "    poetry install --extras enrichex"
)


def __getattr__(name: str):
    """Resolve the plotting/reporting exports on first access."""
    module = _LAZY.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib

    try:
        value = getattr(importlib.import_module(module), name)
    except ModuleNotFoundError as exc:  # pragma: no cover - needs matplotlib absent
        # Same contract as require_scanpy (algorithms/expression/matrix_utils.py) and
        # _require_matplotlib (algorithms/qtlcascade/plot.py): name the extra rather than
        # leave the caller a bare ModuleNotFoundError. Deferring the import must not also
        # degrade the message -- that was the review finding on the scanpy extra in #248.
        if exc.name and exc.name.split(".")[0] in {"matplotlib", "seaborn"}:
            raise ImportError(_PLOTTING_HINT) from exc
        raise
    globals()[name] = value  # cache, so the import cost is paid once
    return value


def __dir__():
    return sorted(set(globals()) | set(_LAZY))

try:
    from hvantk.algorithms.enrichex.simulation import generate_synthetic_burden_cohort
except ModuleNotFoundError:  # pragma: no cover - depends on optional Hail install
    generate_synthetic_burden_cohort = VariantFilter  # type: ignore[assignment]

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
    "build_variant_classes_from_presets",
    "permutation_burden_test",
    "compute_geneset_burden_mt",
    "compute_per_gene_burden_mt",
    "logistic_burden_test",
    "linear_burden_test",
    "run_burden_analysis",
    "run_stratified_burden_analysis",
    # Multiple testing correction
    "apply_correction",
    "fdr_threshold",
    # Plotting
    "plot_enrichment_dotplot",
    "plot_enrichment_barplot",
    "plot_burden_forest",
    "plot_burden_volcano",
    "plot_celltype_burden_heatmap",
    "plot_celltype_forest",
    "encode_figure_to_base64",
    # Pipeline orchestration
    "BurdenConfig",
    "BurdenPipeline",
    "BurdenRunResult",
    # Reporting
    "generate_report",
    # Simulation / Validation
    "generate_synthetic_burden_cohort",
    "check_type_i_error",
    # Constants
    "DEFAULT_ALPHA",
    "DEFAULT_CORRECTION_METHOD",
    "DEFAULT_MAX_AF",
    "DEFAULT_MIN_SCORE",
    "DEFAULT_N_PERMUTATIONS",
    "CORRECTION_METHODS",
    "GENOTYPE_AGGREGATION_METHODS",
    "PHENOTYPE_TYPES",
    "VARIANT_CLASS_PRESETS",
]

__version__ = "0.1.0"
