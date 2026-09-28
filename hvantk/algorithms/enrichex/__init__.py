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

from hvantk.core.utils.lazy_exports import install_lazy_exports

logger = logging.getLogger(__name__)

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

_PLOTTING_HINT = (
    "matplotlib is required for EnrichEx plotting and reporting, and is not part of the "
    "base install. Install the 'enrichex' extra:\n"
    "    pip install 'hvantk[enrichex]'\n"
    "    poetry install --extras enrichex"
)

_HAIL_HINT = (
    "Burden analysis requires Hail, which is part of the base install; if this import "
    "fails the environment is broken -- reinstall hvantk (pip install hvantk) or check "
    "that Java 8/11 is available to Hail."
)


def _missing_hint(exc: ModuleNotFoundError) -> "str | None":
    root = (exc.name or "").split(".")[0]
    if root in {"matplotlib", "seaborn"}:
        return _PLOTTING_HINT
    if root == "hail":
        return _HAIL_HINT
    return None


# Everything below is resolved on ATTRIBUTE ACCESS, not at package import (PEP 562).
#
# It used to be plain `from ... import` lines here, and that broke `hvantk` outright: the
# burden CLI needs `GENOTYPE_AGGREGATION_METHODS` at decorator time (`click.Choice(...)`),
# so it imports `hvantk.algorithms.enrichex.constants`, which runs this module -- and this
# module eagerly imported burden (Hail, ~5 s), overlap (scipy), pipeline, simulation and
# gene_sets, whether or not the command being run needed any of them. `hvantk enrichex
# --help` paid the full cost just to print its own help text (#306).
#
# The names below stay importable and stay in __all__, so this is not an API change; the
# import simply happens on first use, via install_lazy_exports (hvantk/core/utils/lazy_exports.py).
install_lazy_exports(
    globals(),
    {
        # Core data structures + loaders
        "GeneSet": "hvantk.core.utils.gene_sets",
        "GeneSetCollection": "hvantk.core.utils.gene_sets",
        "load_gene_sets_from_dict": "hvantk.core.utils.gene_sets",
        "load_marker_genes": "hvantk.core.utils.gene_sets",
        # Overlap enrichment (scipy at module scope in overlap.py)
        "OverlapResult": "hvantk.algorithms.enrichex.overlap",
        "compute_overlap_enrichment": "hvantk.algorithms.enrichex.overlap",
        "compute_overlap_enrichment_pandas": "hvantk.algorithms.enrichex.overlap",
        # Burden testing (Hail)
        "VariantFilter": "hvantk.algorithms.enrichex.burden",
        "build_variant_classes_from_presets": "hvantk.algorithms.enrichex.burden",
        "permutation_burden_test": "hvantk.algorithms.enrichex.burden",
        "compute_geneset_burden_mt": "hvantk.algorithms.enrichex.burden",
        "compute_per_gene_burden_mt": "hvantk.algorithms.enrichex.burden",
        "logistic_burden_test": "hvantk.algorithms.enrichex.burden",
        "linear_burden_test": "hvantk.algorithms.enrichex.burden",
        "run_burden_analysis": "hvantk.algorithms.enrichex.burden",
        "run_stratified_burden_analysis": "hvantk.algorithms.enrichex.burden",
        # Multiple testing correction
        "apply_correction": "hvantk.algorithms.statistics.correction",
        "fdr_threshold": "hvantk.algorithms.statistics.correction",
        # Plotting / reporting (matplotlib)
        "plot_enrichment_dotplot": "hvantk.algorithms.enrichex.plot",
        "plot_enrichment_barplot": "hvantk.algorithms.enrichex.plot",
        "plot_burden_forest": "hvantk.algorithms.enrichex.plot",
        "plot_burden_volcano": "hvantk.algorithms.enrichex.plot",
        "plot_celltype_burden_heatmap": "hvantk.algorithms.enrichex.plot",
        "plot_celltype_forest": "hvantk.algorithms.enrichex.plot",
        "encode_figure_to_base64": "hvantk.algorithms.visualization.base",
        "generate_report": "hvantk.algorithms.enrichex.report",
        # Pipeline orchestration
        "BurdenConfig": "hvantk.algorithms.enrichex.pipeline",
        "BurdenPipeline": "hvantk.algorithms.enrichex.pipeline",
        "BurdenRunResult": "hvantk.algorithms.enrichex.pipeline",
        # Simulation / validation
        "generate_synthetic_burden_cohort": "hvantk.algorithms.enrichex.simulation",
        "check_type_i_error": "hvantk.algorithms.enrichex.simulation",
    },
    missing_hint=_missing_hint,
)

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
