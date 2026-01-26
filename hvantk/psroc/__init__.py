"""
PSROC (Prediction Score ROC Analysis) module.

This module provides functionality for evaluating variant pathogenicity prediction
scores using ROC curve analysis against ClinVar truth labels.

Example:
    >>> from hvantk.psroc import PSROCConfig, PSROCPipeline
    >>> config = PSROCConfig(
    ...     genes=["BRCA1", "BRCA2"],
    ...     clinvar_ht="/data/clinvar.ht",
    ...     dbnsfp_ht="/data/dbnsfp.ht",
    ...     scores=["CADD_phred", "REVEL_score"],
    ...     output_dir="/results/psroc",
    ... )
    >>> pipeline = PSROCPipeline(config)
    >>> result = pipeline.run()
"""

from hvantk.psroc.roc import (
    ScoreMissingness,
    ROCResult,
    compute_roc_metrics,
    find_optimal_threshold,
    compute_score_missingness,
    filter_scores_by_missingness,
    compute_all_missingness,
)

from hvantk.psroc.plots import (
    plot_roc_curves,
    plot_roc_curve_single,
    plot_auc_comparison,
    plot_missingness_summary,
    plot_psroc_summary_dashboard,
)

from hvantk.psroc.pipeline import (
    PSROCConfig,
    PSROCState,
    PSROCResult,
    PSROCPipeline,
    PSROCStage,
    PATHOGENIC_LABELS,
    BENIGN_LABELS,
)

__all__ = [
    # Dataclasses
    "ScoreMissingness",
    "ROCResult",
    # ROC analysis functions
    "compute_roc_metrics",
    "find_optimal_threshold",
    "compute_score_missingness",
    "filter_scores_by_missingness",
    "compute_all_missingness",
    # Plotting functions
    "plot_roc_curves",
    "plot_roc_curve_single",
    "plot_auc_comparison",
    "plot_missingness_summary",
    "plot_psroc_summary_dashboard",
    # Pipeline classes
    "PSROCConfig",
    "PSROCState",
    "PSROCResult",
    "PSROCPipeline",
    "PSROCStage",
    # Constants
    "PATHOGENIC_LABELS",
    "BENIGN_LABELS",
]
