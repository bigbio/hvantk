"""
Ancestry inference module for hvantk.

This module provides tools for inferring genetic ancestry using PCA-based
methods with Random Forest classification.

Main components:
- merge: MatrixTable merging utilities
- filter: Variant filtering for ancestry analysis
- pca: PCA computation and projection
- classify: Random Forest training and prediction
- pipeline: End-to-end ancestry inference orchestration

Example usage:
    >>> from hvantk.ancestry import run_ancestry_inference
    >>> result = run_ancestry_inference(
    ...     query_mt=hl.read_matrix_table("cohort.mt"),
    ...     reference_mt=hl.read_matrix_table("1kg.mt"),
    ...     ancestry_col="super_pop",
    ... )
    >>> predictions = result.get_predictions_df()
"""

from hvantk.ancestry.constants import (
    # Pipeline defaults
    DEFAULT_MIN_AF,
    DEFAULT_MAX_AF,
    DEFAULT_MIN_CALL_RATE,
    DEFAULT_LD_R2,
    DEFAULT_LD_WINDOW,
    DEFAULT_N_PCS,
    DEFAULT_N_PCS_CLASSIFY,
    DEFAULT_N_ESTIMATORS,
    DEFAULT_MIN_PROB,
    # Column names
    ANCESTRY_COL,
    PREDICTED_ANCESTRY_COL,
    # Population mappings
    SUPERPOP_TO_SUBPOPS,
    SUBPOP_TO_SUPERPOP,
    ALL_1KG_SUPERPOPS,
    ALL_1KG_SUBPOPS,
    POPULATION_NAMES,
    # Colors
    SUPERPOP_COLORS,
    POPULATION_COLORS,
)

from hvantk.ancestry.merge import (
    validate_matrixtable_compatibility,
    get_shared_variants_stats,
    merge_matrixtables,
    check_sample_overlap,
)

from hvantk.ancestry.filter import (
    filter_to_autosomes,
    filter_to_biallelic_snps,
    compute_variant_qc,
    filter_variants_for_ancestry,
    ld_prune_variants,
    prepare_ancestry_variants,
)

from hvantk.ancestry.pca import (
    PCAResult,
    compute_pca,
    project_samples,
)

from hvantk.ancestry.classify import (
    ClassificationResult,
    train_classifier,
    predict_ancestry,
    validate_training_data,
    get_training_samples,
    get_query_samples,
)

from hvantk.ancestry.pipeline import (
    PipelineConfig,
    AncestryInferenceResult,
    run_ancestry_inference,
)

__all__ = [
    # Pipeline defaults
    "DEFAULT_MIN_AF",
    "DEFAULT_MAX_AF",
    "DEFAULT_MIN_CALL_RATE",
    "DEFAULT_LD_R2",
    "DEFAULT_LD_WINDOW",
    "DEFAULT_N_PCS",
    "DEFAULT_N_PCS_CLASSIFY",
    "DEFAULT_N_ESTIMATORS",
    "DEFAULT_MIN_PROB",
    # Column names
    "ANCESTRY_COL",
    "PREDICTED_ANCESTRY_COL",
    # Population mappings
    "SUPERPOP_TO_SUBPOPS",
    "SUBPOP_TO_SUPERPOP",
    "ALL_1KG_SUPERPOPS",
    "ALL_1KG_SUBPOPS",
    "POPULATION_NAMES",
    # Colors
    "SUPERPOP_COLORS",
    "POPULATION_COLORS",
    # Merge functions
    "validate_matrixtable_compatibility",
    "get_shared_variants_stats",
    "merge_matrixtables",
    "check_sample_overlap",
    # Filter functions
    "filter_to_autosomes",
    "filter_to_biallelic_snps",
    "compute_variant_qc",
    "filter_variants_for_ancestry",
    "ld_prune_variants",
    "prepare_ancestry_variants",
    # PCA functions
    "PCAResult",
    "compute_pca",
    "project_samples",
    # Classification functions
    "ClassificationResult",
    "train_classifier",
    "predict_ancestry",
    "validate_training_data",
    "get_training_samples",
    "get_query_samples",
    # Pipeline functions
    "PipelineConfig",
    "AncestryInferenceResult",
    "run_ancestry_inference",
]
