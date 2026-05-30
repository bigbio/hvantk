# ClinVar default wiring for the generic TrainingSetBuilder.
#
# Provides the ClinVar-specific pieces that the source-agnostic
# tools/training_sets/builder.py intentionally does not know about:
#   - default annotator stack,
#   - composite-feature post-processing (migrated from the old enhanced.py),
#   - a convenience entry point that assembles base + annotators + outcome.
#
# tools/ may compose skills/ (ClinVar streamer + label derivation) with
# algorithms/ (annotators), so all imports here are within the allowed layer.

from __future__ import annotations

from typing import Iterable, List, Optional

import hail as hl

from hvantk.algorithms.annotation.annotator import (
    Annotator,
    GeneConstraintAnnotator,
    GeneExpressionAnnotator,
    PopulationFrequencyAnnotator,
    VariantPredictionScoreAnnotator,
)
from hvantk.skills.clinvar.streamers import (
    ClinVarVariantTableStreamer,
    apply_clinvar_training_labels,
)
from hvantk.tools.training_sets.builder import OutcomeSpec, TrainingSetBuilder


def default_clinvar_annotators(
    *,
    tissue_focus: str = "heart",
    include_prediction_scores: bool = True,
    include_expression: bool = True,
    include_constraint: bool = True,
    include_population_freq: bool = True,
) -> List[Annotator]:
    """Build the default ClinVar annotator stack (mirrors the old enhanced.py)."""
    annotators: List[Annotator] = []
    if include_prediction_scores:
        annotators.append(VariantPredictionScoreAnnotator(""))
    if include_expression:
        annotators.append(GeneExpressionAnnotator("", tissue_focus))
    if include_constraint:
        annotators.append(GeneConstraintAnnotator())
    if include_population_freq:
        annotators.append(PopulationFrequencyAnnotator())
    return annotators


def clinvar_training_post_process(ht: hl.Table) -> hl.Table:
    """Add composite features combining multiple annotation sources.

    Migrated verbatim from the old EnhancedClinvarTrainingSetProcessor
    ``_add_final_features``: a composite ``pathogenicity_score`` and a
    ``feature_completeness`` coverage score. The column guards below mean the
    step degrades gracefully when an annotator was disabled and its columns are
    absent.
    """
    row_fields = set(ht.row.dtype.field_names())

    def _present(name: str):
        return hl.is_defined(ht[name]) if name in row_fields else hl.literal(False)

    # Composite pathogenicity score (requires deleteriousness + rarity +
    # constraint columns; missing where any are absent/undefined).
    if {"combined_deleteriousness", "rarity_score", "constraint_score"} <= row_fields:
        pathogenicity_score = (
            hl.case()
            .when(
                hl.is_defined(ht.combined_deleteriousness)
                & hl.is_defined(ht.rarity_score)
                & hl.is_defined(ht.constraint_score),
                (
                    ht.combined_deleteriousness * 0.4
                    + ht.rarity_score / 10.0 * 0.3
                    + ht.constraint_score * 0.3
                ),
            )
            .or_missing()
        )
    else:
        pathogenicity_score = hl.missing(hl.tfloat64)

    # Feature completeness score across the five core feature columns.
    feature_completeness = (
        hl.int32(_present("combined_deleteriousness"))
        + hl.int32(_present("conservation_score"))
        + hl.int32(_present("median_expression"))
        + hl.int32(_present("constraint_score"))
        + hl.int32(_present("AF"))
    ) / 5.0

    return ht.annotate(
        pathogenicity_score=pathogenicity_score,
        feature_completeness=feature_completeness,
    )


def build_clinvar_training_set(
    clinvar_table_path: str,
    *,
    gene_set: Optional[Iterable[str]] = None,
    disease_terms: Optional[Iterable[str]] = None,
    annotators: Optional[List[Annotator]] = None,
    output_path: Optional[str] = None,
    export_tsv: bool = False,
    tissue_focus: str = "heart",
    include_prediction_scores: bool = True,
    include_expression: bool = True,
    include_constraint: bool = True,
    include_population_freq: bool = True,
) -> hl.Table:
    """Build a feature-rich ClinVar training set as a Hail Table.

    Loads a pre-built ClinVar annotation-table artifact from
    ``clinvar_table_path``, derives TP/TN labels, applies the (default or
    supplied) annotator stack, and adds composite features.

    Note (legacy gene semantics): ``gene_set`` only GATES gene-based TP inside
    the label derivation; it does NOT hard-filter the table to those genes.
    """
    base = ClinVarVariantTableStreamer.from_path(clinvar_table_path)

    outcome = OutcomeSpec(
        column="rf_label",
        kind="label",
        derive=lambda ht: apply_clinvar_training_labels(
            ht,
            pathogenic_labels=base.PATHOGENIC_LABELS,
            benign_labels=base.BENIGN_LABELS,
            gene_set=gene_set,
            disease_terms=disease_terms,
        ),
    )

    anns = (
        annotators
        if annotators is not None
        else default_clinvar_annotators(
            tissue_focus=tissue_focus,
            include_prediction_scores=include_prediction_scores,
            include_expression=include_expression,
            include_constraint=include_constraint,
            include_population_freq=include_population_freq,
        )
    )

    return TrainingSetBuilder(
        base,
        anns,
        outcome,
        post_process=clinvar_training_post_process,
        name="ClinvarTrainingSetBuilder",
    ).build(output_path=output_path, export_tsv=export_tsv)
