# Enhanced Training Set Generation with Multi-Source Annotations
#
# Demo driver: builds a feature-rich ClinVar training set by composing the
# ClinVar skill data source with algorithm-level annotators, via the generic
# TrainingSetBuilder (tools/training_sets/builder.py) and its ClinVar wiring
# (tools/training_sets/clinvar.py).
#
# Input is a pre-built ClinVar AnnotationTable artifact (a `.ht` directory
# loaded via VariantTableStreamer.from_path), NOT a raw VCF. Build it first
# with `hvantk reprocess clinvar:variants ...`, then point CLINVAR_HT at it.

import logging
import os

import hail as hl

from hvantk.core.utils.gene_sets import load_sample_chd_gene_set
from hvantk.tools.training_sets.clinvar import build_clinvar_training_set

logger = logging.getLogger(__name__)

DEFAULT_CLINVAR_HT = "./data/clinvar/clinvar.ht"


def main():
    """Generate an enhanced training set with multiple annotation sources."""

    # Configuration
    output_dir = "./data/training_set"
    output_path = f"{output_dir}/enhanced_ts.clinvar.ht"
    clinvar_table_path = os.environ.get("CLINVAR_HT", DEFAULT_CLINVAR_HT)

    logger.info("Starting enhanced Clinvar training set generation")
    logger.info(f"Clinvar table path: {clinvar_table_path}")
    logger.info(f"Output directory: {output_dir}")

    # Load CHD gene set (gates gene-based TP labeling; does not hard-filter).
    logger.info("Loading CHD-associated genes")
    gene_set = load_sample_chd_gene_set()
    logger.info(f"Loaded {len(gene_set)} sample genes")

    try:
        training_set = build_clinvar_training_set(
            clinvar_table_path=clinvar_table_path,
            gene_set=gene_set,
            output_path=output_path,
            export_tsv=True,
            tissue_focus="heart",  # Focus on heart tissue for CHD
            include_prediction_scores=True,  # CADD, SIFT, PolyPhen, REVEL, etc.
            include_expression=True,  # Gene expression across tissues
            include_constraint=True,  # pLI, LOEUF, constraint metrics
            include_population_freq=True,  # gnomAD allele frequencies
        )

        if training_set is not None:
            logger.info(
                f"Successfully generated enhanced training set with "
                f"{training_set.count()} variants"
            )
            _log_statistics(training_set)
        else:
            logger.warning("No enhanced training set generated")

    except Exception as e:
        logger.error(f"Error generating enhanced training set: {e}")
        raise


def _log_statistics(training_set: hl.Table) -> None:
    """Log TP/TN counts and feature-coverage statistics."""

    # Show detailed statistics
    tp_count = training_set.filter(training_set.rf_label == "TP").count()
    tn_count = training_set.filter(training_set.rf_label == "TN").count()

    logger.info("Training set statistics:")
    logger.info(f"  True Positives (TP): {tp_count}")
    logger.info(f"  True Negatives (TN): {tn_count}")
    logger.info(f"  Total: {tp_count + tn_count}")

    # Feature coverage statistics (single-pass; guard absent fields)
    row_fields = set(training_set.row.dtype.field_names())
    feature_stats = training_set.aggregate(
        hl.struct(
            has_prediction_scores=(
                hl.agg.fraction(
                    hl.is_defined(training_set.combined_deleteriousness)
                )
                if "combined_deleteriousness" in row_fields
                else hl.agg.fraction(hl.literal(False))
            ),
            has_expression=(
                hl.agg.fraction(hl.is_defined(training_set.median_expression))
                if "median_expression" in row_fields
                else hl.agg.fraction(hl.literal(False))
            ),
            has_constraint=(
                hl.agg.fraction(hl.is_defined(training_set.constraint_score))
                if "constraint_score" in row_fields
                else hl.agg.fraction(hl.literal(False))
            ),
            has_frequency=(
                hl.agg.fraction(hl.is_defined(training_set.AF))
                if "AF" in row_fields
                else hl.agg.fraction(hl.literal(False))
            ),
            # Mean over all rows (returns missing if field absent or always missing)
            avg_feature_completeness=(
                hl.agg.mean(training_set.feature_completeness)
                if "feature_completeness" in row_fields
                else hl.agg.mean(hl.null(hl.tfloat64))
            ),
        )
    )

    logger.info("Feature coverage statistics:")
    logger.info(
        f"  Prediction scores: {feature_stats['has_prediction_scores']:.2%}"
    )
    logger.info(f"  Gene expression: {feature_stats['has_expression']:.2%}")
    logger.info(f"  Constraint metrics: {feature_stats['has_constraint']:.2%}")
    logger.info(f"  Population frequency: {feature_stats['has_frequency']:.2%}")
    logger.info(
        f"  Average feature completeness: "
        f"{feature_stats['avg_feature_completeness']:.2f}"
    )

    # Show top features by pathogenicity score
    if "pathogenicity_score" in training_set.row.dtype.fields:
        top_pathogenic = (
            training_set.filter(hl.is_defined(training_set.pathogenicity_score))
            .order_by(hl.desc(training_set.pathogenicity_score))
            .take(5)
        )

        logger.info("Top 5 variants by pathogenicity score:")
        for variant in top_pathogenic:
            logger.info(
                f"  {variant.gene}: {variant.pathogenicity_score:.3f} "
                f"(Label: {variant.rf_label})"
            )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    main()
