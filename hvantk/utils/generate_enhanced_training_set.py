# Enhanced Training Set Generation with Multi-Source Annotations
# Generates feature-rich training sets by combining Clinvar labels with variant and gene annotations

import logging
import os
import hail as hl
from hvantk.hgc.constants import VCF_EXTENSION
from hvantk.annotation.annotation_streamer import create_enhanced_clinvar_training_streamer
from hvantk.utils.gene_sets import load_sample_chd_gene_set

logger = logging.getLogger(__name__)

def main():
    """Main function to generate enhanced training set with multiple annotation sources"""

    # Configuration
    output_dir = "./data/training_set"
    clinvar_path = os.environ.get("CLINVAR_VCF", f"./data/clinvar/clinvar_20220403{VCF_EXTENSION}")

    logger.info(f"Starting enhanced Clinvar training set generation")
    logger.info(f"Clinvar path: {clinvar_path}")
    logger.info(f"Output directory: {output_dir}")

    # Load CHD gene set
    logger.info("Loading CHD-associated genes")
    gene_set = load_sample_chd_gene_set()
    logger.info(f"Loaded {len(gene_set)} sample genes")

    # Create the enhanced streaming processor with all annotation sources
    processor = create_enhanced_clinvar_training_streamer(
        clinvar_path=clinvar_path,
        output_dir=output_dir,
        gene_set=gene_set,
        tissue_focus="heart",  # Focus on heart tissue for CHD
        include_prediction_scores=True,  # CADD, SIFT, PolyPhen, REVEL, etc.
        include_expression=True,         # Gene expression across tissues
        include_constraint=True,         # pLI, LOEUF, constraint metrics
        include_population_freq=True     # gnomAD allele frequencies
    )

    # Process the data through annotation pipeline
    try:
        training_set = processor.process()
        if training_set:
            logger.info(f"Successfully generated enhanced training set with {training_set.count()} variants")

            # Show detailed statistics
            tp_count = training_set.filter(training_set.rf_label == "TP").count()
            tn_count = training_set.filter(training_set.rf_label == "TN").count()

            logger.info(f"Training set statistics:")
            logger.info(f"  True Positives (TP): {tp_count}")
            logger.info(f"  True Negatives (TN): {tn_count}")
            logger.info(f"  Total: {tp_count + tn_count}")

            # Feature coverage statistics (single-pass; guard absent fields)
            row_fields = set(training_set.row.dtype.field_names())
            feature_stats = training_set.aggregate(
                hl.struct(
                    has_prediction_scores=(
                        hl.agg.fraction(hl.is_defined(training_set.combined_deleteriousness))
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

            logger.info(f"Feature coverage statistics:")
            logger.info(f"  Prediction scores: {feature_stats['has_prediction_scores']:.2%}")
            logger.info(f"  Gene expression: {feature_stats['has_expression']:.2%}")
            logger.info(f"  Constraint metrics: {feature_stats['has_constraint']:.2%}")
            logger.info(f"  Population frequency: {feature_stats['has_frequency']:.2%}")
            logger.info(f"  Average feature completeness: {feature_stats['avg_feature_completeness']:.2f}")

            # Show top features by pathogenicity score
            if 'pathogenicity_score' in training_set.row.dtype.fields:
                top_pathogenic = training_set.filter(
                    hl.is_defined(training_set.pathogenicity_score)
                ).order_by(hl.desc(training_set.pathogenicity_score)).take(5)

                logger.info("Top 5 variants by pathogenicity score:")
                for variant in top_pathogenic:
                    logger.info(f"  {variant.gene}: {variant.pathogenicity_score:.3f} (Label: {variant.rf_label})")

        else:
            logger.warning("No enhanced training set generated")

    except Exception as e:
        logger.error(f"Error generating enhanced training set: {e}")
        raise


def create_minimal_feature_set():
    """Create a training set with minimal annotations for quick testing"""

    logger.info("Creating minimal feature training set")

    processor = create_enhanced_clinvar_training_streamer(
        clinvar_path=os.environ.get("CLINVAR_VCF"),
        output_dir="./data/training_set",
        include_prediction_scores=True,   # Only include prediction scores
        include_expression=False,         # Skip expression data
        include_constraint=True,          # Include constraint metrics
        include_population_freq=True      # Include frequency data
    )

    return processor.process()


def create_expression_focused_set():
    """Create a training set focused on expression and tissue-specific features"""

    logger.info("Creating expression-focused training set")

    processor = create_enhanced_clinvar_training_streamer(
        clinvar_path=os.environ.get("CLINVAR_VCF"),
        output_dir="./data/training_set",
        tissue_focus="heart",
        include_prediction_scores=False,  # Skip prediction scores
        include_expression=True,          # Focus on expression
        include_constraint=True,          # Include constraint
        include_population_freq=False     # Skip frequency data
    )

    return processor.process()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    main()
