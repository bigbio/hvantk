# Enhanced Training Set Examples
# Demonstrates different annotation combinations and use cases

import logging
from hvantk.annotation.annotation_streamer import (
    create_enhanced_clinvar_training_streamer,
    VariantPredictionScoreStreamer,
    GeneConstraintStreamer,
)
from hvantk.data.data_streamer import StreamProcessor
from hvantk.utils.clinvar_streamer import ClinvarDataStreamer

logger = logging.getLogger(__name__)

def example_full_annotation_pipeline():
    """Example: Full annotation pipeline with all features"""
    print("=== Full Annotation Pipeline ===")

    processor = create_enhanced_clinvar_training_streamer(
        clinvar_path="./data/clinvar/clinvar_sample.vcf.gz",
        output_dir="./data/training_set/full",
        tissue_focus="heart",
        include_prediction_scores=True,  # CADD, SIFT, PolyPhen, REVEL
        include_expression=True,         # Tissue-specific expression
        include_constraint=True,         # pLI, LOEUF scores
        include_population_freq=True     # gnomAD frequencies
    )

    result = processor.process()

    if result:
        print(f"Generated {result.count()} fully annotated training examples")
        print("Features include:")
        print("  - Clinvar P/B labels")
        print("  - Variant prediction scores (CADD, SIFT, etc.)")
        print("  - Gene expression levels")
        print("  - Constraint metrics (pLI, LOEUF)")
        print("  - Population frequencies")
        print("  - Derived composite scores")

    return result


def example_custom_annotation_pipeline():
    """Example: Custom pipeline with specific annotators"""
    print("=== Custom Annotation Pipeline ===")

    # Build custom pipeline step by step
    pipeline = StreamProcessor("CustomAnnotationPipeline")

    # Add base Clinvar data
    chd_genes = {"GATA4", "NKX2-5", "TBX5", "NOTCH1"}
    clinvar_streamer = ClinvarDataStreamer(
        clinvar_path="./data/clinvar/clinvar_sample.vcf.gz",
        chd_genes=chd_genes
    )
    pipeline.add_streamer(clinvar_streamer)

    # Add only specific annotations
    pipeline.add_streamer(VariantPredictionScoreStreamer(""))  # Prediction scores only
    pipeline.add_streamer(GeneConstraintStreamer())           # Constraint metrics only

    result = pipeline.process("./data/training_set/custom.ht")

    if result:
        print(f"Generated {len(result)} custom annotated chunks")
        print("Includes only prediction scores and constraint metrics")

    return result


def example_tissue_specific_analysis():
    """Example: Tissue-specific expression analysis"""
    print("=== Tissue-Specific Analysis ===")

    # Heart-focused analysis
    heart_processor = create_enhanced_clinvar_training_streamer(
        clinvar_path="./data/clinvar/clinvar_sample.vcf.gz",
        output_dir="./data/training_set/heart",
        tissue_focus="heart",
        include_prediction_scores=False,
        include_expression=True,      # Focus on expression
        include_constraint=True,
        include_population_freq=False
    )

    heart_result = heart_processor.process()

    if heart_result:
        print(f"Heart-focused training set: {heart_result.count()} variants")

        # Show heart-enriched genes
        heart_enriched = heart_result.filter(heart_result.heart_enriched).count()
        print(f"Heart-enriched variants: {heart_enriched}")

    return heart_result


def example_prediction_score_analysis():
    """Example: Focus on variant prediction scores"""
    print("=== Prediction Score Analysis ===")

    processor = create_enhanced_clinvar_training_streamer(
        clinvar_path="./data/clinvar/clinvar_sample.vcf.gz",
        output_dir="./data/training_set/prediction",
        include_prediction_scores=True,   # Main focus
        include_expression=False,
        include_constraint=False,
        include_population_freq=True      # Add frequency for context
    )

    result = processor.process()

    if result:
        print(f"Prediction-focused training set: {result.count()} variants")

        # Analyze score distributions
        if hasattr(result, 'combined_deleteriousness'):
            high_impact = result.filter(result.combined_deleteriousness > 0.8).count()
            print(f"High-impact variants (score > 0.8): {high_impact}")

    return result


def example_feature_engineering_showcase():
    """Example: Showcase of derived features"""
    print("=== Feature Engineering Showcase ===")

    processor = create_enhanced_clinvar_training_streamer(
        clinvar_path="./data/clinvar/clinvar_sample.vcf.gz",
        output_dir="./data/training_set/engineered"
    )

    result = processor.process()

    if result:
        print(f"Feature-engineered training set: {result.count()} variants")
        print("\nDerived features include:")
        print("  - combined_deleteriousness: Weighted combination of prediction scores")
        print("  - conservation_score: Evolutionary conservation metrics")
        print("  - tissue_specificity: Tissue-specific expression ratios")
        print("  - haploinsufficiency_score: Gene haploinsufficiency classification")
        print("  - rarity_score: -log10 transformed allele frequency")
        print("  - pathogenicity_score: Composite pathogenicity prediction")
        print("  - feature_completeness: Fraction of features available per variant")

        # Show feature completeness distribution
        if hasattr(result, 'feature_completeness'):
            complete_variants = result.filter(result.feature_completeness > 0.8).count()
            print(f"\nVariants with >80% feature completeness: {complete_variants}")

    return result


def example_machine_learning_ready():
    """Example: Generate ML-ready dataset with balanced features"""
    print("=== Machine Learning Ready Dataset ===")

    processor = create_enhanced_clinvar_training_streamer(
        clinvar_path="./data/clinvar/clinvar_sample.vcf.gz",
        output_dir="./data/training_set/ml_ready"
    )

    result = processor.process()

    if result:
        # Additional ML preprocessing
        ml_ready = result.select(
            'rf_label',  # Target variable
            # Core features
            'combined_deleteriousness',
            'conservation_score',
            'rarity_score',
            'constraint_score',
            # Categorical features
            'functional_consensus',
            'expression_level',
            'haploinsufficiency_score',
            'frequency_category',
            # Derived features
            'pathogenicity_score',
            'feature_completeness'
        )

        # Filter to high-quality examples
        ml_ready = ml_ready.filter(ml_ready.feature_completeness >= 0.6)

        print(f"ML-ready dataset: {ml_ready.count()} high-quality variants")
        print("Features selected for model training:")
        print("  - Prediction scores and conservation")
        print("  - Population frequency and constraint")
        print("  - Functional and expression categories")
        print("  - Composite pathogenicity score")

        # Export as ML-friendly format
        ml_ready.export("./data/training_set/ml_ready.tsv")
        print("Exported ML-ready TSV file")

    return result


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    print("Enhanced Training Set Examples")
    print("=" * 50)

    # Note: These examples assume you have annotation data available
    print("\nAvailable examples:")
    print("1. example_full_annotation_pipeline() - All features")
    print("2. example_custom_annotation_pipeline() - Custom selection")
    print("3. example_tissue_specific_analysis() - Heart-focused")
    print("4. example_prediction_score_analysis() - Prediction scores only")
    print("5. example_feature_engineering_showcase() - Derived features")
    print("6. example_machine_learning_ready() - ML preprocessing")

    print("\nTo run examples:")
    print("1. Ensure annotation data sources are available")
    print("2. Update file paths in examples")
    print("3. Uncomment desired example calls")

    # Uncomment to run specific examples:
    # example_full_annotation_pipeline()
    # example_feature_engineering_showcase()
    # example_machine_learning_ready()
