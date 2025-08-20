# Examples of Flexible Annotation Framework
# Demonstrates how to add custom annotations and future data sources

import logging
from hvantk.annotation.annotation_pipeline import (
    AnnotationConfig,
    FlexibleAnnotationStreamer,
    AnnotationRegistry,
    ConfigurableAnnotationPipeline,
    create_flexible_pipeline,
    add_custom_annotation
)
from hvantk.utils.clinvar_streamer import ClinvarDataStreamer
import hail as hl

logger = logging.getLogger(__name__)

def example_custom_bed_annotation():
    """Example: Adding custom BED file annotations"""
    print("=== Custom BED File Annotation ===")

    # Define custom annotation for regulatory regions
    regulatory_config = AnnotationConfig(
        name="regulatory_regions",
        source_path="./data/annotations/regulatory_regions.bed",
        annotation_type="region",
        feature_mapping={
            "target": "regulatory_target",
            "score": "regulatory_score"
        },
        metadata={
            "description": "ENCODE regulatory regions",
            "version": "v3.0",
            "source": "ENCODE"
        }
    )

    # Create streamer
    regulatory_streamer = FlexibleAnnotationStreamer(regulatory_config)

    print(f"Created streamer for {regulatory_config.name}")
    print(f"Annotation type: {regulatory_config.annotation_type}")
    print(f"Feature mapping: {regulatory_config.feature_mapping}")

    return regulatory_streamer


def example_custom_csv_annotation():
    """Example: Adding custom CSV file with gene annotations"""
    print("=== Custom CSV Gene Annotation ===")

    # Custom preprocessing function
    def preprocess_gene_scores(ht):
        """Custom preprocessing for gene score data"""
        return ht.annotate(
            normalized_score=ht.raw_score / 100.0,
            score_category=hl.case()
            .when(ht.raw_score > 80, "High")
            .when(ht.raw_score > 50, "Medium")
            .default("Low")
        )

    # Define custom gene annotation
    gene_scores_config = AnnotationConfig(
        name="custom_gene_scores",
        source_path="./data/annotations/gene_pathogenicity_scores.csv",
        annotation_type="gene",
        join_key="gene_symbol",
        preprocessing_func=preprocess_gene_scores,
        feature_mapping={
            "raw_score": "pathogenicity_raw",
            "normalized_score": "pathogenicity_norm"
        },
        metadata={
            "description": "Custom pathogenicity scores",
            "algorithm": "RandomForest_v2.1"
        }
    )

    # Create pipeline with custom annotation
    base_streamer = ClinvarDataStreamer("./data/clinvar/clinvar.vcf.gz")
    pipeline = create_flexible_pipeline(base_streamer, "CustomScorePipeline")
    pipeline.add_annotation(config=gene_scores_config)

    print(f"Added custom gene annotation: {gene_scores_config.name}")
    print(f"Preprocessing applied: {gene_scores_config.preprocessing_func is not None}")

    return pipeline


def example_registry_based_annotations():
    """Example: Using registry for organized annotation management"""
    print("=== Registry-Based Annotation Management ===")

    # Create custom registry
    registry = AnnotationRegistry()

    # Register multiple custom annotations
    annotations = [
        AnnotationConfig(
            name="splice_ai_scores",
            source_path="./data/annotations/spliceai_scores.tsv",
            annotation_type="variant",
            feature_mapping={"DS_AG": "splice_acceptor_gain"}
        ),
        AnnotationConfig(
            name="protein_domains",
            source_path="./data/annotations/pfam_domains.bed",
            annotation_type="region",
            metadata={"source": "Pfam", "version": "35.0"}
        ),
        AnnotationConfig(
            name="tissue_expression",
            source_path="./data/annotations/gtex_tissue_expr.tsv",
            annotation_type="gene",
            feature_mapping={"heart_tpm": "heart_expression"}
        )
    ]

    # Register by category
    for config in annotations:
        category = config.metadata.get("source", "custom")
        registry.register(config, category.lower())

    # Create pipeline and add by category
    base_streamer = ClinvarDataStreamer("./data/clinvar/clinvar.vcf.gz")
    pipeline = ConfigurableAnnotationPipeline("RegistryPipeline", base_streamer, registry)

    # Add all custom annotations
    for config in annotations:
        pipeline.add_annotation(config=config)

    print(f"Registered {len(annotations)} custom annotations")
    print(f"Categories: {registry._categories.keys()}")
    print(f"Pipeline has {len(pipeline.streamers)} streamers")

    return pipeline


def example_api_based_annotation():
    """Example: Adding annotation from API or database"""
    print("=== API-Based Annotation ===")

    def load_from_api(source_path):
        """Custom loader that fetches data from API"""
        # Simulate API call - in reality, you'd make HTTP requests
        # and convert response to Hail Table
        import pandas as pd

        # Mock API data
        mock_data = {
            'gene_symbol': ['GATA4', 'NKX2-5', 'TBX5'],
            'disease_association': [0.95, 0.87, 0.72],
            'confidence': ['High', 'High', 'Medium']
        }

        df = pd.DataFrame(mock_data)
        # Convert to Hail Table (simplified)
        return hl.Table.from_pandas(df, key='gene_symbol')

    # API-based annotation config
    api_config = AnnotationConfig(
        name="disease_api_scores",
        source_path="https://api.example.com/gene_disease_scores",
        annotation_type="gene",
        loader_func=load_from_api,
        metadata={
            "source": "DiseaseDB_API",
            "update_frequency": "weekly"
        }
    )

    # Create streamer
    api_streamer = FlexibleAnnotationStreamer(api_config)

    print(f"Created API-based streamer: {api_config.name}")
    print(f"Source: {api_config.source_path}")
    print(f"Custom loader: {api_config.loader_func is not None}")

    return api_streamer


def example_feature_transformations():
    """Example: Adding custom feature transformations"""
    print("=== Custom Feature Transformations ===")

    base_streamer = ClinvarDataStreamer("./data/clinvar/clinvar.vcf.gz")
    pipeline = create_flexible_pipeline(base_streamer, "TransformationPipeline")

    # Add some annotations
    pipeline.add_annotation(annotation_name="dbnsfp_scores")
    pipeline.add_annotation(annotation_name="gnomad_frequencies")

    # Add custom feature transformations
    def create_composite_score(ht):
        """Create composite pathogenicity score"""
        return ht.annotate(
            composite_pathogenicity=hl.case()
            .when(
                hl.is_defined(ht.cadd_score) & hl.is_defined(ht.allele_frequency),
                (ht.cadd_score / 30.0) * (-hl.log10(ht.allele_frequency + 1e-8) / 8.0)
            )
            .or_missing()
        )

    def add_risk_categories(ht):
        """Categorize variants by risk"""
        return ht.annotate(
            risk_category=hl.case()
            .when(ht.composite_pathogenicity > 0.8, "High_Risk")
            .when(ht.composite_pathogenicity > 0.5, "Medium_Risk")
            .when(hl.is_defined(ht.composite_pathogenicity), "Low_Risk")
            .default("Unknown_Risk")
        )

    # Add transformations to pipeline
    pipeline.add_feature_transformation(create_composite_score, "Composite pathogenicity score")
    pipeline.add_feature_transformation(add_risk_categories, "Risk categorization")

    print(f"Pipeline with {len(pipeline._feature_transformations)} custom transformations")

    return pipeline


def example_future_annotation_source():
    """Example: Framework for future/unknown annotation sources"""
    print("=== Future Annotation Source Example ===")

    # Hypothetical future annotation source
    def load_ai_predictions(source_path):
        """Loader for future AI-based variant predictions"""
        # This could be any future format: Parquet, Arrow, custom binary, etc.
        # The framework adapts without code changes

        # Simulate loading future AI prediction format
        mock_predictions = hl.utils.range_table(1000)
        mock_predictions = mock_predictions.annotate(
            locus=hl.locus("chr1", mock_predictions.idx + 100000),
            alleles=["A", "T"],
            ai_pathogenicity=hl.rand_unif(0, 1),
            ai_confidence=hl.rand_unif(0.5, 1.0),
            model_version="GPT-Genomics-v5.0"
        )
        return mock_predictions.key_by('locus', 'alleles')

    # Future AI annotation config
    ai_config = AnnotationConfig(
        name="future_ai_predictions",
        source_path="./data/future_annotations/ai_predictions_v5.parquet",
        annotation_type="variant",
        loader_func=load_ai_predictions,
        feature_mapping={
            "ai_pathogenicity": "ai_path_score",
            "ai_confidence": "ai_confidence_score"
        },
        metadata={
            "model": "GPT-Genomics",
            "version": "5.0",
            "training_date": "2025-01-01"
        }
    )

    # The framework handles it seamlessly
    ai_streamer = FlexibleAnnotationStreamer(ai_config)

    print(f"Future AI annotation ready: {ai_config.name}")
    print(f"Model: {ai_config.metadata.get('model')}")
    print(f"Framework adapts without code changes!")

    return ai_streamer


def example_complete_flexible_workflow():
    """Example: Complete workflow with multiple custom sources"""
    print("=== Complete Flexible Workflow ===")

    # Create base streamer
    base_streamer = ClinvarDataStreamer("./data/clinvar/clinvar.vcf.gz")

    # Create flexible pipeline
    pipeline = create_flexible_pipeline(base_streamer, "CompleteFlexiblePipeline")

    # Add built-in annotations
    pipeline.add_annotation(annotation_name="dbnsfp_scores")
    pipeline.add_annotation(annotation_name="gene_expression")

    # Add multiple custom annotations easily
    custom_sources = [
        ("regulatory_data", "./data/custom/regulatory.bed", "region"),
        ("protein_features", "./data/custom/protein_domains.tsv", "gene"),
        ("population_specific", "./data/custom/population_freq.vcf", "variant")
    ]

    for name, path, ann_type in custom_sources:
        add_custom_annotation(pipeline, name, path, ann_type)

    # Add feature engineering
    def final_score_computation(ht):
        """Compute final integrated score"""
        return ht.annotate(
            integrated_score=hl.case()
            .when(
                hl.is_defined(ht.cadd_score) &
                hl.is_defined(ht.median_expression) &
                hl.is_defined(ht.allele_frequency),
                (ht.cadd_score / 30.0 * 0.4 +
                 hl.log10(ht.median_expression + 0.1) / 3.0 * 0.3 +
                 -hl.log10(ht.allele_frequency + 1e-8) / 8.0 * 0.3)
            )
            .or_missing()
        )

    pipeline.add_feature_transformation(final_score_computation, "Integrated pathogenicity score")

    print(f"Complete pipeline ready with {len(pipeline.streamers)} annotation sources")
    print("Includes: Built-in annotations + Custom sources + Feature engineering")

    return pipeline


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    print("Flexible Annotation Framework Examples")
    print("=" * 50)

    print("\nExample functions available:")
    print("1. example_custom_bed_annotation() - BED file annotation")
    print("2. example_custom_csv_annotation() - CSV with preprocessing")
    print("3. example_registry_based_annotations() - Registry management")
    print("4. example_api_based_annotation() - API/database sources")
    print("5. example_feature_transformations() - Custom feature engineering")
    print("6. example_future_annotation_source() - Future/unknown formats")
    print("7. example_complete_flexible_workflow() - Complete workflow")

    print("\nFramework Benefits:")
    print("✅ No code changes needed for new annotation sources")
    print("✅ Automatic format detection and loading")
    print("✅ Flexible join strategies (variant/gene/region)")
    print("✅ Custom preprocessing and feature mapping")
    print("✅ Registry system for annotation discovery")
    print("✅ Feature transformation pipeline")
    print("✅ Error handling and fallback mechanisms")

    # Uncomment to run examples:
    # example_custom_bed_annotation()
    # example_registry_based_annotations()
    # example_complete_flexible_workflow()
