# Flexible Annotation Framework
# Enhanced architecture for custom and future annotation sources

import hail as hl
from typing import Iterator, Optional, Set, List, Dict, Any, Union, Callable
from hvantk.data.data_streamer import HailDataStreamer, StreamProcessor
from abc import ABC, abstractmethod
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


class AnnotationConfig:
    """Configuration class for annotation sources"""

    def __init__(self,
                 name: str,
                 source_path: str,
                 annotation_type: str = "variant",  # "variant", "gene", "region"
                 join_key: Optional[str] = None,
                 loader_func: Optional[Callable] = None,
                 preprocessing_func: Optional[Callable] = None,
                 feature_mapping: Optional[Dict[str, str]] = None,
                 metadata: Optional[Dict[str, Any]] = None):
        self.name = name
        self.source_path = source_path
        self.annotation_type = annotation_type
        self.join_key = join_key or self._default_join_key()
        self.loader_func = loader_func
        self.preprocessing_func = preprocessing_func
        self.feature_mapping = feature_mapping or {}
        self.metadata = metadata or {}

    def _default_join_key(self) -> str:
        """Default join keys based on annotation type"""
        defaults = {
            "variant": "locus,alleles",
            "gene": "gene_symbol",
            "region": "locus",
            "transcript": "transcript_id"
        }
        return defaults.get(self.annotation_type, "key")


class FlexibleAnnotationStreamer(HailDataStreamer):
    """
    Highly flexible annotation streamer that can handle any data source
    through configuration rather than hard-coded implementations.
    """

    def __init__(self, config: AnnotationConfig, chunk_size: int = 10000):
        super().__init__(f"FlexibleAnnotation_{config.name}", chunk_size)
        self.config = config
        self.annotation_data = None

    def load_annotation_data(self) -> hl.Table:
        """Load annotation data using flexible configuration"""
        self.logger.info(f"Loading {self.config.name} from {self.config.source_path}")

        # Use custom loader if provided
        if self.config.loader_func:
            data = self.config.loader_func(self.config.source_path)
        else:
            # Auto-detect format and load
            data = self._auto_load_data(self.config.source_path)

        # Apply preprocessing if specified
        if self.config.preprocessing_func:
            data = self.config.preprocessing_func(data)

        return data

    def _auto_load_data(self, path: str) -> hl.Table:
        """Automatically detect file format and load appropriately"""
        path_obj = Path(path)
        suffix = path_obj.suffix.lower()

        if suffix == '.ht':
            return hl.read_table(path)
        elif suffix in ['.tsv', '.txt', '.csv']:
            delimiter = '\t' if suffix in ['.tsv', '.txt'] else ','
            return hl.import_table(path, delimiter=delimiter, impute=True)
        elif suffix in ['.vcf', '.vcf.gz']:
            return hl.import_vcf(path).rows()
        elif suffix in ['.bed', '.bed.gz']:
            return hl.import_bed(path)
        else:
            # Default to table import
            return hl.import_table(path, impute=True)

    def setup(self) -> None:
        """Load annotation data with error handling"""
        super().setup()
        try:
            self.annotation_data = self.load_annotation_data()
            self.logger.info(f"Loaded {self.config.name}: {self.annotation_data.count()} records")
        except Exception as e:
            self.logger.error(f"Failed to load {self.config.name}: {e}")
            raise

    def annotate_chunk(self, chunk: hl.Table) -> hl.Table:
        """Flexible annotation based on configuration"""
        if self.annotation_data is None:
            self.logger.warning(f"No annotation data loaded for {self.config.name}")
            return chunk

        try:
            # Determine join strategy based on annotation type
            if self.config.annotation_type == "variant":
                annotated = self._annotate_by_variant(chunk)
            elif self.config.annotation_type == "gene":
                annotated = self._annotate_by_gene(chunk)
            elif self.config.annotation_type == "region":
                annotated = self._annotate_by_region(chunk)
            else:
                annotated = self._annotate_custom(chunk)

            # Apply feature mapping if specified
            if self.config.feature_mapping:
                annotated = self._apply_feature_mapping(annotated)

            return annotated

        except Exception as e:
            self.logger.error(f"Annotation failed for {self.config.name}: {e}")
            return chunk  # Return original chunk on failure

    def _annotate_by_variant(self, chunk: hl.Table) -> hl.Table:
        """Annotate by variant (locus + alleles)"""
        return chunk.annotate(**self.annotation_data[chunk.key])

    def _annotate_by_gene(self, chunk: hl.Table) -> hl.Table:
        """Annotate by gene symbol"""
        if 'gene' not in chunk.row:
            self.logger.warning("No 'gene' field found in chunk for gene annotation")
            return chunk
        return chunk.annotate(**self.annotation_data[chunk.gene])

    def _annotate_by_region(self, chunk: hl.Table) -> hl.Table:
        """Annotate by genomic region overlap"""
        return chunk.annotate(**self.annotation_data[chunk.locus])

    def _annotate_custom(self, chunk: hl.Table) -> hl.Table:
        """Custom annotation logic - can be extended"""
        join_fields = self.config.join_key.split(',')
        if len(join_fields) == 1:
            key = chunk[join_fields[0]]
        else:
            key = hl.struct(**{field: chunk[field] for field in join_fields})
        return chunk.annotate(**self.annotation_data[key])

    def _apply_feature_mapping(self, annotated: hl.Table) -> hl.Table:
        """Apply feature name mapping for standardization"""
        renames = {}
        for old_name, new_name in self.config.feature_mapping.items():
            if old_name in annotated.row:
                renames[new_name] = annotated[old_name]
        return annotated.annotate(**renames) if renames else annotated

    def stream(self) -> Iterator[hl.Table]:
        """Not implemented for annotation streamers"""
        raise NotImplementedError("FlexibleAnnotationStreamer is for processing existing chunks")

    def process_chunk(self, chunk: hl.Table) -> hl.Table:
        """Process chunk by adding annotations"""
        return self.annotate_chunk(chunk)


class AnnotationRegistry:
    """Registry for managing and discovering annotation sources"""

    def __init__(self):
        self._annotations = {}
        self._categories = {}

    def register(self, config: AnnotationConfig, category: str = "general"):
        """Register an annotation source"""
        self._annotations[config.name] = config
        if category not in self._categories:
            self._categories[category] = []
        self._categories[category].append(config.name)
        logger.info(f"Registered annotation: {config.name} (category: {category})")

    def get(self, name: str) -> Optional[AnnotationConfig]:
        """Get annotation config by name"""
        return self._annotations.get(name)

    def list_by_category(self, category: str) -> List[str]:
        """List annotations in a category"""
        return self._categories.get(category, [])

    def list_all(self) -> List[str]:
        """List all registered annotations"""
        return list(self._annotations.keys())

    def create_streamer(self, name: str, **kwargs) -> FlexibleAnnotationStreamer:
        """Create a streamer for a registered annotation"""
        config = self.get(name)
        if not config:
            raise ValueError(f"Annotation '{name}' not registered")
        return FlexibleAnnotationStreamer(config, **kwargs)


class ConfigurableAnnotationPipeline(StreamProcessor):
    """
    Highly configurable annotation pipeline that can be built from
    registry entries or custom configurations.
    """

    def __init__(self,
                 name: str,
                 base_streamer: HailDataStreamer,
                 registry: Optional[AnnotationRegistry] = None):
        super().__init__(name)
        self.registry = registry or AnnotationRegistry()
        self.add_streamer(base_streamer)
        self._feature_transformations = []

    def add_annotation(self,
                      annotation_name: str = None,
                      config: AnnotationConfig = None,
                      **streamer_kwargs) -> 'ConfigurableAnnotationPipeline':
        """Add annotation by name (from registry) or direct config"""
        if config:
            streamer = FlexibleAnnotationStreamer(config, **streamer_kwargs)
        elif annotation_name and self.registry:
            streamer = self.registry.create_streamer(annotation_name, **streamer_kwargs)
        else:
            raise ValueError("Either annotation_name or config must be provided")

        self.add_streamer(streamer)
        return self

    def add_annotations_by_category(self,
                                   category: str,
                                   **streamer_kwargs) -> 'ConfigurableAnnotationPipeline':
        """Add all annotations from a category"""
        annotation_names = self.registry.list_by_category(category)
        for name in annotation_names:
            self.add_annotation(annotation_name=name, **streamer_kwargs)
        return self

    def add_feature_transformation(self,
                                  transform_func: Callable[[hl.Table], hl.Table],
                                  description: str = "Custom transformation"):
        """Add custom feature transformation to be applied after annotations"""
        self._feature_transformations.append((transform_func, description))
        return self

    def process(self, output_path: Optional[str] = None) -> Optional[hl.Table]:
        """Process with feature transformations"""
        # Run standard pipeline
        result = super().process(output_path=None)  # Don't save yet

        if result is None:
            return None

        # Apply feature transformations
        for transform_func, description in self._feature_transformations:
            self.logger.info(f"Applying transformation: {description}")
            if isinstance(result, list):
                result = [transform_func(chunk) for chunk in result]
            else:
                result = transform_func(result)

        # Save if output path specified
        if output_path and result:
            if isinstance(result, list):
                # Union chunks if needed
                final_result = result[0]
                for chunk in result[1:]:
                    final_result = final_result.union(chunk)
                final_result = final_result.checkpoint(output_path, overwrite=True)
            else:
                result = result.checkpoint(output_path, overwrite=True)

        return result


# Built-in annotation configurations for common sources

def create_builtin_registry() -> AnnotationRegistry:
    """Create registry with built-in annotation configurations"""
    registry = AnnotationRegistry()

    # Variant prediction scores
    registry.register(
        AnnotationConfig(
            name="dbnsfp_scores",
            source_path="",  # Will be set by dataset function
            annotation_type="variant",
            loader_func=lambda _: __import__('hvantk.data.dataset', fromlist=['get_dbnsfp_scores_ht']).get_dbnsfp_scores_ht(),
            feature_mapping={
                "CADD_phred": "cadd_score",
                "REVEL_score": "revel_score"
            },
            metadata={"category": "prediction", "data_type": "scores"}
        ),
        category="prediction"
    )

    # Gene expression
    registry.register(
        AnnotationConfig(
            name="gene_expression",
            source_path="",
            annotation_type="gene",
            loader_func=lambda _: __import__('hvantk.data.dataset', fromlist=['get_gene_expression_ht']).get_gene_expression_ht(),
            metadata={"category": "expression", "data_type": "levels"}
        ),
        category="expression"
    )

    # Population frequencies
    registry.register(
        AnnotationConfig(
            name="gnomad_frequencies",
            source_path="",
            annotation_type="variant",
            loader_func=lambda _: __import__('hvantk.data.dataset', fromlist=['get_gnomad_af_ht']).get_gnomad_af_ht(),
            feature_mapping={
                "AF": "allele_frequency",
                "AC": "allele_count"
            },
            metadata={"category": "population", "data_type": "frequencies"}
        ),
        category="population"
    )

    return registry


# Factory functions for easy usage

def create_flexible_pipeline(base_streamer: HailDataStreamer,
                           pipeline_name: str = "FlexibleAnnotationPipeline") -> ConfigurableAnnotationPipeline:
    """Create a flexible annotation pipeline with built-in registry"""
    registry = create_builtin_registry()
    return ConfigurableAnnotationPipeline(pipeline_name, base_streamer, registry)


def add_custom_annotation(pipeline: ConfigurableAnnotationPipeline,
                         name: str,
                         source_path: str,
                         annotation_type: str = "variant",
                         **config_kwargs) -> ConfigurableAnnotationPipeline:
    """Helper to add custom annotation to existing pipeline"""
    config = AnnotationConfig(
        name=name,
        source_path=source_path,
        annotation_type=annotation_type,
        **config_kwargs
    )
    return pipeline.add_annotation(config=config)

