# Flexible Annotation Framework
# Enhanced architecture for custom and future annotation sources

import hail as hl
from typing import Iterator, Optional, List, Dict, Any, Callable
from hvantk.data.data_streamer import HailDataStreamer, StreamProcessor
import logging
import os
import warnings
from pathlib import Path

logger = logging.getLogger(__name__)


class AnnotationConfig:
    """Configuration class for annotation sources"""

    def __init__(
        self,
        name: str,
        source_path: str,
        annotation_type: str = "variant",  # "variant", "gene", "region"
        join_key: Optional[str] = None,
        loader_func: Optional[Callable] = None,
        preprocessing_func: Optional[Callable] = None,
        feature_mapping: Optional[Dict[str, str]] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ):
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
            "transcript": "transcript_id",
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
        suffixes = [s.lower() for s in path_obj.suffixes]
        path_lower = str(path).lower()

        # Helper: warn if large and not .bgz
        def warn_if_large_and_not_bgz(path, path_lower):
            try:
                file_size = os.path.getsize(path)
                if file_size > 100 * 1024 * 1024 and not path_lower.endswith(".bgz"):
                    warnings.warn(
                        f"File {path} is large and not block-compressed (.bgz). Hail will process it with a single CPU, which may be slow. For parallel processing, use bgzip (.bgz) compression.",
                        UserWarning,
                        stacklevel=2,
                    )
            except Exception:
                pass

        # VCF detection (compressed and uncompressed)
        if (
            path_lower.endswith(".vcf")
            or path_lower.endswith(".vcf.gz")
            or path_lower.endswith(".vcf.bgz")
            or (".vcf" in suffixes and (".gz" in suffixes or ".bgz" in suffixes))
        ):
            warn_if_large_and_not_bgz(path, path_lower)
            return hl.import_vcf(path).rows()
        # BED detection (compressed and uncompressed)
        elif (
            path_lower.endswith(".bed")
            or path_lower.endswith(".bed.gz")
            or path_lower.endswith(".bed.bgz")
            or (".bed" in suffixes and (".gz" in suffixes or ".bgz" in suffixes))
        ):
            warn_if_large_and_not_bgz(path, path_lower)
            return hl.import_bed(path)
        # Hail Table
        elif suffix == ".ht":
            return hl.read_table(path)
        # Tabular files (.tsv/.csv/.txt, compressed or not)
        elif any(
            path_lower.endswith(ext)
            or path_lower.endswith(ext + ".gz")
            or path_lower.endswith(ext + ".bgz")
            for ext in [".tsv", ".csv", ".txt"]
        ) or any(
            ext in suffixes and (".gz" in suffixes or ".bgz" in suffixes)
            for ext in [".tsv", ".csv", ".txt"]
        ):
            warn_if_large_and_not_bgz(path, path_lower)
            delimiter = "\t" if any(e in path_lower for e in [".tsv", ".txt"]) else ","
            return hl.import_table(path, delimiter=delimiter, impute=True)
        else:
            # Default to table import
            warn_if_large_and_not_bgz(path, path_lower)
            return hl.import_table(path, impute=True)

    def _ensure_annotation_key(self) -> None:
        """
        Ensure self.annotation_data is keyed on the expected join fields for
        efficient/valid self.annotation_data[key_expr] lookups.
        """
        ht = self.annotation_data
        if not isinstance(ht, hl.Table):
            return

        try:
            if self.config.annotation_type == "variant":
                # Prefer (locus, alleles) if present
                row_fields = set(ht.row.dtype.fields)
                if {"locus", "alleles"} <= row_fields:
                    self.annotation_data = ht.key_by("locus", "alleles")
                    return

            elif self.config.annotation_type == "gene":
                # Try common gene identifiers in order of preference (case-insensitive)
                row_fields = set(ht.row.dtype.fields)
                candidates = ("gene_symbol", "gene", "Gene", "gene_name", "gene_id")
                # Build a mapping of lowercase field names to actual field names
                row_fields_lc = {f.lower(): f for f in row_fields}
                for candidate in candidates:
                    candidate_lc = candidate.lower()
                    if candidate_lc in row_fields_lc:
                        self.annotation_data = ht.key_by(row_fields_lc[candidate_lc])
                        return

            elif self.config.annotation_type == "region":
                # Interval-keyed tables support point lookups with a locus
                row_fields = set(ht.row.dtype.fields)
                if "interval" in row_fields:
                    self.annotation_data = ht.key_by("interval")
                    return
                if "locus" in row_fields:
                    self.annotation_data = ht.key_by("locus")
                    return

            # Fallback: use explicit join_key if provided and fields exist
            join_fields = [
                f.strip() for f in (self.config.join_key or "").split(",") if f.strip()
            ]
            if join_fields:
                row_fields = set(ht.row.dtype.fields)
                if all(f in row_fields for f in join_fields):
                    self.annotation_data = ht.key_by(*join_fields)
                    return

            # If we reach here, we did not rekey. Proceed as-is; _safe_apply will guard.
            self.logger.debug(
                f"No suitable key inferred for {self.config.name}; proceeding without rekey"
            )

        except Exception as ke:
            self.logger.warning(
                f"Failed to set key on annotation data for {self.config.name}: {ke}"
            )

    def setup(self) -> None:
        """Load annotation data with error handling"""
        super().setup()
        try:
            self.annotation_data = self.load_annotation_data()
            # Ensure proper keying for subsequent index lookups; avoid full counts here.
            self._ensure_annotation_key()
            key_descr = [str(k) for k in self.annotation_data.key]  # informational only
            row_fields = list(self.annotation_data.row.dtype.fields)
            self.logger.info(
                f"Loaded {self.config.name} (fields={len(row_fields)}, key={key_descr})"
            )
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

    # -------------------- Internal helpers for safer annotation --------------------
    def _check_fields(self, chunk: hl.Table, required: List[str]) -> bool:
        row_fields = set(chunk.row.dtype.fields)
        missing = [f for f in required if f not in row_fields]
        if missing:
            self.logger.warning(
                f"Missing required field(s) {missing} for annotation type '{self.config.annotation_type}' in {self.config.name}; skipping annotation on this chunk"
            )
            return False
        return True

    def _safe_apply(self, chunk: hl.Table, key_expr) -> hl.Table:
        """Attempt to index annotation_data; if key missing or fails, return chunk unmodified.
        This protects against runtime failures due to absent keys."""
        try:
            ann_row = self.annotation_data[key_expr]
            # Use the annotation row's field names for annotation, skipping join-key fields and existing chunk fields
            field_names = list(self.annotation_data.row.dtype.fields)
            # Determine join-key fields to skip
            if self.config.annotation_type == "variant":
                skip_fields = {"locus", "alleles"}
            elif self.config.annotation_type == "gene":
                skip_fields = {"gene", "gene_symbol"}
            elif self.config.annotation_type == "region":
                skip_fields = {"interval", "locus"}
            else:
                skip_fields = set(f.strip() for f in (self.config.join_key or "").split(",") if f.strip())
            # Also skip any fields already present in chunk
            chunk_fields = set(chunk.row.dtype.fields)
            fields_to_copy = [fname for fname in field_names if fname not in skip_fields and fname not in chunk_fields]
            if not fields_to_copy:
                return chunk
            return chunk.annotate(**{fname: ann_row[fname] for fname in fields_to_copy})
        except Exception as e:
            self.logger.warning(
                f"Safe annotation lookup failed for {self.config.name}: {e}; returning original chunk"
            )
            return chunk

    def _annotate_by_variant(self, chunk: hl.Table) -> hl.Table:
        """Annotate by variant (locus + alleles). Validates key fields exist."""
        if not self._check_fields(chunk, ["locus", "alleles"]):
            return chunk
        key_expr = hl.struct(locus=chunk.locus, alleles=chunk.alleles)
        return self._safe_apply(chunk, key_expr)

    def _annotate_by_gene(self, chunk: hl.Table) -> hl.Table:
        """Annotate by gene symbol. Accepts common gene field names, case-insensitive."""
        accepted_names = {"gene", "gene_symbol", "gene_id", "gene_name"}
        row_fields = list(chunk.row.dtype.fields)
        field_name = None
        for f in row_fields:
            if f.lower() in accepted_names:
                field_name = f
                break
        if not field_name:
            self.logger.warning(
                f"No accepted gene field found in chunk for gene annotation; present fields: {row_fields}"
            )
            return chunk
        return self._safe_apply(chunk, chunk[field_name])

    def _annotate_by_region(self, chunk: hl.Table) -> hl.Table:
        """Annotate by genomic region overlap (expects 'locus')."""
        if not self._check_fields(chunk, ["locus"]):
            return chunk
        return self._safe_apply(chunk, chunk.locus)

    def _annotate_custom(self, chunk: hl.Table) -> hl.Table:
        """Custom annotation logic - validates all join fields."""
        join_fields = [f.strip() for f in self.config.join_key.split(",") if f.strip()]
        if not join_fields:
            self.logger.warning(
                f"No join_key specified for custom annotation {self.config.name}; skipping"
            )
            return chunk
        if not self._check_fields(chunk, join_fields):
            return chunk
        if len(join_fields) == 1:
            key_expr = chunk[join_fields[0]]
        else:
            key_expr = hl.struct(**{field: chunk[field] for field in join_fields})
        return self._safe_apply(chunk, key_expr)

    def _apply_feature_mapping(self, annotated: hl.Table) -> hl.Table:
        """Apply feature name mapping for standardization"""
        renames = {}
        for old_name, new_name in self.config.feature_mapping.items():
            if old_name in annotated.row:
                renames[new_name] = annotated[old_name]
        return annotated.annotate(**renames) if renames else annotated

    def stream(self) -> Iterator[hl.Table]:
        """Not implemented for annotation streamers"""
        raise NotImplementedError(
            "FlexibleAnnotationStreamer is for processing existing chunks"
        )

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

    def __init__(
        self,
        name: str,
        base_streamer: HailDataStreamer,
        registry: Optional[AnnotationRegistry] = None,
    ):
        super().__init__(name)
        self.registry = registry or AnnotationRegistry()
        self.add_streamer(base_streamer)
        self._feature_transformations = []

    def add_annotation(
        self,
        annotation_name: str = None,
        config: AnnotationConfig = None,
        **streamer_kwargs,
    ) -> "ConfigurableAnnotationPipeline":
        """Add annotation by name (from registry) or direct config"""
        if config:
            streamer = FlexibleAnnotationStreamer(config, **streamer_kwargs)
        elif annotation_name and self.registry:
            streamer = self.registry.create_streamer(annotation_name, **streamer_kwargs)
        else:
            raise ValueError("Either annotation_name or config must be provided")

        self.add_streamer(streamer)
        return self

    def add_annotations_by_category(
        self, category: str, **streamer_kwargs
    ) -> "ConfigurableAnnotationPipeline":
        """Add all annotations from a category"""
        annotation_names = self.registry.list_by_category(category)
        for name in annotation_names:
            self.add_annotation(annotation_name=name, **streamer_kwargs)
        return self

    def add_feature_transformation(
        self,
        transform_func: Callable[[hl.Table], hl.Table],
        description: str = "Custom transformation",
    ):
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
                result = final_result  # Ensure the returned value is the unioned, checkpointed table
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
            loader_func=lambda _: __import__(
                "hvantk.data.dataset", fromlist=["get_dbnsfp_scores_ht"]
            ).get_dbnsfp_scores_ht(),
            feature_mapping={"CADD_phred": "cadd_score", "REVEL_score": "revel_score"},
            metadata={"category": "prediction", "data_type": "scores"},
        ),
        category="prediction",
    )

    # Gene expression
    registry.register(
        AnnotationConfig(
            name="gene_expression",
            source_path="",
            annotation_type="gene",
            loader_func=lambda _: __import__(
                "hvantk.data.dataset", fromlist=["get_gene_expression_ht"]
            ).get_gene_expression_ht(),
            metadata={"category": "expression", "data_type": "levels"},
        ),
        category="expression",
    )

    # Population frequencies
    registry.register(
        AnnotationConfig(
            name="gnomad_frequencies",
            source_path="",
            annotation_type="variant",
            loader_func=lambda _: __import__(
                "hvantk.data.dataset", fromlist=["get_gnomad_af_ht"]
            ).get_gnomad_af_ht(),
            feature_mapping={"AF": "allele_frequency", "AC": "allele_count"},
            metadata={"category": "population", "data_type": "frequencies"},
        ),
        category="population",
    )

    return registry


# Factory functions for easy usage


def create_flexible_pipeline(
    base_streamer: HailDataStreamer, pipeline_name: str = "FlexibleAnnotationPipeline"
) -> ConfigurableAnnotationPipeline:
    """Create a flexible annotation pipeline with built-in registry"""
    registry = create_builtin_registry()
    return ConfigurableAnnotationPipeline(pipeline_name, base_streamer, registry)


def add_custom_annotation(
    pipeline: ConfigurableAnnotationPipeline,
    name: str,
    source_path: str,
    annotation_type: str = "variant",
    **config_kwargs,
) -> ConfigurableAnnotationPipeline:
    """Helper to add custom annotation to existing pipeline"""
    config = AnnotationConfig(
        name=name,
        source_path=source_path,
        annotation_type=annotation_type,
        **config_kwargs,
    )
    return pipeline.add_annotation(config=config)
