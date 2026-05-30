# Multi-Source Annotators
# Annotation pipeline steps that add variant and gene-level annotations to
# existing data chunks. These are NOT DataModel streamers (the "Streamer"
# concept is reserved for the Artifact-wrapping ABCs in
# hvantk.core.streamers); they are pipeline-step "Annotators" with a small,
# self-contained base that intentionally does not depend on the legacy
# chunked-IO base in hvantk.core.utils.streaming.

import hail as hl
from abc import ABC, abstractmethod
from typing import Iterator, Optional
from hvantk.core.utils.hail_context import init_hail, hail_initialized
import logging

logger = logging.getLogger(__name__)

# Default chunk size for annotators. Kept local (value 10000) so this module
# does not depend on the legacy streaming module; matches the historical
# default that the annotators inherited from HailDataStreamer.
DEFAULT_CHUNK_SIZE = 10000


class Annotator(ABC):
    """
    Base class for annotators that add features to existing data chunks.

    Self-contained: provides the small surface the annotator family relies on
    (name, chunk_size, logger, idempotent Hail-initializing setup, no-op
    teardown, and the annotation hooks) without inheriting the legacy
    chunked-IO streaming base.
    """

    def __init__(self, name: str, annotation_source: str, chunk_size: int = DEFAULT_CHUNK_SIZE):
        self.name = name
        self.annotation_source = annotation_source
        self.chunk_size = chunk_size
        self.annotation_data = None
        self.logger = logging.getLogger(f"{__name__}.{name}")

    def load_annotation_data(self) -> hl.Table:
        """
        Load annotation data. Must be implemented by subclasses.
        """
        raise NotImplementedError("Subclasses must implement load_annotation_data")

    def setup(self) -> None:
        """Initialize Hail (idempotently) and load annotation data."""
        if not hail_initialized():
            init_hail()
        self.logger.info(f"Loading annotation data from {self.annotation_source}")
        self.annotation_data = self.load_annotation_data()

    def teardown(self) -> None:
        """No-op for global Hail lifecycle (do not stop shared Hail context)."""
        pass

    @abstractmethod
    def annotate_chunk(self, chunk: hl.Table) -> hl.Table:
        """
        Add annotations to a chunk. Must be implemented by subclasses.
        """
        raise NotImplementedError("Subclasses must implement annotate_chunk")

    def stream(self) -> Iterator[hl.Table]:
        """
        This is used when the annotator is the first in pipeline.
        Usually annotators are used to process existing chunks.
        """
        raise NotImplementedError(
            "Annotators typically process existing data chunks"
        )

    def process_chunk(self, chunk: hl.Table) -> hl.Table:
        """Process a chunk by adding annotations"""
        return self.annotate_chunk(chunk)

    # ---------------- Defensive helpers ----------------
    def _safe_gene_annotate(
        self, chunk: hl.Table, gene_field_candidates=("gene", "gene_symbol", "Gene")
    ) -> Optional[hl.Table]:
        """Safely annotate a chunk using a gene field.

        Steps:
        1. Find first existing gene field from candidates.
        2. If none found, log warning and return original chunk.
        3. Attempt to index annotation_data with that field expression.
        4. On any failure (missing annotation_data, missing key, runtime error) log and return original chunk.
        """
        if self.annotation_data is None:
            self.logger.warning(
                f"Annotation data not loaded for {self.name}; skipping gene annotation"
            )
            return chunk
        # Find available gene field using dtype.fields
        gene_fields = getattr(chunk.row.dtype, "fields", None)
        if gene_fields is None:
            gene_fields = list(chunk.row.dtype)
        gene_field = next((f for f in gene_field_candidates if f in gene_fields), None)
        if gene_field is None:
            self.logger.warning(
                f"No gene field ({gene_field_candidates}) found in chunk for {self.name}; skipping"
            )
            return chunk
        try:
            gene_expr = chunk[gene_field]
            ann_row = self.annotation_data[gene_expr]
            # Use annotation_data row dtype fields for valid annotation names
            ann_fields = getattr(self.annotation_data.row.dtype, "fields", None)
            if ann_fields is None:
                ann_fields = list(self.annotation_data.row.dtype)
            annotate_kwargs = {fname: ann_row[fname] for fname in ann_fields}
            annotated = chunk.annotate(**annotate_kwargs)
            return annotated
        except Exception as e:
            self.logger.warning(
                f"Gene annotation failed for field '{gene_field}' in {self.name}: {e}; returning original chunk"
            )
            return chunk


class VariantPredictionScoreAnnotator(Annotator):
    """
    Adds variant prediction scores (CADD, SIFT, PolyPhen, etc.) from dbNSFP.
    """

    def __init__(self, dbnsfp_path: str, chunk_size: int = DEFAULT_CHUNK_SIZE):
        super().__init__("VariantPredictionScores", dbnsfp_path, chunk_size)
        self.dbnsfp_path = dbnsfp_path

    def load_annotation_data(self) -> hl.Table:
        """Load dbNSFP scores"""
        from hvantk.core.io.legacy_artifacts import load_legacy_table

        return load_legacy_table("dbnsfp_scores")

    def annotate_chunk(self, chunk: hl.Table) -> hl.Table:
        """Add variant prediction scores"""
        self.logger.debug(
            f"Adding variant prediction scores to {chunk.count()} variants"
        )

        # Annotate with prediction scores
        annotated = chunk.annotate(**self.annotation_data[chunk.key])

        # Add derived features
        annotated = annotated.annotate(
            # Combined deleteriousness score
            combined_deleteriousness=hl.case()
            .when(
                hl.is_defined(annotated.CADD_phred)
                & hl.is_defined(annotated.REVEL_score),
                (annotated.CADD_phred / 30.0 + annotated.REVEL_score) / 2.0,
            )
            .when(hl.is_defined(annotated.CADD_phred), annotated.CADD_phred / 30.0)
            .when(hl.is_defined(annotated.REVEL_score), annotated.REVEL_score)
            .or_missing(),
            # Conservation score
            conservation_score=hl.case()
            .when(
                hl.is_defined(annotated.phyloP100way_vertebrate),
                annotated.phyloP100way_vertebrate,
            )
            .when(
                hl.is_defined(annotated.phastCons100way_vertebrate),
                annotated.phastCons100way_vertebrate,
            )
            .or_missing(),
            # Functional prediction consensus
            functional_consensus=hl.case()
            .when(
                (annotated.SIFT_pred == "D")
                & (annotated.Polyphen2_HDIV_pred.contains("D")),
                "Damaging",
            )
            .when(
                (annotated.SIFT_pred == "T")
                & (annotated.Polyphen2_HDIV_pred.contains("B")),
                "Tolerated",
            )
            .default("Unknown"),
        )

        return annotated


class GeneExpressionAnnotator(Annotator):
    """
    Adds gene expression data from various tissue/cell types.
    """

    def __init__(
        self, expression_path: str, tissue_focus: str = "heart", chunk_size: int = DEFAULT_CHUNK_SIZE
    ):
        super().__init__("GeneExpression", expression_path, chunk_size)
        self.expression_path = expression_path
        self.tissue_focus = tissue_focus

    def load_annotation_data(self) -> hl.Table:
        """Load gene expression data"""
        from hvantk.core.io.legacy_artifacts import (
            load_legacy_gene_expression_table,
            load_legacy_table,
        )

        # Load multiple expression datasets
        expr_ht = load_legacy_gene_expression_table()  # General expression
        hca_ht = load_legacy_table("hca")  # Human Cell Atlas
        deg_ht = load_legacy_table("deg")  # Differentially expressed genes

        return expr_ht.join(hca_ht, how="outer").join(deg_ht, how="outer")

    def annotate_chunk(self, chunk: hl.Table) -> hl.Table:
        """Add gene expression annotations"""
        # Avoid triggering full count (expensive); rely on lazy logging context.
        self.logger.debug("Adding gene expression data to chunk")

        annotated = self._safe_gene_annotate(chunk)
        if annotated is chunk:
            # Skip derived features if we couldn't annotate
            return chunk

        # Add derived expression features
        annotated = annotated.annotate(
            # Tissue specificity score
            tissue_specificity=hl.case()
            .when(
                hl.is_defined(annotated.heart_expression)
                & hl.is_defined(annotated.median_expression),
                hl.log10(
                    annotated.heart_expression / annotated.median_expression + 0.1
                ),
            )
            .or_missing(),
            # Expression level category
            expression_level=hl.case()
            .when(annotated.median_expression > 10, "High")
            .when(annotated.median_expression > 1, "Medium")
            .when(annotated.median_expression > 0.1, "Low")
            .default("Very_Low"),
            # Heart-specific expression
            heart_enriched=hl.case()
            .when(
                hl.is_defined(annotated.heart_expression)
                & (annotated.heart_expression > 5),
                True,
            )
            .default(False),
        )

        return annotated


class GeneConstraintAnnotator(Annotator):
    """
    Adds gene constraint metrics (pLI, LOEUF, etc.) and evolutionary metrics.
    """

    def __init__(self, chunk_size: int = DEFAULT_CHUNK_SIZE):
        super().__init__("GeneConstraint", "gnomad_metrics", chunk_size)

    def load_annotation_data(self) -> hl.Table:
        """Load gene constraint data"""
        from hvantk.core.io.legacy_artifacts import load_legacy_table

        # Load constraint metrics
        gnomad_ht = load_legacy_table("gnomad_metrics")
        gevir_ht = load_legacy_table("gevir")

        return gnomad_ht.join(gevir_ht, how="outer")

    def annotate_chunk(self, chunk: hl.Table) -> hl.Table:
        """Add gene constraint annotations"""
        self.logger.debug("Adding gene constraint data to chunk")

        annotated = self._safe_gene_annotate(chunk)
        if annotated is chunk:
            return chunk

        # Add derived constraint features
        annotated = annotated.annotate(
            # Haploinsufficiency score
            haploinsufficiency_score=hl.case()
            .when(hl.is_defined(annotated.pLI) & (annotated.pLI >= 0.9), "High")
            .when(hl.is_defined(annotated.pLI) & (annotated.pLI >= 0.5), "Medium")
            .when(hl.is_defined(annotated.pLI), "Low")
            .default("Unknown"),
            # Loss-of-function tolerance
            lof_tolerance=hl.case()
            .when(
                hl.is_defined(annotated.oe_lof_upper)
                & (annotated.oe_lof_upper <= 0.35),
                "Intolerant",
            )
            .when(
                hl.is_defined(annotated.oe_lof_upper)
                & (annotated.oe_lof_upper <= 0.65),
                "Moderate",
            )
            .when(hl.is_defined(annotated.oe_lof_upper), "Tolerant")
            .default("Unknown"),
            # Combined constraint score
            constraint_score=hl.case()
            .when(
                hl.is_defined(annotated.pLI) & hl.is_defined(annotated.oe_lof_upper),
                annotated.pLI * (1 - annotated.oe_lof_upper),
            )
            .when(hl.is_defined(annotated.pLI), annotated.pLI)
            .or_missing(),
        )

        return annotated


class PopulationFrequencyAnnotator(Annotator):
    """
    Adds population frequency data from gnomAD.
    """

    def __init__(self, chunk_size: int = DEFAULT_CHUNK_SIZE):
        super().__init__("PopulationFrequency", "gnomad_af", chunk_size)

    def load_annotation_data(self) -> hl.Table:
        """Load population frequency data"""
        from hvantk.core.io.legacy_artifacts import load_legacy_table

        return load_legacy_table("gnomad_af")

    def annotate_chunk(self, chunk: hl.Table) -> hl.Table:
        """Add population frequency annotations"""
        self.logger.debug(
            f"Adding population frequency data to {chunk.count()} variants"
        )

        # Join on locus and alleles
        annotated = chunk.annotate(**self.annotation_data[chunk.key])

        # Add derived frequency features
        annotated = annotated.annotate(
            # Rarity score (higher = rarer)
            rarity_score=hl.case()
            .when(hl.is_defined(annotated.AF), -hl.log10(annotated.AF + 1e-8))
            .default(8.0),  # Very rare if no frequency data
            # Frequency category
            frequency_category=hl.case()
            .when(annotated.AF > 0.05, "Common")
            .when(annotated.AF > 0.01, "Low_frequency")
            .when(annotated.AF > 0.001, "Rare")
            .when(hl.is_defined(annotated.AF), "Very_rare")
            .default("Novel"),
            # Maximum population frequency
            max_pop_af=hl.max(
                hl.array(
                    [
                        annotated.AF_afr,
                        annotated.AF_amr,
                        annotated.AF_asj,
                        annotated.AF_eas,
                        annotated.AF_fin,
                        annotated.AF_nfe,
                        annotated.AF_oth,
                        annotated.AF_sas,
                    ]
                ).filter(lambda x: hl.is_defined(x))
            ),
        )

        return annotated


