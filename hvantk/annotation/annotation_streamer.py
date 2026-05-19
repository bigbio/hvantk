# Multi-Source Annotation Streamer
# Extends the basic streamer to include variant and gene-level annotations

import hail as hl
from typing import Iterator, Optional, Set
from hvantk.core.streamers.base import HailDataStreamer, StreamProcessor
from hvantk.core.streamers.clinvar import ClinvarDataStreamer
import logging

logger = logging.getLogger(__name__)


class AnnotationStreamer(HailDataStreamer):
    """
    Base class for annotation streamers that add features to variant data.
    """

    def __init__(self, name: str, annotation_source: str, chunk_size: int = 10000):
        super().__init__(name, chunk_size)
        self.annotation_source = annotation_source
        self.annotation_data = None

    def load_annotation_data(self) -> hl.Table:
        """
        Load annotation data. Must be implemented by subclasses.
        """
        raise NotImplementedError("Subclasses must implement load_annotation_data")

    def setup(self) -> None:
        """Load annotation data"""
        super().setup()
        self.logger.info(f"Loading annotation data from {self.annotation_source}")
        self.annotation_data = self.load_annotation_data()

    def annotate_chunk(self, chunk: hl.Table) -> hl.Table:
        """
        Add annotations to a chunk. Must be implemented by subclasses.
        """
        raise NotImplementedError("Subclasses must implement annotate_chunk")

    def stream(self) -> Iterator[hl.Table]:
        """
        This is used when the annotator is the first in pipeline.
        Usually annotation streamers are used to process existing chunks.
        """
        raise NotImplementedError(
            "Annotation streamers typically process existing data chunks"
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


class VariantPredictionScoreStreamer(AnnotationStreamer):
    """
    Adds variant prediction scores (CADD, SIFT, PolyPhen, etc.) from dbNSFP.
    """

    def __init__(self, dbnsfp_path: str, chunk_size: int = 10000):
        super().__init__("VariantPredictionScores", dbnsfp_path, chunk_size)
        self.dbnsfp_path = dbnsfp_path

    def load_annotation_data(self) -> hl.Table:
        """Load dbNSFP scores"""
        from hvantk.core.models.dataset import get_dbnsfp_scores_ht

        return get_dbnsfp_scores_ht()

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


class GeneExpressionStreamer(AnnotationStreamer):
    """
    Adds gene expression data from various tissue/cell types.
    """

    def __init__(
        self, expression_path: str, tissue_focus: str = "heart", chunk_size: int = 10000
    ):
        super().__init__("GeneExpression", expression_path, chunk_size)
        self.expression_path = expression_path
        self.tissue_focus = tissue_focus

    def load_annotation_data(self) -> hl.Table:
        """Load gene expression data"""
        from hvantk.core.models.dataset import get_gene_expression_ht, get_hca_ht, get_deg_ht

        # Load multiple expression datasets
        expr_ht = get_gene_expression_ht()  # General expression
        hca_ht = get_hca_ht()  # Human Cell Atlas
        deg_ht = get_deg_ht()  # Differentially expressed genes

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


class GeneConstraintStreamer(AnnotationStreamer):
    """
    Adds gene constraint metrics (pLI, LOEUF, etc.) and evolutionary metrics.
    """

    def __init__(self, chunk_size: int = 10000):
        super().__init__("GeneConstraint", "gnomad_metrics", chunk_size)

    def load_annotation_data(self) -> hl.Table:
        """Load gene constraint data"""
        from hvantk.core.models.dataset import get_gnomad_metrics_ht, get_gevir_ht

        # Load constraint metrics
        gnomad_ht = get_gnomad_metrics_ht()
        gevir_ht = get_gevir_ht()

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


class PopulationFrequencyStreamer(AnnotationStreamer):
    """
    Adds population frequency data from gnomAD.
    """

    def __init__(self, chunk_size: int = 10000):
        super().__init__("PopulationFrequency", "gnomad_af", chunk_size)

    def load_annotation_data(self) -> hl.Table:
        """Load population frequency data"""
        from hvantk.core.models.dataset import get_gnomad_af_ht

        return get_gnomad_af_ht()

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


class EnhancedClinvarTrainingSetProcessor(StreamProcessor):
    """
    Enhanced processor that generates feature-rich training sets by combining
    Clinvar labels with multiple annotation sources.
    """

    def __init__(
        self,
        clinvar_path: str,
        gene_set: Optional[Set[str]] = None,
        output_dir: str = "./data/training_set",
        include_prediction_scores: bool = True,
        include_expression: bool = True,
        include_constraint: bool = True,
        include_population_freq: bool = True,
        tissue_focus: str = "heart",
        filter_to_gene_set: bool = False,
    ):
        super().__init__("EnhancedClinvarTrainingSet")
        self.output_dir = output_dir
        self.tissue_focus = tissue_focus

        # Base ClinVar streamer (reuse generic gene_set)
        clinvar_streamer = ClinvarDataStreamer(
            clinvar_path=clinvar_path,
            gene_set=gene_set,
            filter_to_gene_set=filter_to_gene_set,
            preserve_variant_keys=True,  # Ensure variant-identifying columns are preserved
        )
        self.add_streamer(clinvar_streamer)

        if include_prediction_scores:
            self.add_streamer(VariantPredictionScoreStreamer(""))
        if include_expression:
            self.add_streamer(GeneExpressionStreamer("", tissue_focus))
        if include_constraint:
            self.add_streamer(GeneConstraintStreamer())
        if include_population_freq:
            self.add_streamer(PopulationFrequencyStreamer())

    def process(self, output_path: Optional[str] = None) -> Optional[hl.Table]:
        """
        Process through the entire annotation pipeline.
        """
        self.logger.info("Starting enhanced Clinvar training set generation")

        # Setup all streamers
        for streamer in self.streamers:
            streamer.setup()

        try:
            # Start with Clinvar data
            clinvar_streamer = self.streamers[0]
            annotation_streamers = self.streamers[1:]

            all_chunks = []

            for base_chunk in clinvar_streamer.stream():
                if base_chunk.count() == 0:
                    continue

                # Apply each annotation streamer sequentially
                annotated_chunk = base_chunk
                for annotator in annotation_streamers:
                    annotated_chunk = annotator.process_chunk(annotated_chunk)

                all_chunks.append(annotated_chunk)

            if not all_chunks:
                self.logger.warning("No training data generated")
                return None

            # Union all chunks
            self.logger.info(f"Combining {len(all_chunks)} annotated chunks")
            final_ht = all_chunks[0]
            for chunk in all_chunks[1:]:
                final_ht = final_ht.union(chunk)

            # Trim to only ['gene', 'rf_label'] after all annotation
            final_ht = final_ht.select("gene", "rf_label")

            # Add final feature engineering
            final_ht = self._add_final_features(final_ht)

            # Save results
            output_path = output_path or f"{self.output_dir}/enhanced_ts.clinvar.ht"
            self.logger.info(f"Saving enhanced training set to {output_path}")
            final_ht = final_ht.checkpoint(output=output_path, overwrite=True)

            # Export as TSV with all features
            final_ht.export(f"{output_path}.tsv")

            self.logger.info(
                f"Enhanced training set complete. Features: {len(final_ht.row.dtype.keys())}"
            )
            return final_ht

        finally:
            for streamer in reversed(self.streamers):
                streamer.teardown()

    def _add_final_features(self, ht: hl.Table) -> hl.Table:
        """Add final derived features combining multiple annotation sources"""

        # Composite pathogenicity score
        ht = ht.annotate(
            pathogenicity_score=hl.case()
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
            .or_missing(),
            # Feature completeness score
            feature_completeness=(
                hl.int32(hl.is_defined(ht.combined_deleteriousness))
                + hl.int32(hl.is_defined(ht.conservation_score))
                + hl.int32(hl.is_defined(ht.median_expression))
                + hl.int32(hl.is_defined(ht.constraint_score))
                + hl.int32(hl.is_defined(ht.AF))
            )
            / 5.0,
        )

        return ht


def create_enhanced_clinvar_training_streamer(
    clinvar_path: str,
    output_dir: str = "./data/training_set",
    gene_set: Optional[Set[str]] = None,
    gene_set_path: Optional[str] = None,
    tissue_focus: str = "heart",
    filter_to_gene_set: bool = False,
    **annotation_flags,
) -> "EnhancedClinvarTrainingSetProcessor":
    """Factory creating an enhanced ClinVar training set processor.

    If no gene_set is provided, falls back to sample CHD-oriented set for sandboxing.
    """
    from hvantk.core.utils.gene_sets import load_gene_set, load_sample_chd_gene_set

    if gene_set is None and gene_set_path is None:
        gene_set = load_sample_chd_gene_set()
    else:
        combined: Set[str] = set()
        if gene_set_path is not None:
            combined |= load_gene_set(gene_set_path)
        if gene_set is not None:
            combined |= set(gene_set)
        gene_set = combined

    return EnhancedClinvarTrainingSetProcessor(
        clinvar_path=clinvar_path,
        gene_set=gene_set,
        output_dir=output_dir,
        tissue_focus=tissue_focus,
        filter_to_gene_set=filter_to_gene_set,
        **annotation_flags,
    )
