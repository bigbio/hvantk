# Enhanced Training Set Generation Pipeline
# Composes a ClinVar skill data source with algorithm-level annotators.
#
# Moved from hvantk/algorithms/annotation/annotation_streamer.py (issue #121).
# Only tools/ may import both skills/ (ClinvarDataStreamer) and algorithms/
# (VariantPredictionScoreStreamer, etc.) -- that is the correct layer for
# code that crosses the boundary between data sources and annotators.

import hail as hl
from typing import Optional, Set

from hvantk.core.utils.streaming import StreamProcessor
from hvantk.skills.clinvar.pipelines.training_set import ClinvarDataStreamer
from hvantk.algorithms.annotation.annotation_streamer import (
    VariantPredictionScoreStreamer,
    GeneExpressionStreamer,
    GeneConstraintStreamer,
    PopulationFrequencyStreamer,
)

import logging

logger = logging.getLogger(__name__)


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

            # Add final feature engineering. NOTE: do not trim to
            # ['gene', 'rf_label'] here — _add_final_features() consumes the
            # annotation columns (combined_deleteriousness, rarity_score,
            # constraint_score) added by the annotation streamers above, and
            # the TSV export below is meant to carry all features. (A stray
            # select() trim, copied from the basic processor, previously
            # dropped those columns before this step.)
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
