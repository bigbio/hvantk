# Clinvar Training Set Data Streamer
# Implements streaming processing for generating training sets from Clinvar data

import hail as hl
from typing import Iterator, Optional, Set
from hvantk.data.data_streamer import HailDataStreamer, StreamProcessor
import logging

logger = logging.getLogger(__name__)


class ClinvarDataStreamer(HailDataStreamer):
    """
    Streams Clinvar data for training set generation.
    """

    # Pathogenic labels from Clinvar
    PATHOGENIC_LABELS = [
        "Pathogenic/Likely_pathogenic",
        "Likely_pathogenic",
        "Pathogenic",
    ]

    # Disease-specific labels from Clinvar
    CHD_LABELS = ["Congenital_heart_disease", "Congenital_heart_defect"]

    # Benign labels from Clinvar
    BENIGN_LABELS = ["Benign/Likely_benign", "Likely_benign", "Benign"]

    def __init__(self,
                 clinvar_path: str,
                 chd_genes: Optional[Set[str]] = None,
                 chunk_size: int = 10000):
        super().__init__("ClinvarTrainingSet", chunk_size)
        self.clinvar_path = clinvar_path
        self.chd_genes = chd_genes or set()
        self.clinvar_ht = None

    def setup(self) -> None:
        """Load and prepare Clinvar data"""
        super().setup()
        self.logger.info(f"Loading Clinvar data from {self.clinvar_path}")

        # Load Clinvar VCF
        self.clinvar_ht = hl.import_vcf(self.clinvar_path, reference_genome='GRCh38')

        # Add gene annotation
        self.logger.info("Adding gene column to Clinvar table")
        self.clinvar_ht = self.clinvar_ht.annotate(
            gene=self.clinvar_ht.info.GENEINFO.split("[:]")[0]
        )

        # Add consequence annotation
        self.logger.info("Adding consequence annotations")
        self.clinvar_ht = self.clinvar_ht.annotate(
            Consequence=self.clinvar_ht.info.MC.map(lambda x: x.split("[|]")[1])
        )

    def stream(self) -> Iterator[hl.Table]:
        """
        Stream Clinvar data in chunks for processing.
        """
        if self.clinvar_ht is None:
            raise ValueError("Clinvar data not loaded. Call setup() first.")

        # Filter to non-synonymous variants
        self.logger.info("Filtering to non-synonymous variants")
        filtered_ht = self.clinvar_ht.filter(
            ~self.clinvar_ht.Consequence.any(lambda x: x == "synonymous_variant")
        )

        # Get total count for chunking
        total_variants = filtered_ht.count()
        self.logger.info(f"Processing {total_variants} non-synonymous variants")

        # Process in chunks
        for i in range(0, total_variants, self.chunk_size):
            chunk_ht = filtered_ht.head(min(self.chunk_size, total_variants - i))
            yield self.process_chunk(chunk_ht)

    def process_chunk(self, chunk_ht: hl.Table) -> hl.Table:
        """
        Process a chunk of Clinvar data to generate training labels.
        """
        self.logger.debug(f"Processing chunk with {chunk_ht.count()} variants")

        # Annotate TP/FP labels based on pathogenic/benign annotations
        ts_ann_expr = {
            "is_tp_site": hl.case()
            .when(
                chunk_ht.info.CLNSIG.any(
                    lambda x: hl.set(self.PATHOGENIC_LABELS).contains(x)
                )
                & hl.set(self.chd_genes).contains(chunk_ht.gene),
                True,
            )
            .when(
                chunk_ht.info.CLNDN.any(
                    lambda x: hl.set(self.CHD_LABELS).contains(x)
                ),
                True
            )
            .default(False),

            "is_tn_site": hl.case()
            .when(
                chunk_ht.info.CLNSIG.any(
                    lambda x: hl.set(self.BENIGN_LABELS).contains(x)
                ),
                True,
            )
            .default(False),
        }

        # Apply annotations
        annotated_ht = chunk_ht.annotate(**ts_ann_expr)

        # Filter to variants that are either TP or TN (but not both)
        annotated_ht = annotated_ht.filter(
            annotated_ht.is_tp_site != annotated_ht.is_tn_site
        )

        # Add final RF label
        annotated_ht = annotated_ht.annotate(
            rf_label=hl.case()
            .when(annotated_ht.is_tp_site, "TP")
            .when(annotated_ht.is_tn_site, "TN")
            .or_missing()
        )

        # Filter to defined labels and select relevant columns
        result_ht = annotated_ht.filter(hl.is_defined(annotated_ht.rf_label))
        result_ht = result_ht.select("gene", "rf_label")

        return result_ht


class ClinvarTrainingSetProcessor(StreamProcessor):
    """
    Complete processor for generating Clinvar training sets.
    """

    def __init__(self,
                 clinvar_path: str,
                 chd_genes: Optional[Set[str]] = None,
                 output_dir: str = "./data/training_set"):
        super().__init__("ClinvarTrainingSetGeneration")
        self.output_dir = output_dir

        # Create and add the Clinvar streamer
        clinvar_streamer = ClinvarDataStreamer(
            clinvar_path=clinvar_path,
            chd_genes=chd_genes
        )
        self.add_streamer(clinvar_streamer)

    def process(self, output_path: Optional[str] = None) -> Optional[hl.Table]:
        """
        Process Clinvar data and generate training set.

        Returns:
            Final training set as Hail Table, or None if no data generated
        """
        self.logger.info("Starting Clinvar training set generation")

        # Setup streamers
        for streamer in self.streamers:
            streamer.setup()

        try:
            # Collect all processed chunks
            all_chunks = []
            clinvar_streamer = self.streamers[0]

            for chunk in clinvar_streamer.stream():
                if chunk.count() > 0:  # Only keep non-empty chunks
                    all_chunks.append(chunk)

            if not all_chunks:
                self.logger.warning("No training data generated")
                return None

            # Union all chunks into final table
            self.logger.info(f"Combining {len(all_chunks)} chunks into final training set")
            final_ht = all_chunks[0]
            for chunk in all_chunks[1:]:
                final_ht = final_ht.union(chunk)

            # Save results
            output_path = output_path or f"{self.output_dir}/ts.clinvar.ht"
            self.logger.info(f"Saving training set to {output_path}")
            final_ht = final_ht.checkpoint(output=output_path, overwrite=True)

            # Export as TSV
            final_ht.export(f"{output_path}.tsv")

            self.logger.info(f"Training set generation complete. Final count: {final_ht.count()}")
            return final_ht

        finally:
            # Teardown
            for streamer in reversed(self.streamers):
                streamer.teardown()


def load_chd_gene_set() -> Set[str]:
    """
    Load CHD-associated genes.
    This is a placeholder - implement based on your data source.
    """
    # TODO: Implement actual CHD gene loading logic
    # For now, return a sample set
    return {
        "GATA4", "NKX2-5", "TBX5", "NOTCH1", "CHD7",
        "TBX1", "MYH6", "ACTC1", "MYH7", "TNNT2"
    }


def create_clinvar_training_set_streamer(
    clinvar_path: str,
    output_dir: str = "./data/training_set",
    chd_genes: Optional[Set[str]] = None
) -> ClinvarTrainingSetProcessor:
    """
    Factory function to create a Clinvar training set processor.

    Args:
        clinvar_path: Path to Clinvar VCF file
        output_dir: Output directory for training set
        chd_genes: Set of CHD-associated genes (optional)

    Returns:
        Configured ClinvarTrainingSetProcessor
    """
    if chd_genes is None:
        chd_genes = load_chd_gene_set()

    return ClinvarTrainingSetProcessor(
        clinvar_path=clinvar_path,
        chd_genes=chd_genes,
        output_dir=output_dir
    )
