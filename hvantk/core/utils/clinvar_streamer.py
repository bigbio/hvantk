# Clinvar Training Set Data Streamer
# Implements streaming processing for generating training sets from Clinvar data
#
# Caller contract: pass a pre-built Hail Table via the `table` parameter.
# Build the table with hvantk.skills.clinvar.builder.create_clinvar_tb() first,
# then construct ClinvarDataStreamer(table=ht, ...).  This keeps core/ free of
# any dependency on skills/.

import hail as hl
from typing import Iterator, Optional, Set, Iterable
from hvantk.core.utils.streaming import HailDataStreamer, StreamProcessor
import logging

from hvantk.core.utils.gene_sets import load_gene_set

logger = logging.getLogger(__name__)


class ClinvarDataStreamer(HailDataStreamer):
    """Streams ClinVar data for training set generation using optional gene and disease filters.

    The caller is responsible for building the Hail Table and passing it via
    the ``table`` parameter.  Example::

        from hvantk.skills.clinvar.builder import create_clinvar_tb
        ht = create_clinvar_tb(input_path="/data/clinvar.vcf.bgz",
                               output_path="/tmp/clinvar.ht")
        streamer = ClinvarDataStreamer(table=ht, gene_set={"GATA4"})

    TP criteria:
      1. Pathogenic/Likely pathogenic CLNSIG + gene in provided gene_set (if any)
      2. Disease name (CLNDN) matches any provided disease_terms (case-insensitive)
    """

    PATHOGENIC_LABELS = [
        "Pathogenic/Likely_pathogenic",
        "Likely_pathogenic",
        "Pathogenic",
    ]
    # Removed CHD_LABELS (legacy) – users must supply disease_terms explicitly
    BENIGN_LABELS = ["Benign/Likely_benign", "Likely_benign", "Benign"]

    def __init__(
        self,
        clinvar_path: Optional[str] = None,
        gene_set: Optional[Set[str]] = None,
        disease_terms: Optional[Set[str]] = None,
        chunk_size: int = 10000,
        use_table_builder: bool = True,
        table_output_path: Optional[str] = None,
        overwrite_table: bool = False,
        filter_to_gene_set: bool = False,
        preserve_variant_keys: bool = False,
        table: Optional[hl.Table] = None,
    ):
        clinvar_name = "ClinvarTrainingSet"
        super().__init__(clinvar_name, chunk_size)
        self.gene_set = gene_set or set()
        # Normalize disease terms (case-insensitive, replace spaces with underscores) if provided
        self.disease_terms = disease_terms or set()
        self._normalized_disease_terms = {
            self._normalize_disease_term(d) for d in self.disease_terms
        }
        self.clinvar_path = clinvar_path
        # Accept a pre-built table to decouple core/ from skills/
        self.clinvar_ht = table
        self.use_table_builder = use_table_builder
        self.table_output_path = table_output_path
        self.overwrite_table = overwrite_table
        self.filter_to_gene_set = filter_to_gene_set
        self.preserve_variant_keys = preserve_variant_keys

    @staticmethod
    def _normalize_disease_term(term: str) -> str:
        return term.replace(" ", "_").lower()

    def setup(self) -> None:
        super().setup()
        # If a pre-built table was provided at construction time, use it directly.
        if self.clinvar_ht is not None:
            self.logger.info("Using pre-built ClinVar Hail Table")
        elif self.clinvar_path:
            self.logger.info(f"Loading ClinVar data from {self.clinvar_path}")
            self.clinvar_ht = hl.import_vcf(
                self.clinvar_path, reference_genome="GRCh38"
            ).rows()
        else:
            raise ValueError(
                "ClinvarDataStreamer requires either 'table' (a pre-built hl.Table) "
                "or 'clinvar_path' (a VCF path) to be provided."
            )
        info_fields = self.clinvar_ht.row.dtype["info"].fields
        if "GENEINFO" in info_fields:
            self.clinvar_ht = self.clinvar_ht.annotate(
                gene=hl.if_else(
                    hl.is_defined(self.clinvar_ht.info.GENEINFO)
                    & (hl.len(self.clinvar_ht.info.GENEINFO) > 0),
                    self.clinvar_ht.info.GENEINFO.split(":")[0],
                    hl.missing(hl.tstr),
                )
            )
        else:
            self.logger.warning(
                "info.GENEINFO field absent; gene set filtering may be ineffective"
            )
            self.clinvar_ht = self.clinvar_ht.annotate(gene=hl.missing(hl.tstr))
        if "MC" in info_fields:
            self.clinvar_ht = self.clinvar_ht.annotate(
                Consequence=self.clinvar_ht.info.MC.map(
                    lambda x: hl.if_else(
                        hl.is_defined(x) & (hl.len(x.split("|")) > 1),
                        x.split("|")[1],
                        "",
                    )
                )
            )
        else:
            self.clinvar_ht = self.clinvar_ht.annotate(
                Consequence=hl.empty_array(hl.tstr)
            )

    def stream(self) -> Iterator[hl.Table]:
        if self.clinvar_ht is None:
            raise ValueError("Clinvar data not loaded. Call setup() first.")
        ht = self.clinvar_ht
        if self.filter_to_gene_set and self.gene_set:
            self.logger.info(f"Filtering to {len(self.gene_set)} genes")
            gs = hl.literal(self.gene_set)
            ht = ht.filter(gs.contains(ht.gene))
        filtered_ht = ht.filter(
            hl.is_missing(ht.Consequence)
            | ~ht.Consequence.any(lambda x: x == "synonymous_variant")
        )
        total = filtered_ht.count()
        self.logger.info(f"Processing {total} variants after filtering")
        if total == 0:
            return
        filtered_ht = filtered_ht.add_index("row_idx")
        n_chunks = (total + self.chunk_size - 1) // self.chunk_size
        for i in range(n_chunks):
            start = i * self.chunk_size
            end = min(start + self.chunk_size, total)
            chunk = filtered_ht.filter(
                (filtered_ht.row_idx >= start) & (filtered_ht.row_idx < end)
            ).drop("row_idx")
            yield self.process_chunk(chunk)

    def process_chunk(self, chunk_ht: hl.Table) -> hl.Table:
        # Build disease-based TP condition if disease_terms provided
        if self._normalized_disease_terms:
            disease_set = hl.literal(self._normalized_disease_terms)
            clndn = hl.or_else(chunk_ht.info.CLNDN, "")
            tokens = (
                hl.str(clndn).split(r"\|").map(lambda t: t.replace(" ", "_").lower())
            )
            disease_tp = tokens.any(lambda t: disease_set.contains(t))
        else:
            disease_tp = hl.literal(False)

        # If gene_set is empty, do not filter by gene (treat as True)
        gene_filter = (
            hl.literal(True)
            if not self.gene_set
            else hl.literal(self.gene_set).contains(chunk_ht.gene)
        )
        gene_tp = (
            chunk_ht.info.CLNSIG.any(
                lambda x: hl.set(self.PATHOGENIC_LABELS).contains(x)
            )
            & gene_filter
        )

        ts_ann_expr = {
            "is_tp_site": gene_tp | disease_tp,
            "is_tn_site": hl.case()
            .when(
                chunk_ht.info.CLNSIG.any(
                    lambda x: hl.set(self.BENIGN_LABELS).contains(x)
                ),
                True,
            )
            .default(False),
        }
        annotated = chunk_ht.annotate(**ts_ann_expr)
        annotated = annotated.filter(annotated.is_tp_site != annotated.is_tn_site)
        annotated = annotated.annotate(
            rf_label=hl.case()
            .when(annotated.is_tp_site, "TP")
            .when(annotated.is_tn_site, "TN")
            .or_missing()
        )
        result = annotated.filter(hl.is_defined(annotated.rf_label))
        if not self.preserve_variant_keys:
            result = result.select("gene", "rf_label")
        return result


class ClinvarTrainingSetProcessor(StreamProcessor):
    def __init__(
        self,
        clinvar_path: str,
        gene_set: Optional[Set[str]] = None,
        disease_terms: Optional[Set[str]] = None,
        output_dir: str = "./data/training_set",
        filter_to_gene_set: bool = False,
    ):
        super().__init__("ClinvarTrainingSetGeneration")
        self.output_dir = output_dir
        self.add_streamer(
            ClinvarDataStreamer(
                clinvar_path=clinvar_path,
                gene_set=gene_set,
                disease_terms=disease_terms,
                filter_to_gene_set=filter_to_gene_set,
            )
        )

    def process(self, output_path: Optional[str] = None) -> Optional[hl.Table]:
        self.logger.info("Starting ClinVar training set generation")
        for streamer in self.streamers:
            streamer.setup()
        try:
            chunks = []
            s = self.streamers[0]
            for chunk in s.stream():
                if chunk.count() > 0:
                    chunks.append(chunk)
            if not chunks:
                self.logger.warning("No training data generated")
                return None
            final_ht = chunks[0]
            for c in chunks[1:]:
                final_ht = final_ht.union(c)
            output_path = output_path or f"{self.output_dir}/ts.clinvar.ht"
            final_ht = final_ht.checkpoint(output=output_path, overwrite=True)
            final_ht = final_ht.select(
                gene=final_ht.gene, rf_label=final_ht.rf_label
            )  # Keep variant keys for downstream annotation compatibility
            final_ht.export(f"{output_path}.tsv")
            self.logger.info(
                f"Training set generation complete. Final count: {final_ht.count()}"
            )
            return final_ht
        finally:
            for streamer in reversed(self.streamers):
                streamer.teardown()


def create_clinvar_training_set_streamer(
    clinvar_path: str,
    output_dir: str = "./data/training_set",
    gene_set: Optional[Iterable[str]] = None,
    gene_set_path: Optional[str] = None,
    disease_terms: Optional[Iterable[str]] = None,
    filter_to_gene_set: bool = False,
) -> ClinvarTrainingSetProcessor:
    """Factory to create a ClinVar training set processor.

    gene_set (optional): iterable of gene identifiers.
    disease_terms (optional): iterable of disease names (from CLNDN) to mark as TP.
    No implicit defaults are applied; unspecified sets are empty.
    """
    combined_genes: Set[str] = set()
    if gene_set_path is not None:
        combined_genes |= load_gene_set(gene_set_path)
    if gene_set is not None:
        combined_genes |= {g for g in gene_set}
    disease_terms_set = set(disease_terms) if disease_terms else None
    return ClinvarTrainingSetProcessor(
        clinvar_path=clinvar_path,
        gene_set=combined_genes if combined_genes else None,
        disease_terms=disease_terms_set,
        output_dir=output_dir,
        filter_to_gene_set=filter_to_gene_set,
    )
