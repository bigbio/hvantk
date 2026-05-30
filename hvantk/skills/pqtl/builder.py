"""Hail Table builder for the pQTL (protein quantitative trait loci) resource.

Owns the Phase B ``build_pqtl_metrics`` builder. The transform (GTEx
variant-ID parsing, SE derivation via ``|BETA / STAT|``, gene-symbol →
Ensembl-ID mapping via a GeneCatalogStreamer) is implemented directly here.

Shared GTEx variant-ID parsing helpers live in
``hvantk.core.utils.qtl_helpers`` (also used by the eQTL builder).
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import hail as hl

if TYPE_CHECKING:
    from hvantk.core.streamers.gene_catalog import GeneCatalogStreamer

from hvantk.core.utils.qtl_helpers import (
    parse_gtex_variant_id,
    scan_tissue_files,
)

logger = logging.getLogger(__name__)


def _import_gtex_fang(input_path, tissue):
    """Import Fang et al. (2025) pQTL allpairs (space-delimited gzip).

    Columns: ``gene_name SNP CHR BP A1 NMISS BETA STAT P``.
    Rows where ``STAT = 0`` are removed (cannot derive SE).
    """
    files = scan_tissue_files(input_path, [".txt.gz", ".tsv.gz"])

    tables = []
    for fp, tname in files:
        if tissue and tname != tissue:
            continue
        logger.info("Importing pQTL allpairs: %s (tissue: %s)", fp, tname)
        ht_part = hl.import_table(
            fp,
            delimiter=" ",
            force=True,
            types={
                "BETA": hl.tfloat64,
                "STAT": hl.tfloat64,
                "P": hl.tfloat64,
            },
        )
        # STAT = 0 → SE undefined
        ht_part = ht_part.filter(ht_part.STAT != 0.0)
        ht_part = ht_part.select(
            gene_symbol=ht_part.gene_name,
            variant_id=ht_part.SNP,
            beta=ht_part.BETA,
            stat=ht_part.STAT,  # kept for SE derivation in transform
            p_value=ht_part.P,
            tissue=tname,
        )
        tables.append(ht_part)

    if not tables:
        raise FileNotFoundError(f"No pQTL allpairs files matched (tissue={tissue})")
    return tables[0].union(*tables[1:]) if len(tables) > 1 else tables[0]


def build_pqtl_metrics(
    parsed_input,
    ctx,
    *,
    reference_genome: str = "GRCh38",
    source: str = "gtex_fang",
    tissue: str | None = None,
    gene_catalog: GeneCatalogStreamer | None = None,
    no_gene_map: bool = False,
    p_threshold: float | None = None,
    fields: list[str] | None = None,
):
    """Phase B builder — returns an AnnotationTable keyed by
    ``(locus, alleles, gene_id)``.

    For ``source='gtex_fang'``: Fang et al. (2025) allpairs files
    (space-delimited gzip, TMT mass spectrometry, 5 tissues). SE is derived as
    ``|BETA / STAT|`` (Fang files lack an SE column).

    Gene symbols are mapped to Ensembl gene IDs via the HGNC table and
    :class:`~hvantk.core.utils.gene_mapper.GeneMapper`. This is **required**
    because the cascade join uses ``(locus, alleles, gene_id)`` with Ensembl
    IDs on the eQTL side; raw gene symbols would produce zero matches.

    Pass ``no_gene_map=True`` to opt out of mapping for non-cascade use cases
    (the table will be keyed by raw gene symbol and will NOT join with eQTL
    tables in cascade analysis).
    """
    from hvantk.core.models import AnnotationTable
    from hvantk.skills.pqtl.shared.constants import PQTL_SOURCES

    if source not in PQTL_SOURCES:
        raise ValueError(
            f"Unknown pQTL source: {source!r}. Supported: {PQTL_SOURCES}"
        )
    if source != "gtex_fang":
        raise NotImplementedError(
            f"pQTL source {source!r} is not yet implemented. "
            "Only 'gtex_fang' (Fang et al. 2025) is currently supported."
        )

    if gene_catalog is None and not no_gene_map:
        raise ValueError(
            "Ensembl gene mapping is required for cascade-compatible pQTL "
            "tables. Provide --plugin-arg hgnc_ht=<path> (HGNC Hail Table "
            "built by 'hvantk reprocess hgnc:lookup'). If you intentionally "
            "want a symbol-keyed table for non-cascade use, pass "
            "--plugin-arg no_gene_map=true."
        )

    ht = _import_gtex_fang(str(parsed_input), tissue)
    ht = parse_gtex_variant_id(ht, "variant_id", reference_genome)

    # SE = |BETA / STAT| (Fang allpairs lack an SE column).
    ht = ht.annotate(se=hl.abs(ht.beta / ht.stat))
    ht = ht.drop("stat", "variant_id")

    if gene_catalog is not None:
        logger.info(
            "Mapping gene symbols → Ensembl IDs via gene catalog"
        )
        symbols = set(ht.aggregate(hl.agg.collect_as_set(ht.gene_symbol)))
        ensembl_mapping = gene_catalog.map_ids(
            list(symbols), source_type="gene_symbol", target_type="ensembl_gene_id"
        )
        mapping_literal = hl.literal(ensembl_mapping)
        ht = ht.annotate(
            gene_id=hl.or_else(
                mapping_literal.get(ht.gene_symbol), ht.gene_symbol
            )
        )
    else:
        logger.warning(
            "no_gene_map=True: using gene symbols as gene_id. "
            "This table will NOT join with eQTL tables in cascade analysis."
        )
        ht = ht.annotate(gene_id=ht.gene_symbol)

    if p_threshold is not None and p_threshold > 0:
        ht = ht.filter(ht.p_value <= p_threshold)

    ht = ht.annotate(source=source, is_cis=True)
    ht = ht.key_by("locus", "alleles", "gene_id")

    if fields is not None:
        ht = ht.select(*fields)

    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="pqtl-v1")
    )
