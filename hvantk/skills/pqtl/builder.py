"""Hail Table builder for the pQTL (protein quantitative trait loci) resource.

Owns the ``build_pqtl_metrics`` builder. The transform (GTEx
variant-ID parsing, SE derivation via ``|BETA / STAT|``, gene-symbol →
Ensembl-ID mapping via a GeneCatalogStreamer) is implemented directly here.

Shared GTEx variant-ID parsing helpers live in
``hvantk.core.utils.qtl_helpers`` (also used by the eQTL builder).
"""

from __future__ import annotations

import logging
from collections import Counter
from typing import TYPE_CHECKING

import hail as hl

if TYPE_CHECKING:
    from hvantk.core.streamers.gene_catalog import GeneCatalogStreamer

from hvantk.core.utils.qtl_helpers import (
    parse_gtex_variant_id,
    scan_tissue_files,
)
from hvantk.core.utils.table_utils import select_fields

logger = logging.getLogger(__name__)


# HGNC symbol lists searched, in this order, for a symbol that is not an approved
# symbol: a previous symbol was once the gene's approved symbol, while an alias is an
# informal name that unrelated genes can share.
_SYMBOL_LISTS = ("prev_symbols", "alias_symbols")


def _resolve_by_symbol_lists(gene_catalog, symbols):
    """Resolve symbols through previous symbols, then aliases, to Ensembl gene IDs.

    A symbol is looked up among aliases only when no gene has it as a previous
    symbol, and it resolves only when exactly one gene lists it there.

    Returns ``(resolved, ambiguous)``: ``resolved`` maps each resolved symbol to
    ``(list it was found in, Ensembl gene ID)``, and ``ambiguous`` holds the symbols
    that more than one gene lists. A symbol whose one gene has no Ensembl ID is in
    neither.
    """
    matched, ambiguous = {}, set()
    pending = sorted(symbols)
    for field in _SYMBOL_LISTS:
        genes = (
            gene_catalog.get_genes_listing_symbols(pending, field) if pending else {}
        )
        unmatched = []
        for symbol in pending:
            hgnc_ids = genes.get(symbol, ())
            if len(hgnc_ids) == 1:
                matched[symbol] = (field, next(iter(hgnc_ids)))
            elif hgnc_ids:
                ambiguous.add(symbol)
            else:
                unmatched.append(symbol)
        pending = unmatched
    if not matched:
        return {}, ambiguous

    to_ensembl = gene_catalog.map_ids(
        sorted({hgnc_id for _, hgnc_id in matched.values()}),
        source_type="hgnc_id",
        target_type="ensembl_gene_id",
    )
    resolved = {
        symbol: (field, to_ensembl[hgnc_id])
        for symbol, (field, hgnc_id) in matched.items()
        if to_ensembl.get(hgnc_id)
    }
    return resolved, ambiguous


def _find_gene_id_collisions(approved, gene_ids, symbols_by_tissue):
    """Find, tissue by tissue, the source symbols that would share a gene_id.

    ``gene_ids`` maps each mapped symbol to its Ensembl ID, and ``approved`` holds the
    approved symbols among them. Two different source symbols can map to the same
    ``gene_id``: an approved symbol and a previous or alias symbol of the same gene
    (``_resolve_by_symbol_lists`` only checks a symbol against *other* symbols in the
    *same* list, so it cannot see that an approved symbol already claims the gene),
    two previous/alias symbols of one gene, or two approved symbols that HGNC gives
    one Ensembl ID. Rows of one tissue labelled by both would share ``(locus, alleles,
    gene_id)``. Rows of different tissues are told apart by ``tissue``, so a source
    that names a protein by an old symbol in one tissue and by its approved symbol in
    another duplicates nothing, and both stay mapped.

    Within a tissue, a gene_id that several symbols map to goes to the approved symbol
    among them when there is exactly one; every other symbol, or all of them when no
    single approved symbol arbitrates, is demoted in that tissue.

    Returns ``{tissue: demoted symbols}`` for the tissues with a collision. A demoted
    symbol keeps its source symbol as ``gene_id`` in that tissue only.
    """
    demoted: dict[str, set[str]] = {}
    for tissue, present in sorted(symbols_by_tissue.items()):
        by_gene_id: dict[str, list[str]] = {}
        for symbol in sorted(present & gene_ids.keys()):
            by_gene_id.setdefault(gene_ids[symbol], []).append(symbol)
        for gene_id, symbols in by_gene_id.items():
            if len(symbols) < 2:
                continue
            arbiters = [symbol for symbol in symbols if symbol in approved]
            if len(arbiters) == 1:
                losers = [symbol for symbol in symbols if symbol != arbiters[0]]
                logger.warning(
                    "In tissue %s, pQTL gene symbol(s) %s resolve to the same gene "
                    "(%s) as the approved symbol %r and keep their source symbol as "
                    "gene_id there",
                    tissue,
                    ", ".join(losers),
                    gene_id,
                    arbiters[0],
                )
            else:
                losers = symbols
                logger.warning(
                    "In tissue %s, pQTL gene symbols %s all resolve to the same gene "
                    "(%s) with no single approved symbol among them to arbitrate; "
                    "all keep their source symbol as gene_id there",
                    tissue,
                    ", ".join(symbols),
                    gene_id,
                )
            demoted.setdefault(tissue, set()).update(losers)
    return demoted


def _map_gene_symbols(gene_catalog, symbols_by_tissue):
    """Map each tissue's pQTL gene symbols to Ensembl gene IDs and log how they resolved.

    ``symbols_by_tissue`` maps each tissue to the gene symbols its rows carry.

    Approved symbols are mapped first, with ``map_ids``. If the catalog can list the
    genes behind a symbol (``get_genes_listing_symbols``, as the HGNC catalog can),
    each remaining symbol that is not an approved symbol is looked up among previous
    symbols and, only when no gene lists it there, among aliases. A match counts only
    when exactly one gene lists the symbol in that list: a symbol that several genes
    list there is ambiguous and stays unmapped, instead of going to whichever of
    those genes comes first.

    Two source symbols can still map to the same gene_id -- typically an approved
    symbol and a previous symbol of the same gene. Where both label rows of one
    tissue, those rows would share ``(locus, alleles, gene_id)``, so in that tissue
    the gene_id goes to the one approved symbol among them, or to none of them
    (``_find_gene_id_collisions``). In a tissue where it collides with nothing, a
    symbol stays mapped.

    The summary log line counts a symbol as mapped when it keeps its Ensembl ID in at
    least one tissue, and as colliding when it is demoted in every tissue it occurs in.

    Returns ``{(tissue, symbol): Ensembl ID}``. Unmapped and ambiguous symbols, and a
    colliding symbol in its tissue, are left out of it, so they keep their source
    symbol as ``gene_id``. Raises ``ValueError`` when no symbol keeps an Ensembl ID.
    """
    # A missing gene_name arrives as None; it cannot be mapped, sorted or printed.
    named = {
        s for present in symbols_by_tissue.values() for s in present if s is not None
    }
    if any(None in present for present in symbols_by_tissue.values()):
        logger.warning("Some pQTL rows have no gene name; their gene_id stays missing")

    approved = {
        symbol: gene_id
        for symbol, gene_id in gene_catalog.map_ids(
            sorted(named), source_type="gene_symbol", target_type="ensembl_gene_id"
        ).items()
        if gene_id
    }

    resolved, ambiguous = {}, set()
    if not hasattr(gene_catalog, "get_genes_listing_symbols"):
        logger.warning(
            "The gene catalog cannot look up previous symbols or aliases; "
            "only approved gene symbols are mapped"
        )
    else:
        # An approved symbol whose gene has no Ensembl ID still names that gene, so
        # it is not handed to another gene that once had, or shares, the name.
        resolved, ambiguous = _resolve_by_symbol_lists(
            gene_catalog,
            [
                s
                for s in named
                if s not in approved and not gene_catalog.is_canonical(s)
            ],
        )
    gene_ids = {**approved, **{s: gene_id for s, (_, gene_id) in resolved.items()}}

    demoted = _find_gene_id_collisions(approved, gene_ids, symbols_by_tissue)
    mapping = {
        (tissue, symbol): gene_ids[symbol]
        for tissue, present in symbols_by_tissue.items()
        for symbol in present
        if symbol in gene_ids and symbol not in demoted.get(tissue, ())
    }

    kept = {symbol for _, symbol in mapping}
    colliding = gene_ids.keys() - kept
    via = Counter("approved" if s in approved else resolved[s][0] for s in kept)
    unmapped = sorted(named - gene_ids.keys() - ambiguous)
    logger.log(
        logging.WARNING if ambiguous or colliding or unmapped else logging.INFO,
        "Mapped %d of %d pQTL gene symbol(s) to Ensembl gene IDs: %d approved, "
        "%d via a previous symbol, %d via an alias; %d ambiguous, %d colliding, "
        "%d unmapped",
        len(kept),
        len(named),
        via["approved"],
        via["prev_symbols"],
        via["alias_symbols"],
        len(ambiguous),
        len(colliding),
        len(unmapped),
    )
    if ambiguous:
        logger.warning(
            "%d pQTL gene symbol(s) match more than one gene and stay unmapped, "
            "keeping their source symbol as gene_id. Examples: %s",
            len(ambiguous),
            ", ".join(sorted(ambiguous)[:10]),
        )
    if unmapped:
        logger.warning(
            "Could not map %d pQTL gene symbol(s) to Ensembl gene IDs; they keep "
            "their source symbol as gene_id. Examples: %s",
            len(unmapped),
            ", ".join(unmapped[:10]),
        )
    if named and not mapping:
        raise ValueError(
            f"None of the {len(named)} pQTL gene symbols mapped to an Ensembl gene "
            "ID, so the table could not join eQTL tables in cascade analysis. Check "
            "that the gene catalog is the HGNC lookup table built by 'hvantk "
            "reprocess hgnc:lookup' and has ensembl_gene_id values. For a "
            "symbol-keyed table, build without the catalog and pass "
            "--plugin-arg no_gene_map=true."
        )
    return mapping


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
    """Plugin builder — returns an AnnotationTable keyed by
    ``(locus, alleles, gene_id)``.

    For ``source='gtex_fang'``: Fang et al. (2025) allpairs files
    (space-delimited gzip, TMT mass spectrometry, 5 tissues). SE is derived as
    ``|BETA / STAT|`` (Fang files lack an SE column).

    Gene symbols are mapped to Ensembl gene IDs through ``gene_catalog``, normally
    the HGNC lookup table: approved symbols first, then previous symbols, then
    aliases, accepting a previous symbol or alias only when exactly one gene lists
    it in that list (``_map_gene_symbols``). Within one tissue, a gene_id that
    several source symbols map to goes to the one approved symbol among them, or to
    none of them, and the others keep their source symbol as ``gene_id``
    (``_find_gene_id_collisions``), so two source symbols of one tissue never map to
    the same Ensembl ID. Unmapped and ambiguous symbols keep their source symbol as
    ``gene_id``, and so does a colliding symbol in its tissue; the log counts them and
    names each collision. If no symbol keeps an Ensembl ID, the build raises
    ``ValueError``. Mapping is **required** because the cascade
    join uses ``(locus, alleles, gene_id)`` with Ensembl IDs on the eQTL side; raw
    gene symbols would produce zero matches.

    Pass ``no_gene_map=True`` to opt out of mapping for non-cascade use cases
    (the table will be keyed by raw gene symbol and will NOT join with eQTL
    tables in cascade analysis).
    """
    from hvantk.core.models import AnnotationTable
    from hvantk.skills.pqtl.shared.constants import PQTL_SOURCES

    if source not in PQTL_SOURCES:
        raise ValueError(f"Unknown pQTL source: {source!r}. Supported: {PQTL_SOURCES}")
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
        logger.info("Mapping gene symbols → Ensembl IDs via gene catalog")
        symbols_by_tissue = ht.aggregate(
            hl.agg.group_by(ht.tissue, hl.agg.collect_as_set(ht.gene_symbol))
        )
        ensembl_mapping = _map_gene_symbols(gene_catalog, symbols_by_tissue)
        # Typed, because Hail cannot impute the type of an empty mapping (no rows).
        mapping_literal = hl.literal(
            ensembl_mapping, hl.tdict(hl.ttuple(hl.tstr, hl.tstr), hl.tstr)
        )
        ht = ht.annotate(
            gene_id=hl.or_else(
                mapping_literal.get(hl.tuple([ht.tissue, ht.gene_symbol])),
                ht.gene_symbol,
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
        ht = select_fields(ht, fields)

    return AnnotationTable.from_hail(ht, provenance=ctx.provenance(schema_id="pqtl-v1"))
