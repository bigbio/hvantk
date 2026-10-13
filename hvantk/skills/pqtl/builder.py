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


def _drop_gene_id_collisions(mapping, resolved):
    """Demote a previous/alias match whose gene_id is already claimed elsewhere.

    Two different source symbols can resolve to the same Ensembl ``gene_id``: an
    approved symbol and a previous symbol of the same gene (``_resolve_by_symbol_lists``
    only checks a symbol against *other* symbols in the *same* list, so it cannot see
    that an approved symbol already claims the gene), or two previous/alias symbols of
    a gene that has no approved-symbol row in this data at all. Left alone, both would
    merge into ``mapping`` and duplicate ``(locus, alleles, gene_id)`` within a tissue.

    An approved symbol (already in ``mapping``) always keeps its claim; a
    previous/alias symbol that collides with it is demoted, as if it were ambiguous.
    When two previous/alias symbols collide with each other and no approved symbol
    arbitrates, neither keeps the gene_id -- both are demoted.

    Returns ``(resolved, collided)``: ``resolved`` with every demoted symbol removed,
    and the set of demoted symbols (the caller folds this into ``ambiguous``).
    """
    claimed_by_approved = {gene_id: symbol for symbol, gene_id in mapping.items()}
    by_gene_id: dict[str, list[str]] = {}
    for symbol, (_, gene_id) in resolved.items():
        by_gene_id.setdefault(gene_id, []).append(symbol)

    collided: set[str] = set()
    for gene_id, symbols in by_gene_id.items():
        approved_symbol = claimed_by_approved.get(gene_id)
        if approved_symbol is not None:
            collided.update(symbols)
            for symbol in symbols:
                logger.warning(
                    "pQTL gene symbol %r resolves to the same gene (%s) as the "
                    "approved symbol %r already does; %r stays unmapped",
                    symbol,
                    gene_id,
                    approved_symbol,
                    symbol,
                )
        elif len(symbols) > 1:
            collided.update(symbols)
            logger.warning(
                "pQTL gene symbols %s all resolve to the same gene (%s) with no "
                "approved symbol to arbitrate between them; all stay unmapped",
                ", ".join(sorted(symbols)),
                gene_id,
            )

    if not collided:
        return resolved, collided
    return {s: v for s, v in resolved.items() if s not in collided}, collided


def _map_gene_symbols(gene_catalog, symbols):
    """Map pQTL gene symbols to Ensembl gene IDs and log how they resolved.

    Approved symbols are mapped first, with ``map_ids``. If the catalog can list the
    genes behind a symbol (``get_genes_listing_symbols``, as the HGNC catalog can),
    each remaining symbol that is not an approved symbol is looked up among previous
    symbols and, only when no gene lists it there, among aliases. A match counts only
    when exactly one gene lists the symbol in that list: a symbol that several genes
    list there is ambiguous and stays unmapped, instead of going to whichever of
    those genes comes first.

    Two different source symbols can still resolve to the same gene_id -- typically
    an approved symbol and a previous symbol of the same gene -- which would otherwise
    duplicate ``(locus, alleles, gene_id)`` within a tissue. ``_drop_gene_id_collisions``
    closes that: a previous/alias match whose gene_id is already claimed by another
    source symbol is demoted instead (an approved symbol always keeps its claim;
    between two colliding previous/alias symbols with no approved symbol to arbitrate,
    neither keeps it). A demoted symbol is logged and counted separately from
    "ambiguous": it matched exactly one gene in its list, unlike a genuinely ambiguous
    symbol, which matched more than one; the only thing it has in common with ambiguous
    is that it too stays unmapped.

    Returns the symbol -> Ensembl ID mapping. Unmapped, ambiguous and colliding symbols
    are all left out of it, so they keep their source symbol as ``gene_id``. Raises
    ``ValueError`` when no symbol maps.
    """
    # A missing gene_name arrives as None; it cannot be mapped, sorted or printed.
    named = {symbol for symbol in symbols if symbol is not None}
    if len(named) < len(symbols):
        logger.warning("Some pQTL rows have no gene name; their gene_id stays missing")

    mapping = {
        symbol: gene_id
        for symbol, gene_id in gene_catalog.map_ids(
            sorted(named), source_type="gene_symbol", target_type="ensembl_gene_id"
        ).items()
        if gene_id
    }
    n_approved = len(mapping)

    resolved, ambiguous, collided = {}, set(), set()
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
            [s for s in named if s not in mapping and not gene_catalog.is_canonical(s)],
        )
        resolved, collided = _drop_gene_id_collisions(mapping, resolved)
    mapping.update({symbol: gene_id for symbol, (_, gene_id) in resolved.items()})
    via = Counter(field for field, _ in resolved.values())

    unmapped = sorted(named - mapping.keys() - ambiguous - collided)
    logger.log(
        logging.WARNING if ambiguous or collided or unmapped else logging.INFO,
        "Mapped %d of %d pQTL gene symbol(s) to Ensembl gene IDs: %d approved, "
        "%d via a previous symbol, %d via an alias; %d ambiguous, %d colliding, "
        "%d unmapped",
        len(mapping),
        len(named),
        n_approved,
        via["prev_symbols"],
        via["alias_symbols"],
        len(ambiguous),
        len(collided),
        len(unmapped),
    )
    if ambiguous:
        logger.warning(
            "%d pQTL gene symbol(s) match more than one gene and stay unmapped, "
            "keeping their source symbol as gene_id. Examples: %s",
            len(ambiguous),
            ", ".join(sorted(ambiguous)[:10]),
        )
    if collided:
        logger.warning(
            "%d pQTL gene symbol(s) each matched exactly one gene, but a gene_id "
            "another source symbol already claims, and stay unmapped instead, keeping "
            "their source symbol as gene_id (see the per-collision warning above for "
            "which other symbol and gene). Examples: %s",
            len(collided),
            ", ".join(sorted(collided)[:10]),
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
    it in that list (``_map_gene_symbols``). A previous/alias symbol that would
    resolve to a gene_id another source symbol already claims is demoted to
    ambiguous instead -- an approved symbol always keeps its claim, and two
    colliding previous/alias symbols with no approved symbol to arbitrate both stay
    unmapped (``_drop_gene_id_collisions``) -- so two source symbols never produce
    the same ``(locus, alleles, gene_id)`` within a tissue. Unmapped and ambiguous
    symbols keep their source symbol as ``gene_id`` and are counted in the log; if
    no symbol maps, the build raises ``ValueError``. Mapping is **required** because
    the cascade join uses
    ``(locus, alleles, gene_id)`` with Ensembl IDs on the eQTL side; raw gene
    symbols would produce zero matches.

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
        symbols = set(ht.aggregate(hl.agg.collect_as_set(ht.gene_symbol)))
        ensembl_mapping = _map_gene_symbols(gene_catalog, symbols)
        # Typed, because Hail cannot impute the type of an empty mapping (no rows).
        mapping_literal = hl.literal(ensembl_mapping, hl.tdict(hl.tstr, hl.tstr))
        ht = ht.annotate(
            gene_id=hl.or_else(mapping_literal.get(ht.gene_symbol), ht.gene_symbol)
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

    return AnnotationTable.from_hail(ht, provenance=ctx.provenance(schema_id="pqtl-v1"))
