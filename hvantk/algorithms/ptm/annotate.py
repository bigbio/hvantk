"""Variant-PTM cross-reference annotation.

Annotates a variant Hail Table with PTM site proximity information using
position expansion and locus-based joins. This avoids Hail's interval join
limitations with overlapping intervals (common for adjacent PTM sites).

Also exposes a SYMBOL+chrom annotator (``annotate_variants_by_symbol``) for
variant tables held in pandas, with a gene symbol, chromosome and position
per row rather than a Hail locus.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pandas as pd

from hvantk.algorithms.ptm.constants import PROXIMAL_BP

if TYPE_CHECKING:
    import hail as hl

logger = logging.getLogger(__name__)

#: The genomic positions of a site's codon bases, in transcript order.
CODON_BASE_FIELDS = ("codon_pos1", "codon_pos2", "codon_pos3")


def _require_codon_bases(fields, what: str) -> None:
    """Refuse a sites table that records only ``codon_start``/``codon_end``.

    Such a table cannot say which bases of a codon split by an intron are the codon's,
    and its sites were placed without regard to the isoform they are numbered on.
    """
    missing = [f for f in CODON_BASE_FIELDS if f not in fields]
    if missing:
        raise ValueError(
            f"{what} has no {', '.join(missing)}: it was built by an hvantk version "
            "that recorded only codon_start/codon_end. Rebuild it with "
            "`hvantk ptm build` or `hvantk reprocess uniprot-ptm:sites`."
        )


def annotate_variants_with_ptm(
    variants_ht: hl.Table,
    ptm_ht: hl.Table,
    flanking_codons: int = 5,
) -> hl.Table:
    """Annotate a variant table with PTM site information.

    For each variant, adds:
    - is_ptm_site: bool (variant falls within a PTM codon)
    - is_ptm_proximal: bool (variant within flanking window but not at codon)
    - ptm_types: set<str> (PTM categories at overlapping/proximal sites)
    - ptm_distance: int (approximate distance in residues to nearest PTM site)
    - ptm_evidence: array<struct> (per-site evidence metadata for the nearest
      PTM sites, with fields: source_db, evidence_type, n_observations,
      uniprot_id, gene_symbol, residue_pos, amino_acid, ptm_type)

    A site is its codon's three bases (``codon_pos1``..``codon_pos3``); for a codon
    split by an intron the intron between them is not part of the site. Its window is
    every position within ``flanking_codons * 3`` bp of one of those bases, counted in
    genomic bases.

    The approach expands each PTM site's flanking window to individual genomic
    positions, aggregates by locus, then joins with the variant table. This
    handles overlapping flanking intervals from adjacent PTM sites correctly.

    Parameters
    ----------
    variants_ht : hl.Table
        Variant table keyed by locus (and optionally alleles).
    ptm_ht : hl.Table
        PTM sites table built via ``hvantk reprocess uniprot_ptm:sites``,
        keyed by locus. Must have fields: codon_pos1, codon_pos2, codon_pos3,
        ptm_category. Optional evidence fields:
        source_db, evidence_type, n_observations, uniprot_id, gene_symbol,
        residue_pos, amino_acid, ptm_type.
    flanking_codons : int
        Number of flanking codons for the proximal window (default: 5).

    Returns
    -------
    hl.Table
        Input table with PTM annotation fields added.
    """
    import hail as hl

    ref_genome = variants_ht.locus.dtype.reference_genome.name
    flank_bp = flanking_codons * 3

    logger.info(
        f"Annotating variants with PTM sites "
        f"(flanking_codons={flanking_codons}, flank_bp={flank_bp})"
    )

    _require_codon_bases(set(ptm_ht.row), "The PTM sites table")
    ptm = ptm_ht.annotate(_codon=hl.array([ptm_ht[f] for f in CODON_BASE_FIELDS]))

    # Expand each PTM site to every position within flank_bp of one of its codon's bases
    ptm = ptm.annotate(
        _positions=hl.array(
            hl.set(
                ptm._codon.flatmap(
                    lambda p: hl.range(hl.max(1, p - flank_bp), p + flank_bp + 1)
                )
            )
        )
    )
    ptm = ptm.explode("_positions")

    # Filter out positions beyond contig boundaries
    contig_lengths = hl.literal(hl.get_reference(ref_genome).lengths)
    ptm = ptm.filter(ptm._positions <= contig_lengths.get(ptm.locus.contig, 0))

    # Classify each position: inside codon vs flanking
    ptm = ptm.annotate(
        _in_codon=ptm._codon.contains(ptm._positions),
        _bp_dist=hl.min(ptm._codon.map(lambda p: hl.abs(ptm._positions - p))),
        _join_locus=hl.locus(ptm.locus.contig, ptm._positions, ref_genome),
    )

    # Detect which evidence fields are available in the PTM table
    ptm_fields = set(ptm.row)
    evidence_fields = [
        "source_db",
        "evidence_type",
        "n_observations",
        "uniprot_id",
        "gene_symbol",
        "residue_pos",
        "amino_acid",
        "ptm_type",
    ]
    available_evidence = [f for f in evidence_fields if f in ptm_fields]
    has_evidence = len(available_evidence) > 0

    if has_evidence:
        logger.info(f"Evidence fields available: {available_evidence}")
        # Build a per-row evidence struct from available fields
        ptm = ptm.annotate(
            _evidence=hl.struct(**{f: ptm[f] for f in available_evidence})
        )

    # Aggregate by genomic position (handles multiple PTM sites at same position)
    agg_exprs = dict(
        _any_codon=hl.agg.any(ptm._in_codon),
        ptm_categories=hl.agg.collect_as_set(ptm.ptm_category),
        _min_bp_dist=hl.agg.min(ptm._bp_dist),
    )
    if has_evidence:
        # Collect evidence structs, keeping only the closest PTM sites
        agg_exprs["_evidence_all"] = hl.agg.collect(
            hl.struct(
                _bp_dist=ptm._bp_dist, _in_codon=ptm._in_codon, evidence=ptm._evidence
            )
        )

    ptm_by_pos = ptm.group_by(locus=ptm._join_locus).aggregate(**agg_exprs)

    if has_evidence:
        # Keep only evidence from the nearest PTM site(s) at each position
        ptm_by_pos = ptm_by_pos.annotate(
            _nearest_evidence=ptm_by_pos._evidence_all.filter(
                lambda x: x._bp_dist == ptm_by_pos._min_bp_dist
            ).map(lambda x: x.evidence)
        )
        ptm_by_pos = ptm_by_pos.drop("_evidence_all")

    # Join with variant table
    ann = ptm_by_pos[variants_ht.locus]

    result = variants_ht.annotate(
        is_ptm_site=hl.is_defined(ann) & ann._any_codon,
        is_ptm_proximal=hl.is_defined(ann) & ~ann._any_codon,
        ptm_types=ann.ptm_categories,
        ptm_distance=hl.if_else(
            hl.is_defined(ann),
            (ann._min_bp_dist + 2) // 3,
            hl.missing(hl.tint32),
        ),
        **(
            {
                "ptm_evidence": hl.if_else(
                    hl.is_defined(ann),
                    ann._nearest_evidence,
                    hl.missing(ann._nearest_evidence.dtype),
                )
            }
            if has_evidence
            else {}
        ),
    )

    logger.info("PTM annotation complete")
    return result


def annotate_variants_by_symbol(
    variants_df: pd.DataFrame,
    ptm_df: pd.DataFrame,
    proximal_bp: int = PROXIMAL_BP,
    variant_gene_col: str = "SYMBOL",
    variant_chrom_col: str = "chrom",
    variant_pos_col: str = "pos",
    atlas_gene_col: str = "gene_symbol",
    atlas_chrom_col: str = "chrom",
) -> pd.DataFrame:
    """Per-variant ``is_ptm_site`` / ``is_ptm_proximal`` via SYMBOL + chrom merge.

    Both flags can be True simultaneously (``is_ptm_proximal`` is NOT
    exclusive of ``is_ptm_site``).

    The routine:

    1. Strips ``^chr`` from ``variant_chrom_col`` and ``atlas_chrom_col``.
    2. Drops atlas rows with a null ``atlas_gene_col`` or codon base
       (``codon_pos1``..``codon_pos3``), and casts the bases to ``int``.
    3. Does an inner merge of the variants with the codon bases on
       ``(variant_gene_col, variant_chrom_col) == (atlas_gene_col, atlas_chrom_col)``.
    4. Flags ``is_ptm_site`` when ``pos`` is one of a codon's three bases, and
       ``is_ptm_proximal`` when it lies within ``proximal_bp`` of one. The intron
       inside a codon split by an intron is neither.

    Parameters
    ----------
    variants_df : pandas.DataFrame
        Variants with at minimum ``variant_gene_col``, ``variant_chrom_col``,
        and ``variant_pos_col`` columns.
    ptm_df : pandas.DataFrame
        PTM-site rows from ``hvantk ptm build`` (``ptm_sites_mapped.tsv.bgz``, or
        ``ptm_sites_combined.tsv.bgz`` when PeptideAtlas or CPTAC sites were
        added) with ``atlas_gene_col``, ``atlas_chrom_col`` and ``codon_pos1``..
        ``codon_pos3`` columns; the ``atlas_*`` names refer to this table.
    proximal_bp : int
        Distance (in base pairs) from the nearest codon base within which a variant
        is ``is_ptm_proximal``. Default: :data:`PROXIMAL_BP` (21 bp).
    variant_gene_col, variant_chrom_col, variant_pos_col : str
        Column names in ``variants_df``.
    atlas_gene_col, atlas_chrom_col : str
        Column names in ``ptm_df``; codon columns are always named
        ``codon_pos1``..``codon_pos3`` to match the TSV written by
        ``hvantk ptm build``
        (:func:`~hvantk.algorithms.ptm.pipeline.ptm_build_pipeline_core`).

    Returns
    -------
    pandas.DataFrame
        A copy of ``variants_df`` with ``is_ptm_site`` and ``is_ptm_proximal``
        boolean columns added. Existing columns are preserved unchanged.
    """
    required_variant_cols = {variant_gene_col, variant_chrom_col, variant_pos_col}
    missing = required_variant_cols - set(variants_df.columns)
    if missing:
        raise KeyError(f"variants_df is missing required columns: {sorted(missing)}")
    required_atlas_cols = {atlas_gene_col, atlas_chrom_col}
    missing = required_atlas_cols - set(ptm_df.columns)
    if missing:
        raise KeyError(f"ptm_df is missing required columns: {sorted(missing)}")
    _require_codon_bases(set(ptm_df.columns), "ptm_df")

    out = variants_df.copy()

    # --- Normalize chromosome notation: strip 'chr' prefix on both sides ---
    out[variant_chrom_col] = (
        out[variant_chrom_col].astype(str).str.replace(r"^chr", "", regex=True)
    )

    atlas = ptm_df.copy()
    atlas[atlas_chrom_col] = (
        atlas[atlas_chrom_col].astype(str).str.replace(r"^chr", "", regex=True)
    )

    # One row per codon base, dropping rows with a null key or base.
    atlas = atlas.dropna(subset=[atlas_gene_col, *CODON_BASE_FIELDS])
    bases = atlas.melt(
        id_vars=[atlas_gene_col, atlas_chrom_col],
        value_vars=list(CODON_BASE_FIELDS),
        value_name="_base",
    )[[atlas_gene_col, atlas_chrom_col, "_base"]]
    bases["_base"] = bases["_base"].astype(int)

    left = (
        out[[variant_gene_col, variant_chrom_col, variant_pos_col]]
        .reset_index()
        .rename(columns={"index": "_var_ix"})
    )
    joined = left.merge(
        bases,
        how="inner",
        left_on=[variant_gene_col, variant_chrom_col],
        right_on=[atlas_gene_col, atlas_chrom_col],
    )
    dist = (joined[variant_pos_col] - joined["_base"]).abs()

    # --- is_ptm_site (one of a codon's bases) ---
    out["is_ptm_site"] = out.index.isin(set(joined.loc[dist == 0, "_var_ix"]))

    # --- is_ptm_proximal (within proximal_bp of a codon base) ---
    out["is_ptm_proximal"] = out.index.isin(
        set(joined.loc[dist <= int(proximal_bp), "_var_ix"])
    )

    return out
