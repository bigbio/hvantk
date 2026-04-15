"""Variant-PTM cross-reference annotation.

Annotates a variant Hail Table with PTM site proximity information using
position expansion and locus-based joins. This avoids Hail's interval join
limitations with overlapping intervals (common for adjacent PTM sites).

Also exposes a pandas-based SYMBOL+chrom annotator
(``annotate_variants_by_symbol``) that reproduces notebook N's Cell 4 / 11
semantics for the CHD case-control workflow.
"""

import hail as hl
import logging

import pandas as pd

from hvantk.ptm.constants import PROXIMAL_BP

logger = logging.getLogger(__name__)


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

    The approach expands each PTM site's flanking window to individual genomic
    positions, aggregates by locus, then joins with the variant table. This
    handles overlapping flanking intervals from adjacent PTM sites correctly.

    Parameters
    ----------
    variants_ht : hl.Table
        Variant table keyed by locus (and optionally alleles).
    ptm_ht : hl.Table
        PTM sites table from create_ptm_sites_tb, keyed by locus. Must have
        fields: codon_start, codon_end, ptm_category. Optional evidence fields:
        source_db, evidence_type, n_observations, uniprot_id, gene_symbol,
        residue_pos, amino_acid, ptm_type.
    flanking_codons : int
        Number of flanking codons for the proximal window (default: 5).

    Returns
    -------
    hl.Table
        Input table with PTM annotation fields added.
    """
    ref_genome = variants_ht.locus.dtype.reference_genome.name
    flank_bp = flanking_codons * 3

    logger.info(
        f"Annotating variants with PTM sites "
        f"(flanking_codons={flanking_codons}, flank_bp={flank_bp})"
    )

    # Cap codon_end for split codons (exon boundary): treat as contiguous 3bp
    ptm = ptm_ht.annotate(_eff_end=hl.min(ptm_ht.codon_end, ptm_ht.codon_start + 2))

    # Expand each PTM site to all positions in its flanking window
    ptm = ptm.annotate(
        _positions=hl.range(
            hl.max(1, ptm.codon_start - flank_bp),
            ptm._eff_end + flank_bp + 1,
        )
    )
    ptm = ptm.explode("_positions")

    # Filter out positions beyond contig boundaries
    contig_lengths = hl.literal(hl.get_reference(ref_genome).lengths)
    ptm = ptm.filter(ptm._positions <= contig_lengths.get(ptm.locus.contig, 0))

    # Classify each position: inside codon vs flanking
    ptm = ptm.annotate(
        _in_codon=(
            (ptm._positions >= ptm.codon_start) & (ptm._positions <= ptm._eff_end)
        ),
        _bp_dist=hl.max(
            0, hl.max(ptm.codon_start - ptm._positions, ptm._positions - ptm._eff_end)
        ),
        _join_locus=hl.locus(ptm.locus.contig, ptm._positions, ref_genome),
    )

    # Detect which evidence fields are available in the PTM table
    ptm_fields = set(ptm.row)
    evidence_fields = [
        "source_db", "evidence_type", "n_observations",
        "uniprot_id", "gene_symbol", "residue_pos", "amino_acid", "ptm_type",
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
            hl.struct(_bp_dist=ptm._bp_dist, _in_codon=ptm._in_codon,
                      evidence=ptm._evidence)
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
        **({
            "ptm_evidence": hl.if_else(
                hl.is_defined(ann),
                ann._nearest_evidence,
                hl.missing(ann._nearest_evidence.dtype),
            )
        } if has_evidence else {}),
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

    Reproduces notebook_n Cell 4 / Cell 11 semantics exactly. Both flags can
    be True simultaneously (``is_ptm_proximal`` is NOT exclusive of
    ``is_ptm_site``).

    The routine:

    1. Strips ``^chr`` from ``variant_chrom_col`` and ``atlas_chrom_col``.
    2. Drops atlas rows with null ``atlas_gene_col``, ``codon_start``, or
       ``codon_end``.
    3. Casts ``codon_start`` / ``codon_end`` to ``int``.
    4. Does an inner merge on ``(variant_gene_col, variant_chrom_col) ==
       (atlas_gene_col, atlas_chrom_col)``, then filters by
       ``codon_start <= pos <= codon_end`` to flag ``is_ptm_site``.
    5. Repeats the merge against the expanded window
       ``[codon_start - proximal_bp, codon_end + proximal_bp]`` to flag
       ``is_ptm_proximal``.

    Parameters
    ----------
    variants_df : pandas.DataFrame
        Variants with at minimum ``variant_gene_col``, ``variant_chrom_col``,
        and ``variant_pos_col`` columns.
    ptm_df : pandas.DataFrame
        PTM atlas rows (e.g. parsed from ``ptm_sites_combined.tsv.bgz``) with
        ``atlas_gene_col``, ``atlas_chrom_col``, ``codon_start``, and
        ``codon_end`` columns.
    proximal_bp : int
        Flank (in base pairs) applied to both ends of the codon interval for
        ``is_ptm_proximal``. Default: :data:`PROXIMAL_BP` (21 bp).
    variant_gene_col, variant_chrom_col, variant_pos_col : str
        Column names in ``variants_df``.
    atlas_gene_col, atlas_chrom_col : str
        Column names in ``ptm_df``; codon columns are always named
        ``codon_start`` / ``codon_end`` to match the Phase-2 atlas TSV.

    Returns
    -------
    pandas.DataFrame
        A copy of ``variants_df`` with ``is_ptm_site`` and ``is_ptm_proximal``
        boolean columns added. Existing columns are preserved unchanged.
    """
    required_variant_cols = {variant_gene_col, variant_chrom_col, variant_pos_col}
    missing = required_variant_cols - set(variants_df.columns)
    if missing:
        raise KeyError(
            f"variants_df is missing required columns: {sorted(missing)}"
        )
    required_atlas_cols = {atlas_gene_col, atlas_chrom_col, "codon_start", "codon_end"}
    missing = required_atlas_cols - set(ptm_df.columns)
    if missing:
        raise KeyError(
            f"ptm_df is missing required columns: {sorted(missing)}"
        )

    out = variants_df.copy()

    # --- Normalize chromosome notation: strip 'chr' prefix on both sides ---
    out[variant_chrom_col] = (
        out[variant_chrom_col].astype(str).str.replace(r"^chr", "", regex=True)
    )

    atlas = ptm_df.copy()
    atlas[atlas_chrom_col] = (
        atlas[atlas_chrom_col].astype(str).str.replace(r"^chr", "", regex=True)
    )

    # Drop atlas rows with null keys / codon bounds, cast codon cols to int.
    atlas = atlas.dropna(subset=[atlas_gene_col, "codon_start", "codon_end"])
    atlas["codon_start"] = atlas["codon_start"].astype(int)
    atlas["codon_end"] = atlas["codon_end"].astype(int)

    # --- Cell 4: is_ptm_site (inside codon interval) ---
    left = (
        out[[variant_gene_col, variant_chrom_col, variant_pos_col]]
        .reset_index()
        .rename(columns={"index": "_var_ix"})
    )
    joined_site = left.merge(
        atlas[[atlas_gene_col, atlas_chrom_col, "codon_start", "codon_end"]],
        how="inner",
        left_on=[variant_gene_col, variant_chrom_col],
        right_on=[atlas_gene_col, atlas_chrom_col],
    )
    hit_site = joined_site[
        (joined_site[variant_pos_col] >= joined_site["codon_start"])
        & (joined_site[variant_pos_col] <= joined_site["codon_end"])
    ]
    ptm_var_ix = set(hit_site["_var_ix"].unique())
    out["is_ptm_site"] = out.index.isin(ptm_var_ix)

    # --- Cell 11: is_ptm_proximal (inside expanded codon interval) ---
    atlas_prox = atlas.copy()
    atlas_prox["prox_start"] = atlas_prox["codon_start"] - int(proximal_bp)
    atlas_prox["prox_end"] = atlas_prox["codon_end"] + int(proximal_bp)
    joined_prox = left.merge(
        atlas_prox[[atlas_gene_col, atlas_chrom_col, "prox_start", "prox_end"]],
        how="inner",
        left_on=[variant_gene_col, variant_chrom_col],
        right_on=[atlas_gene_col, atlas_chrom_col],
    )
    hit_prox = joined_prox[
        (joined_prox[variant_pos_col] >= joined_prox["prox_start"])
        & (joined_prox[variant_pos_col] <= joined_prox["prox_end"])
    ]
    prox_var_ix = set(hit_prox["_var_ix"].unique())
    out["is_ptm_proximal"] = out.index.isin(prox_var_ix)

    return out
