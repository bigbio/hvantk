"""Variant-PTM cross-reference annotation.

Annotates a variant Hail Table with PTM site proximity information using
position expansion and locus-based joins. This avoids Hail's interval join
limitations with overlapping intervals (common for adjacent PTM sites).
"""

import hail as hl
import logging

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

    The approach expands each PTM site's flanking window to individual genomic
    positions, aggregates by locus, then joins with the variant table. This
    handles overlapping flanking intervals from adjacent PTM sites correctly.

    Parameters
    ----------
    variants_ht : hl.Table
        Variant table keyed by locus (and optionally alleles).
    ptm_ht : hl.Table
        PTM sites table from create_ptm_sites_tb, keyed by locus. Must have
        fields: codon_start, codon_end, ptm_category.
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
    ptm = ptm_ht.annotate(
        _eff_end=hl.min(ptm_ht.codon_end, ptm_ht.codon_start + 2)
    )

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

    # Aggregate by genomic position (handles multiple PTM sites at same position)
    ptm_by_pos = ptm.group_by(locus=ptm._join_locus).aggregate(
        _any_codon=hl.agg.any(ptm._in_codon),
        ptm_categories=hl.agg.collect_as_set(ptm.ptm_category),
        _min_bp_dist=hl.agg.min(ptm._bp_dist),
    )

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
    )

    logger.info("PTM annotation complete")
    return result
