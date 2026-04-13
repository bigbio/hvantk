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
