"""Hail Table builder for the UniProt post-translational modification (PTM) sites resource.

Owns the Phase B ``build_uniprot_ptm_sites`` builder. Turns the mapped PTM
coordinates TSV (produced by
:mod:`hvantk.algorithms.ptm.pipeline.map_ptm_sites`) into a Hail Table keyed
by locus, wrapped in an ``AnnotationTable`` with Provenance.
"""

from __future__ import annotations

import logging

import hail as hl

logger = logging.getLogger(__name__)


def build_uniprot_ptm_sites(
    parsed_input,
    ctx,
    *,
    reference_genome: str = "GRCh38",
    flanking_codons: int = 5,
    fields=None,
):
    """Phase B builder — returns an AnnotationTable.

    Reads the mapped PTM sites TSV (output of
    hvantk.algorithms.ptm.pipeline.map_ptm_sites) and wraps the lazy Hail Table
    with Provenance.
    """
    from hvantk.core.models import AnnotationTable

    if flanking_codons < 0:
        raise ValueError(f"flanking_codons must be >= 0, got {flanking_codons}")

    ht = hl.import_table(
        paths=str(parsed_input),
        impute=False,
        min_partitions=16,
        types={
            "codon_start": hl.tstr,
            "codon_end": hl.tstr,
            "residue_pos": hl.tstr,
            "n_observations": hl.tstr,
        },
    )

    contig_remap = hl.dict({"MT": "M"})
    ht = ht.annotate(
        _contig=hl.str("chr") + contig_remap.get(ht.chrom, ht.chrom),
    )

    valid_contigs = hl.set(hl.literal(hl.get_reference(reference_genome).contigs))
    ht = ht.filter(valid_contigs.contains(ht._contig))

    ht = ht.annotate(
        locus=hl.locus(
            ht._contig,
            hl.int32(ht.codon_start),
            reference_genome=reference_genome,
        ),
    )

    ht = ht.annotate(
        codon_start=hl.int32(ht.codon_start),
        codon_end=hl.int32(ht.codon_end),
        residue_pos=hl.int32(ht.residue_pos),
        n_observations=hl.int32(ht.n_observations),
    )

    flank_bp = flanking_codons * 3
    ref = hl.get_reference(reference_genome)
    chrom_lengths = hl.dict(hl.literal({c: ref.lengths[c] for c in ref.contigs}))
    ht = ht.annotate(
        flanking_interval=hl.locus_interval(
            ht._contig,
            hl.max(1, ht.codon_start - flank_bp),
            hl.min(chrom_lengths.get(ht._contig), ht.codon_end + flank_bp),
            reference_genome=reference_genome,
            includes_end=True,
        ),
    )
    ht = ht.drop("_contig")
    ht = ht.key_by("locus")

    if fields is not None:
        ht = ht.select(*fields)

    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="uniprot-ptm-sites-v1")
    )
