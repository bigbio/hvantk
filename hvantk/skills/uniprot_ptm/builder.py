"""Hail Table builder for the UniProt post-translational modification (PTM) sites resource.

This module owns ``create_ptm_sites_tb``, the canonical builder that turns the
mapped PTM coordinates TSV (produced by
:mod:`hvantk.algorithms.ptm.pipeline.map_ptm_sites`) into a Hail Table keyed by locus. It
was migrated out of :mod:`hvantk.core.builders.table` so that everything
UniProt-PTM-specific (builder, downloader, dataset class, tests, fixtures,
SKILL) lives under the plugin folder at :mod:`hvantk.skills.uniprot_ptm`.

The shared helper ``create_table_base`` intentionally stays in its existing
module because it is reused by other builders.
"""

from __future__ import annotations

import logging
from typing import List, Optional

import hail as hl

from hvantk.core.builders.table import create_table_base

logger = logging.getLogger(__name__)


def create_ptm_sites_tb(
    input_path: str,
    output_path: str,
    reference_genome: str = "GRCh38",
    flanking_codons: int = 5,
    overwrite: bool = False,
    export_tsv: bool = False,
    fields: Optional[List[str]] = None,
) -> "hl.Table":
    """
    Create a Hail Table of PTM sites in genomic coordinates.

    Input is a TSV produced by the PTM coordinate mapper
    (see ``hvantk.algorithms.ptm.pipeline.map_ptm_sites``) with columns: chrom,
    codon_start, codon_end, strand, uniprot_id, gene_symbol, residue_pos,
    amino_acid, ptm_type, ptm_category, source_db, evidence_type,
    n_observations, tissue_type.

    ``tissue_type`` carries sample provenance for sources that distinguish
    it (e.g. CPTAC ``"normal"``/``"tumor"``); curated or bulk-MS sources
    emit an empty string.

    The table is keyed by locus (codon start position), with a
    ``flanking_interval`` field for proximity-based annotation joins.

    Parameters
    ----------
    input_path : str
        Path to the mapped PTM sites TSV file.
    output_path : str
        Path to write the output Hail Table.
    reference_genome : str, optional
        Reference genome (default: "GRCh38").
    flanking_codons : int, optional
        Number of flanking codons for proximal window (default: 5).
    overwrite : bool, optional
        Whether to overwrite existing file (default: False).
    export_tsv : bool, optional
        If True, also export TSV version (default: False).
    fields : list of str, optional
        List of fields to select (default: None, keeps all).

    Returns
    -------
    hl.Table
        The checkpointed Hail Table keyed by locus.
    """
    if flanking_codons < 0:
        raise ValueError(f"flanking_codons must be >= 0, got {flanking_codons}")

    def transform(ht):
        # Remap contig names to match GRCh38 (e.g., MT -> M for chrM)
        contig_remap = hl.dict({"MT": "M"})
        ht = ht.annotate(
            _contig=hl.str("chr") + contig_remap.get(ht.chrom, ht.chrom),
        )

        # Filter to valid contigs in the reference genome
        valid_contigs = hl.set(hl.literal(hl.get_reference(reference_genome).contigs))
        ht = ht.filter(valid_contigs.contains(ht._contig))

        # Parse locus from contig + codon_start
        ht = ht.annotate(
            locus=hl.locus(
                ht._contig,
                hl.int32(ht.codon_start),
                reference_genome=reference_genome,
            ),
        )

        # Cast numeric fields
        ht = ht.annotate(
            codon_start=hl.int32(ht.codon_start),
            codon_end=hl.int32(ht.codon_end),
            residue_pos=hl.int32(ht.residue_pos),
            n_observations=hl.int32(ht.n_observations),
        )

        # Add flanking interval (codon ± flanking_codons * 3 bp)
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
        return ht

    return create_table_base(
        source_name="PTM sites",
        input_path=input_path,
        output_path=output_path,
        import_func=lambda: hl.import_table(
            paths=input_path,
            impute=False,
            min_partitions=16,
            types={
                "codon_start": hl.tstr,
                "codon_end": hl.tstr,
                "residue_pos": hl.tstr,
                "n_observations": hl.tstr,
            },
        ),
        transform_func=transform,
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )


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
