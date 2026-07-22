"""The gene spine -- the single table every annotation joins onto.

One row per protein-coding gene, keyed on Ensembl ``gene_id``, carrying the identifiers
and structural covariates that later stages need. Identifier reconciliation happens here
and nowhere else, so a downstream source only ever has to map onto ``gene_id``.

``hgnc_id`` is a LEFT join on purpose: a protein-coding Ensembl gene with no HGNC record
must stay in the spine with a missing ``hgnc_id``, not vanish. Dropping it would make the
universe depend on HGNC's coverage, which is the class of silent loss this rebuild exists
to remove.
"""
from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

SPINE_FIELDS = [
    "gene_name",
    "chromosome",
    "gene_start",
    "gene_end",
    "gene_biotype",
    "cds_length",
    "n_coding_exons",
    "n_transcripts",
    "mane_select",
    "hgnc_id",
    "gene_group",
]


def build_spine(genes_ht, structure_ht, hgnc_ht, *, biotype: str = "protein_coding"):
    """Join Ensembl genes, GTF structure and HGNC into the gene spine.

    Parameters
    ----------
    genes_ht : hail.Table
        ``ensembl-gene:genes``, keyed on ``gene_id``.
    structure_ht : hail.Table
        ``ensembl-gene:structure``, keyed on ``gene_id``.
    hgnc_ht : hail.Table
        ``hgnc:lookup``, keyed on ``hgnc_id``, carrying ``ensembl_gene_id``.
    biotype : str
        Gene biotype to keep. Pass ``None`` to keep every biotype.

    Returns
    -------
    hail.Table
        Keyed on ``gene_id``, one row per gene, fields as in ``SPINE_FIELDS``.
    """
    import hail as hl

    ht = genes_ht.join(structure_ht, how="inner")

    if biotype is not None:
        ht = ht.filter(ht.gene_biotype == biotype)

    # Re-key HGNC by its Ensembl ID so the spine can left-join on gene_id. Rows with no
    # Ensembl ID cannot participate and are dropped from the lookup, not from the spine.
    hgnc_by_gene = hgnc_ht.filter(
        hl.is_defined(hgnc_ht.ensembl_gene_id) & (hgnc_ht.ensembl_gene_id != "")
    )

    # HGNC can carry more than one record for the same Ensembl gene. Left-joining a
    # non-unique right side would MULTIPLY spine rows -- and every feature joined on
    # afterwards with them -- so collapse to one record per gene first. The lowest
    # hgnc_id wins, which makes the pick deterministic across rebuilds rather than
    # whichever row the engine happened to emit first.
    hgnc_by_gene = hgnc_by_gene.group_by(
        gene_id=hgnc_by_gene.ensembl_gene_id
    ).aggregate(
        _records=hl.sorted(
            hl.agg.collect(
                hl.struct(
                    hgnc_id=hgnc_by_gene.hgnc_id,
                    gene_group=hgnc_by_gene.gene_group,
                )
            ),
            key=lambda r: r.hgnc_id,
        )
    )
    # Only hgnc_id and gene_group are pulled from HGNC. gene_group is array<str> upstream
    # (a gene can be in several groups) and is carried through unreduced -- a downstream
    # consumer sees the real shape rather than a lossy scalar. mane_select is taken from
    # `structure` (a plain str), NOT from HGNC's own array<str> mane_select field, so this
    # select must not start pulling HGNC's mane_select or it would collide with structure's.
    hgnc_by_gene = hgnc_by_gene.select(
        hgnc_id=hgnc_by_gene._records[0].hgnc_id,
        gene_group=hgnc_by_gene._records[0].gene_group,
    )

    before = ht.count()
    ht = ht.join(hgnc_by_gene, how="left")
    after = ht.count()
    if before != after:
        raise ValueError(
            f"HGNC join changed the row count ({before} -> {after}); the lookup is "
            "not unique per gene_id"
        )

    ht = ht.select(*SPINE_FIELDS)

    logger.info("Spine: %d genes (biotype=%s)", after, biotype)
    return ht


def spine_mapping_rate(spine_ht) -> float:
    """Fraction of spine genes carrying an HGNC record.

    Reported rather than enforced: a gene with no HGNC record is a real Ensembl gene, not
    an error. The rate is a health signal for the identifier layer, and P2's per-source
    gates are where mapping loss becomes fatal.
    """
    import hail as hl

    total = spine_ht.count()
    if total == 0:
        return 1.0
    mapped = spine_ht.aggregate(hl.agg.count_where(hl.is_defined(spine_ht.hgnc_id)))
    return mapped / total
