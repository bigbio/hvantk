"""Hail Table builder for the Ensembl BioMart gene annotation resource.

Owns the Phase B ``build_ensembl_gene_genes`` builder. Imports the Ensembl
BioMart TSV, renames fields, filters to canonical transcripts, groups by
``gene_id``, and wraps with Provenance.
"""
from __future__ import annotations

import logging

import hail as hl

from hvantk.skills.ensembl_gene.shared.constants import ENSEMBL_BIOMART_FIELDS

logger = logging.getLogger(__name__)


def build_ensembl_gene_genes(
    parsed_input,
    ctx,
    **params,
):
    """Phase B builder — returns an AnnotationTable.

    Imports the Ensembl BioMart TSV, renames fields, filters to canonical
    transcripts, groups by gene_id, and wraps with Provenance.

    Parameters
    ----------
    parsed_input : str | Path
        Path to the Ensembl BioMart gene annotation TSV/BGZ file.
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    **params
        Optional: canonical (bool, default True), fields (list of str).
    """
    from hvantk.core.models import AnnotationTable

    canonical = params.get("canonical", True)
    fields = params.get("fields", None)

    ht = hl.import_table(
        paths=str(parsed_input),
        min_partitions=50,
        impute=True,
    )

    logger.info("Replacing field names")
    ht = ht.rename(ENSEMBL_BIOMART_FIELDS)

    if canonical:
        logger.info("Filtering canonical transcripts")
        ht = ht.filter(ht.canonical == "1", keep=True)

    logger.info("Grouping by gene_id")
    ht = (
        ht.group_by(ht.gene_id)
        .aggregate(
            transcript_id=hl.agg.collect_as_set(ht.transcript_id).filter(
                lambda x: x != ""
            ),
            protein_id=hl.agg.collect_as_set(ht.protein_id).filter(
                lambda x: x != ""
            ),
            gene_synonym=hl.agg.collect_as_set(ht.gene_synonym).filter(
                lambda x: x != ""
            ),
            gene_name=hl.agg.take(hl.or_else(ht.gene_name, ""), 1)[0],
            chromosome=hl.agg.take(hl.or_else(ht.chromosome, ""), 1)[0],
            gene_start=hl.agg.take(
                hl.or_else(ht.gene_start, hl.missing(ht.gene_start.dtype)), 1
            )[0],
            gene_end=hl.agg.take(
                hl.or_else(ht.gene_end, hl.missing(ht.gene_end.dtype)), 1
            )[0],
            gene_type=hl.agg.take(hl.or_else(ht.gene_type, ""), 1)[0],
        )
        .key_by("gene_id")
    )

    # 3. Optional field selection
    if fields is not None:
        logger.info("Selecting fields: %s", fields)
        ht = ht.select(*fields)

    # 4. Wrap with provenance
    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="ensembl-gene-v1")
    )
