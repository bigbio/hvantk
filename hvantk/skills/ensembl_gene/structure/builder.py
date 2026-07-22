"""Builder for ``ensembl-gene:structure`` -- per-gene structural summary from the GTF.

``ensembl-gene:genes`` (the BioMart dataset) carries only coordinates, name and biotype.
This dataset adds the structural covariates that downstream annotation depends on:
CDS length (needed to normalise PTM site counts), coding-exon count, transcript count and
the MANE Select transcript. Gene length is also a mechanical confounder of rare-variant
burden counts, so it must be available as a nuisance covariate rather than omitted.
"""
from __future__ import annotations

import logging
import os

from hvantk.resources.ensembl_release import ENSEMBL_GTF_FILENAME

logger = logging.getLogger(__name__)

SCHEMA_ID = "ensembl-gene-structure-v1"


def _resolve_gtf_path(parsed_input) -> str:
    """Resolve the builder's input to the GTF file itself.

    This dataset declares ``lifecycle.download`` but no ``lifecycle.parse``, so
    ``hvantk reprocess`` hands the builder the raw *directory* (the download wrote the
    GTF inside it), not the file. Sibling plugins survive this because ``hl.import_table``
    tolerates a directory path; this builder reads the GTF with plain Python ``open()``,
    which raises ``IsADirectoryError`` on a directory. So resolve a directory to the known
    committed filename here. A path that already points at a file is returned unchanged,
    which keeps the direct ``build(input_path=<file>)`` calls (tests, snapshots) working.
    """
    path = str(parsed_input)
    if os.path.isdir(path):
        return os.path.join(path, ENSEMBL_GTF_FILENAME)
    return path


def build_ensembl_gene_structure(parsed_input, ctx, **params):
    """Phase B builder -- returns an AnnotationTable keyed on ``gene_id``.

    Parameters
    ----------
    parsed_input : str | Path
        Path to an Ensembl GTF (gzipped or plain).
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    **params
        protein_coding_only : bool, default False
            Keep only ``gene_biotype == "protein_coding"`` rows.

    Returns
    -------
    hvantk.core.models.AnnotationTable
    """
    import hail as hl

    from hvantk.core.models import AnnotationTable
    from hvantk.skills.ensembl_gene.structure.parse import parse_gtf_structure

    df = parse_gtf_structure(_resolve_gtf_path(parsed_input))

    if params.get("protein_coding_only", False):
        df = df[df.gene_biotype == "protein_coding"].reset_index(drop=True)
        logger.info("Filtered to %d protein-coding genes", len(df))

    ht = hl.Table.from_pandas(df, key=["gene_id"])

    return AnnotationTable.from_hail(ht, provenance=ctx.provenance(schema_id=SCHEMA_ID))
