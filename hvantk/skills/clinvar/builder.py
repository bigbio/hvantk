"""Hail Table builder for the ClinVar VCF resource.

Owns the Phase B ``build_clinvar`` builder. Turns the ClinVar VCF into an
``AnnotationTable`` keyed by ``(locus, alleles)``.
"""

from __future__ import annotations

import logging

import hail as hl

from hvantk.core.utils.genome import contig_recoding

logger = logging.getLogger(__name__)


def build_clinvar(
    parsed_input,
    ctx,
    *,
    reference_genome: str = "GRCh38",
):
    """Phase B builder — returns an AnnotationTable.

    Builds the lazy Hail Table from the VCF and wraps it with provenance.
    The platform's run_builder_for_spec calls artifact.save() to materialize
    the table to disk via ht.write().

    Parameters
    ----------
    parsed_input : str | Path
        Path to the ClinVar VCF (.vcf.gz or .vcf.bgz). When the plugin has
        no parse_fn, this is the raw download path from download_fn.
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    reference_genome : str, optional
        Reference genome to use for parsing (default "GRCh38").

    Returns
    -------
    hvantk.core.models.AnnotationTable
        The ClinVar variant table wrapped with Provenance.
    """
    from hvantk.core.models import AnnotationTable

    recode = contig_recoding()
    ht = (
        hl.import_vcf(
            path=str(parsed_input),
            force=True,
            reference_genome=reference_genome,
            contig_recoding=recode,
            skip_invalid_loci=True,
        )
        .rows()
        .repartition(100)
        .key_by("locus", "alleles")
    )
    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="clinvar-variants-v1")
    )
