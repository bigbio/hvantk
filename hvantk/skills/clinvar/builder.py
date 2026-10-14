"""Hail Table builder for the ClinVar VCF resource.

Owns the ``build_clinvar`` builder. Turns the ClinVar VCF into an
``AnnotationTable`` keyed by ``(locus, alleles)``.
"""

from __future__ import annotations

import logging

import hail as hl

from hvantk.core.utils.genome import contig_recoding
from hvantk.core.utils.raw_inputs import resolve_raw_file

logger = logging.getLogger(__name__)

# clinvar declares no lifecycle.parse, so `hvantk reprocess` hands the builder the raw
# download directory (see resolve_raw_file). The downloader writes clinvar.vcf.gz (or
# clinvar_<date>.vcf.gz) and, by default, a .tbi index alongside it; the committed test
# fixture is a .vcf.bgz. These patterns match the VCF under any of those names while
# never matching the .tbi index.
_RAW_FILE_PATTERNS = ("*.vcf.gz", "*.vcf.bgz", "*.vcf")


def build_clinvar(
    parsed_input,
    ctx,
    *,
    reference_genome: str = "GRCh38",
):
    """Plugin builder — returns an AnnotationTable.

    Builds the lazy Hail Table from the VCF and wraps it with provenance.
    The platform's run_builder_for_spec calls artifact.save() to materialize
    the table to disk via ht.write().

    Parameters
    ----------
    parsed_input : str | Path
        Path to the ClinVar VCF (.vcf.gz or .vcf.bgz), or the raw download
        directory produced by ``download_fn`` -- clinvar declares no
        ``lifecycle.parse``, so ``hvantk reprocess`` hands the builder the raw
        directory as-is; ``resolve_raw_file`` picks the one VCF out of it.
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
            path=resolve_raw_file(parsed_input, patterns=_RAW_FILE_PATTERNS),
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
