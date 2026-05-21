"""Hail Table builder for the ClinVar VCF resource.

This module owns ``create_clinvar_tb``, the canonical builder that turns the
ClinVar VCF into a Hail Table keyed by ``(locus, alleles)``. It was migrated
out of :mod:`hvantk.core.builders.table` so that everything ClinVar-specific
(builder, downloader, dataset class, tests, fixtures, SKILL) lives under the
plugin folder at :mod:`hvantk.skills.clinvar`.

The shared helper ``create_table_base`` and the ``contig_recoding`` utility
intentionally stay in their existing modules because they are reused by other
builders.
"""

from __future__ import annotations

import logging

import hail as hl

from hvantk.core.builders.table import create_table_base
from hvantk.core.utils.genome import contig_recoding

logger = logging.getLogger(__name__)


def create_clinvar_tb(
    input_path: str,
    output_path: str,
    overwrite: bool = False,
    export_tsv: bool = False,
    reference_genome: str = "GRCh38",
) -> "hl.Table":
    """
    Create a Hail Table from a ClinVar VCF file keyed by (locus, alleles).

    Example usage:
        ht = create_clinvar_tb(
            input_path="/path/to/clinvar.vcf.gz",
            output_path="/path/to/output.ht"
        )

    Parameters
    ----------
    input_path : str
        Path to the ClinVar VCF input file.
    output_path : str
        Path to write the output Hail Table.
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a flattened TSV version (default: False).
    reference_genome : str, optional
        Reference genome to use for parsing variants (default: "GRCh38").

    Returns
    -------
    hl.Table
        Hail Table keyed by (locus, alleles) with ClinVar annotations.

    Notes
    -----
    Uses force=True for import which may impact performance (single-threaded processing).
    ClinVar's TSV export is special: the row schema contains a deeply nested
    ``info`` struct, so when ``export_tsv=True`` the table is flattened before
    export. This is why we pass ``export_tsv=False`` to ``create_table_base``
    and run ``flatten().export(...)`` ourselves afterwards.
    """
    recode = contig_recoding()

    clinvar_tb = create_table_base(
        source_name="ClinVar",
        input_path=input_path,
        output_path=output_path,
        import_func=lambda: hl.import_vcf(
            path=input_path,
            force=True,
            reference_genome=reference_genome,
            contig_recoding=recode,
            skip_invalid_loci=True,
        ).rows(),
        transform_func=lambda ht: ht.repartition(100).key_by("locus", "alleles"),
        overwrite=overwrite,
        export_tsv=False,  # Handle custom export below
    )

    if export_tsv:
        logger.info(f"Exporting flattened table to {output_path}.tsv.bgz")
        clinvar_tb.flatten().export(f"{output_path}.tsv.bgz")

    return clinvar_tb


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
