"""Hail Table builder for GTEx (and compatible) cis-eQTL summary statistics.

This module owns ``create_eqtl_tb``, the canonical builder that turns a
directory of per-tissue GTEx v11 (parquet), GTEx v8 (TSV), or eQTLGen (TSV)
cis-eQTL summary statistics into a Hail Table triple-keyed by
``(locus, alleles, gene_id)``. It was migrated out of
:mod:`hvantk.tables.table_builders` so that everything gtex-eqtl-specific
(builder, drift probe, tests, fixtures, SKILL) lives under the plugin folder
at :mod:`hvantk.skills.gtex_eqtl`.

The shared helper ``_create_table_base`` and the source-specific import
helpers (``_import_eqtl_gtex_parquet`` / ``_import_eqtl_gtex_tsv`` /
``_import_eqtl_eqtlgen``) intentionally stay in
``hvantk.tables.table_builders`` because they touch internal parsing
machinery shared across multiple QTL builders.
"""

from __future__ import annotations

import logging
from typing import List, Optional

import hail as hl

from hvantk.tables.table_builders import (
    _create_table_base,
    _import_eqtl_eqtlgen,
    _import_eqtl_gtex_parquet,
    _import_eqtl_gtex_tsv,
    _parse_gtex_variant_id,
    _strip_ensembl_version,
)

logger = logging.getLogger(__name__)


def create_eqtl_tb(
    input_path: str,
    output_path: str,
    reference_genome: str = "GRCh38",
    source: str = "gtex_v11",
    tissue: Optional[str] = None,
    p_threshold: float = 5e-8,
    overwrite: bool = False,
    export_tsv: bool = False,
    fields: Optional[List[str]] = None,
) -> "hl.Table":
    """Build an eQTL Hail Table keyed by ``(locus, alleles, gene_id)``.

    One variant can be an eQTL for multiple genes; the ``gene_id`` key
    prevents information loss and enables correct cascade joins.

    Input can be a single file or a directory of per-tissue files.  Tissue
    name is inferred from the filename prefix before the first dot (e.g.
    ``Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz`` -> ``Brain_Cortex``).

    Supported sources:

    * ``gtex_v11`` -- Parquet signif_pairs (``spark.read.parquet`` ->
      ``hl.Table.from_spark``).
    * ``gtex_v8``  -- TSV ``signif_variant_gene_pairs.txt.gz``
      (``hl.import_table``).
    * ``eqtlgen`` -- TSV cis-eQTLs (single file, different column names).

    Gene-ID version suffixes are stripped for cross-table compatibility
    (``ENSG00000000003.15`` -> ``ENSG00000000003``).

    Parameters
    ----------
    input_path : str
        Single file or directory of per-tissue eQTL files.
    output_path : str
        Output Hail Table path.
    reference_genome : str
        ``GRCh38`` or ``GRCh37``.
    source : str
        Data-source identifier.
    tissue : str, optional
        Restrict import to files matching this tissue name.
    p_threshold : float
        P-value cutoff.  Set to ``0`` to keep all pairs (for coloc).
    overwrite, export_tsv, fields
        Standard builder parameters.
    """
    from hvantk.qtlcascade.constants import EQTL_SOURCES

    if source not in EQTL_SOURCES:
        raise ValueError(f"Unknown eQTL source: {source!r}. Supported: {EQTL_SOURCES}")

    def import_func():
        if source == "gtex_v11":
            return _import_eqtl_gtex_parquet(input_path, tissue, reference_genome)
        if source == "gtex_v8":
            return _import_eqtl_gtex_tsv(input_path, tissue)
        return _import_eqtl_eqtlgen(input_path, reference_genome)

    def transform(ht):
        ht = _parse_gtex_variant_id(ht, "variant_id", reference_genome)
        ht = ht.annotate(gene_id=_strip_ensembl_version(ht.gene_id_raw))
        ht = ht.drop("gene_id_raw", "variant_id")
        if p_threshold > 0:
            ht = ht.filter(ht.p_value <= p_threshold)
        ht = ht.annotate(source=source, is_cis=True)
        ht = ht.key_by("locus", "alleles", "gene_id")
        return ht

    return _create_table_base(
        source_name=f"eQTL ({source})",
        input_path=input_path,
        output_path=output_path,
        import_func=import_func,
        transform_func=transform,
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
