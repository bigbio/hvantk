"""Hail Table builder for the Interactome Insider (INSIDER) per-residue BED.

This module owns ``create_interactome_tb``, the canonical builder that turns
the INSIDER ``Whole_Human_Interactome_Interface_hg38.bed`` (UCSC-style BED
with ``track name=<P1>_ppi_<P2>`` headers identifying each PPI) into a Hail
Table keyed by ``interval<locus<rg>>`` with a ``ppi_ids: array<str>`` field
that preserves the PPI identity carried by the track headers. It was
migrated out of :mod:`hvantk.core.builders.table` so that everything
INSIDER-specific (builder, drift probe, tests, fixtures, SKILL) lives under
the plugin folder at :mod:`hvantk.skills.insider`.

The shared helpers ``_create_table_base``, ``_parse_insider_bed_to_temp_tsv``,
and ``_cleanup_temp_file`` intentionally stay in
``hvantk.core.builders.table`` for now -- ``_create_table_base`` and
``_cleanup_temp_file`` are reused across many builders, and the INSIDER BED
helper is small enough that moving it alongside is not worth the churn while
the cleanup-tests still load ``table_builders.py`` directly with stubs.
"""

from __future__ import annotations

import logging

import hail as hl

from hvantk.core.builders.table import (
    _cleanup_temp_file,
    _create_table_base,
    _parse_insider_bed_to_temp_tsv,
)

logger = logging.getLogger(__name__)


def create_interactome_tb(
    input_path: str,
    output_path: str,
    overwrite: bool = False,
    export_tsv: bool = False,
    reference_genome: str = "GRCh38",
) -> "hl.Table":
    """
    Create a Hail Table from an Interactome Insider per-residue BED file.

    The BED is segmented by `track name=<P1>_ppi_<P2>` directives; this builder
    parses those headers and preserves PPI identity as a `ppi_ids: array<str>`
    field per interval. Intervals appearing in multiple PPI tracks are
    aggregated (collected as a sorted, deduplicated array).

    Example usage:
        ht = create_interactome_tb(
            input_path="/path/to/Whole_Human_Interactome_Interface_hg38.bed",
            output_path="/path/to/output.ht"
        )

    Parameters
    ----------
    input_path : str
        Path to the INSIDER BED input file (must contain `track name=...`
        directives to identify PPIs; plain BEDs without tracks produce
        empty output).
    output_path : str
        Path to write the output Hail Table.
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).
    reference_genome : str, optional
        Reference genome to use for parsing intervals (default: "GRCh38").

    Returns
    -------
    hl.Table
        Hail Table keyed by `interval<locus<rg>>` with field
        `ppi_ids: array<str>` carrying the PPI identifiers from track headers.
    """
    tsv_path = None

    def _import() -> "hl.Table":
        nonlocal tsv_path
        tsv_path = _parse_insider_bed_to_temp_tsv(input_path)
        ht = hl.import_table(
            tsv_path,
            types={"start": hl.tint32, "end": hl.tint32},
            min_partitions=4,
        )
        # BED is 0-based half-open; Hail loci are 1-based. Match hl.import_bed's
        # conversion by shifting both endpoints by +1 (so a BED row [100, 200)
        # becomes Hail interval [chr:101, chr:201)).
        ht = ht.annotate(
            interval=hl.locus_interval(
                ht.contig,
                ht.start + 1,
                ht.end + 1,
                reference_genome=reference_genome,
            )
        )
        return ht.select("interval", "ppi_id")

    def _transform(ht: "hl.Table") -> "hl.Table":
        grouped = ht.group_by(ht.interval).aggregate(
            ppi_ids=hl.agg.collect_as_set(ht.ppi_id)
        )
        return grouped.annotate(ppi_ids=hl.sorted(hl.array(grouped.ppi_ids)))

    try:
        return _create_table_base(
            source_name="interactome",
            input_path=input_path,
            output_path=output_path,
            import_func=_import,
            transform_func=_transform,
            overwrite=overwrite,
            export_tsv=export_tsv,
        )
    finally:
        if tsv_path is not None:
            _cleanup_temp_file(tsv_path)


def build_insider_interactome(
    parsed_input,
    ctx,
    *,
    reference_genome: str = "GRCh38",
):
    """Phase B builder — returns an AnnotationTable.

    Parses the INSIDER BED via the existing _parse_insider_bed_to_temp_tsv,
    imports + aggregates, and returns the lazy Hail Table wrapped with
    Provenance. The temp TSV is intentionally NOT cleaned up here (Hail's
    lazy evaluation holds a reference until artifact.save() materializes
    the table); the OS cleans the Hail temp dir at session end.
    """
    from hvantk.core.models import AnnotationTable

    tsv_path = _parse_insider_bed_to_temp_tsv(str(parsed_input))
    try:
        ht = hl.import_table(
            tsv_path,
            types={"start": hl.tint32, "end": hl.tint32},
            min_partitions=4,
        )
        ht = ht.annotate(
            interval=hl.locus_interval(
                ht.contig,
                ht.start + 1,
                ht.end + 1,
                reference_genome=reference_genome,
            )
        )
        ht = ht.select("interval", "ppi_id")
        grouped = ht.group_by(ht.interval).aggregate(
            ppi_ids=hl.agg.collect_as_set(ht.ppi_id)
        )
        grouped = grouped.annotate(ppi_ids=hl.sorted(hl.array(grouped.ppi_ids)))
        return AnnotationTable.from_hail(
            grouped, provenance=ctx.provenance(schema_id="insider-variants-v1")
        )
    except Exception:
        _cleanup_temp_file(tsv_path)
        raise
    # Success: temp file intentionally retained for lazy Hail materialization.
