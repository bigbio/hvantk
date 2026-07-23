"""Hail Table builder for the Interactome Insider (INSIDER) per-residue BED.

Owns the Phase B ``build_insider_interactome`` builder. Turns the INSIDER
``Whole_Human_Interactome_Interface_hg38.bed`` (UCSC-style BED with
``track name=<P1>_ppi_<P2>`` headers identifying each PPI) into an
``AnnotationTable`` keyed by ``interval<locus<rg>>`` with a
``ppi_ids: array<str>`` field that preserves the PPI identity carried by
the track headers.

The BED preprocessor lives here (only consumer is this plugin). The shared
``cleanup_temp_file`` lives in ``hvantk.core.utils.hail_helpers``.
"""

from __future__ import annotations

import logging
import re

import hail as hl

from hvantk.core.utils.hail_helpers import cleanup_temp_file

logger = logging.getLogger(__name__)

_TRACK_NAME_RE = re.compile(r'name=([^\s]+)')
_FILE_URI_PREFIX = "file://"


def _normalize_hadoop_path(path: str) -> str:
    """Normalize local file URIs for filesystem APIs."""
    return path[len(_FILE_URI_PREFIX):] if path.startswith(_FILE_URI_PREFIX) else path


def _parse_insider_bed_to_temp_tsv(input_path: str) -> str:
    """Pre-process an Interactome Insider BED into a TSV with ppi_id column.

    The INSIDER BED is segmented by ``track name=<P1>_ppi_<P2> ...``
    directives, each followed by per-residue BED data rows.
    ``hl.import_bed`` silently skips the track headers, dropping PPI identity.
    This helper iterates the BED line-by-line, tracks the current PPI from
    each ``track name=...`` header, and writes a 4-column TSV
    (``contig\\tstart\\tend\\tppi_id``) for downstream ``hl.import_table``.

    Filters applied:
      - ``browser`` lines are ignored.
      - Track headers with no parseable ``name=...`` are skipped
        (``current_ppi_id`` becomes None, so subsequent rows until the next
        valid track are dropped).
      - Degenerate BED rows (``start >= end``, i.e. zero- or negative-length)
        are dropped here — they cannot form a valid interval. Loci that fall
        outside the reference bounds or on a non-reference contig are handled
        separately, downstream in the builder (see
        :func:`build_insider_interactome`).

    Returns the path to a Hail-managed temp file (extension ``tsv``).
    """
    import hailtop.fs as hfs

    out_path = hl.utils.new_temp_file(extension="tsv")
    current_ppi_id: str | None = None
    n_rows_written = 0
    with hfs.open(_normalize_hadoop_path(input_path), "r") as src:
        with hfs.open(_normalize_hadoop_path(out_path), "w") as dst:
            dst.write("contig\tstart\tend\tppi_id\n")
            for line in src:
                if line.startswith("browser"):
                    continue
                if line.startswith("track"):
                    match = _TRACK_NAME_RE.search(line)
                    current_ppi_id = match.group(1) if match else None
                    continue
                if current_ppi_id is None:
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError:
                    continue
                if start >= end:
                    continue
                dst.write(f"{fields[0]}\t{start}\t{end}\t{current_ppi_id}\n")
                n_rows_written += 1
    logger.info(
        "Parsed INSIDER BED %s into %s (%d data rows after filtering)",
        input_path,
        out_path,
        n_rows_written,
    )
    return out_path


def build_insider_interactome(
    parsed_input,
    ctx,
    *,
    reference_genome: str = "GRCh38",
):
    """Phase B builder — returns an AnnotationTable.

    Parses the INSIDER BED via :func:`_parse_insider_bed_to_temp_tsv`,
    imports + aggregates, and returns the lazy Hail Table wrapped with
    Provenance. The temp TSV is intentionally NOT cleaned up on success
    (Hail's lazy evaluation holds a reference until ``artifact.save()``
    materializes the table); the OS cleans the Hail temp dir at session end.
    """
    from hvantk.core.models import AnnotationTable

    tsv_path = _parse_insider_bed_to_temp_tsv(str(parsed_input))
    try:
        ht = hl.import_table(
            tsv_path,
            types={"start": hl.tint32, "end": hl.tint32},
            min_partitions=4,
        )
        # BED is 0-based half-open; with locus_interval's default
        # includes_start=True/includes_end=False, [start+1, end+1) maps the BED
        # span onto Hail's 1-based coordinates. invalid_missing=True sets loci
        # outside the reference bounds or on a non-reference contig to NA rather
        # than raising (the default invalid_missing=False aborts the whole build
        # on a single bad row); the is_defined filter below then drops them,
        # matching the prior builder's skip_invalid_intervals=True intent.
        ht = ht.annotate(
            interval=hl.locus_interval(
                ht.contig,
                ht.start + 1,
                ht.end + 1,
                reference_genome=reference_genome,
                invalid_missing=True,
            )
        )
        ht = ht.filter(hl.is_defined(ht.interval))
        ht = ht.select("interval", "ppi_id")
        grouped = ht.group_by(ht.interval).aggregate(
            ppi_ids=hl.agg.collect_as_set(ht.ppi_id)
        )
        grouped = grouped.annotate(ppi_ids=hl.sorted(hl.array(grouped.ppi_ids)))
        return AnnotationTable.from_hail(
            grouped, provenance=ctx.provenance(schema_id="insider-variants-v1")
        )
    except Exception:
        cleanup_temp_file(tsv_path)
        raise
    # Success: temp file intentionally retained for lazy Hail materialization.
