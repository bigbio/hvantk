"""Download + derive lifecycle for ``gevir:metrics`` (issue #386).

The public GeVIR artifact is the paper's Supplementary Table 2 -- an Excel
workbook, not the bgzipped TSV ``build_gevir_metrics`` reads (see
``hvantk/skills/gevir/drift_probe.py`` for the workbook URL/fingerprint and
``hvantk/skills/gevir/SKILL.md`` section 4 for the two-files-not-one background).
This module is the missing extraction step: it downloads the workbook,
parses sheet ``table_2`` with the stdlib (``zipfile`` + ``xml.etree`` --
openpyxl is deliberately not a dependency), and writes the derived TSV the
builder expects, BGZF-compressed.

Layout (frozen -- this is a published supplementary table, not a live feed):
header on row 2, data on rows 3-19363 (19,361 genes). Columns read: A
``gnomad_gene_name``, B ``gene_id``, C ``canonical_transcript``, D
``approved_symbol`` (read only to resolve the restoration below; never
written to the output TSV), F ``gevir_percentile`` (written as ``gevir_pct``),
H ``virlof_percentile`` (written as ``virlof_pct``). The header row's text is
validated before any extraction so a layout change fails loudly rather than
silently misreading a shifted column.

``gevir_pct``/``virlof_pct`` are written at FULL precision: the cell's raw
``<v>`` text, verbatim, with no rounding/reformatting. (The pre-#386 manually
derived TSV rounded to Excel's ~11-character text-export width; this
downloader does not reproduce that rounding.)

Gene-name restoration
----------------------
29 ``gnomad_gene_name`` cells in the live workbook are not strings: Excel's
autocorrect silently reinterpreted the original gene symbol as a date, so the
cell holds a 1900-date-system serial instead of text. This is not a MARCH/SEPT-
only phenomenon -- empirically (live workbook checked 2026-10-13, cross-
referenced against HGNC's REST API and the Ensembl GRCh37 archive) the 29
cells split into five distinct symbol families:

  * MARCHF1-11 (today's HGNC symbol) <- originally ``MARCH<n>``       (23 cells... see below)
  * SEPTIN1-14  (today's HGNC symbol) <- originally ``SEPT<n>``
  * MTARC1/MTARC2                     <- originally ``MARC<n>`` (NOT ``MARCH<n>`` --
    "MARC" alone is close enough to the month name for Excel's parser to
    collide with it too, at a different symbol than the MARCH ring-finger
    family)
  * SELENOF (no numbered suffix)      <- originally ``SEP15`` (NOT ``SEPT15``;
    the other classic Excel-gene-mangling example from Ziemann et al. 2016,
    alongside MARCH/SEPT)
  * DELEC1 ("deleted in esophageal cancer 1"), approved_symbol blank in this
    row <- originally ``DEC1``, confirmed via HGNC's live record
    (prev_symbol=["DEC1"], rest.genenames.org/fetch/symbol/DELEC1)
  * one GRCh37-only retired duplicate gene ID for the septin-1 locus,
    approved_symbol also blank in this row <- originally ``SEPT1``,
    confirmed via the Ensembl GRCh37 archive (display_name "SEPT1" for
    ENSG00000270466, grch37.rest.ensembl.org)

Decoding the serial alone (month + trailing two digits of the decoded year --
the day is uniformly 1; Excel parsed the typed text as a month+2-digit-year
pattern, not day+month, which is also why the pre-existing TSV already showed
these as "Mar 06"-style strings rather than raw numbers) cannot by itself
distinguish MARCH from MARC, or SEPT from SEP -- both pairs parse to the
identical date. Restoration therefore primarily reads the row's OWN
``approved_symbol`` (column D) and mechanically reverses one of three known,
systematic HGNC rename patterns (``_FAMILY_PATTERNS`` below); the decoded
(month, year%100) is still computed and cross-checked against the number
recovered from ``approved_symbol`` as a consistency guard. The two rows whose
``approved_symbol`` is blank cannot be resolved from in-row data at all, so
they are pinned via ``_GENE_ID_OVERRIDES``, each backed by the external
citation in its comment -- not a guess. Any numeric/date-typed name cell that
matches none of these (a sixth family, or a genuinely new one in a future
revision of the workbook) raises a clear, diagnostic error rather than
emitting an unverified symbol, per the explicit "fail loudly" requirement in
issue #386.
"""

from __future__ import annotations

import logging
import os
import re
import tempfile
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime, timedelta
from pathlib import Path

import click

from hvantk.core.utils.bgzf import BgzfWriter
from hvantk.core.utils.file_utils import download_file
from hvantk.skills.gevir.drift_probe import GEVIR_SUPPLEMENTARY_URL

logger = logging.getLogger(__name__)

# Canonical raw-file name, matching the catalog entry's declared file
# (hvantk/skills/gevir/catalog/datasets.json, accession GeVIR_v1.0).
GEVIR_RAW_FILENAME = "gevir_metrics_pmid31873297.tsv.bgz"

# Output TSV column order -- matches the pre-existing derived TSV's header
# (hvantk/tests/testdata/raw/gevir/gevir_metrics_pmid31873297.tsv.bgz) and the
# builder's field names (hvantk/skills/gevir/builder.py).
_OUTPUT_HEADER = (
    "gnomad_gene_name",
    "gene_id",
    "canonical_transcript",
    "gevir_pct",
    "virlof_pct",
)

_MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_NS = {"m": _MAIN_NS}
_R_ID_ATTR = f"{{{_REL_NS}}}id"
_COL_RE = re.compile(r"[A-Z]+")

_SHEET_NAME = "table_2"
_HEADER_ROW_NUM = 2
_FIRST_DATA_ROW_NUM = 3
_EXPECTED_DATA_ROWS = 19361  # rows 3-19363 inclusive

# Row-2 header text for every column this module reads. Column D
# (approved_symbol) is read only to resolve name restoration below and is
# never written to the output TSV.
_EXPECTED_HEADERS = {
    "A": "gnomad_gene_name",
    "B": "gene_id",
    "C": "canonical_transcript",
    "D": "approved_symbol",
    "F": "gevir_percentile",
    "H": "virlof_percentile",
}

# Excel's 1900 date system: day-1 of this epoch is serial 1. Using Dec-30
# (not Dec-31) absorbs the fictitious Feb-29-1900 leap-year bug for any
# serial >= 61, which covers every serial in this workbook by a wide margin.
_EXCEL_EPOCH = datetime(1899, 12, 30)


# ---------------------------------------------------------------------------
# Gene-name restoration
# ---------------------------------------------------------------------------

# Approved-symbol patterns for the three systematic, numbered HGNC rename
# families this workbook's gnomad_gene_name column collides with. Each
# recovers the pre-rename symbol mechanically from the ROW'S OWN
# approved_symbol (column D) -- not from the decoded date, which cannot tell
# MARCH from MARC or SEPT from SEP apart (both members of a pair parse to the
# same date). (expected_month, approved_pattern, original_template).
_FAMILY_PATTERNS = (
    (3, re.compile(r"^MARCHF(\d+)$"), "MARCH{}"),  # MARCHF<n> -> MARCH<n>
    (9, re.compile(r"^SEPTIN(\d+)$"), "SEPT{}"),  # SEPTIN<n> -> SEPT<n>
    (3, re.compile(r"^MTARC(\d+)$"), "MARC{}"),  # MTARC<n>  -> MARC<n>
)

# SELENOF carries no numbered suffix, so it cannot join _FAMILY_PATTERNS.
_SELENOF_APPROVED_SYMBOL = "SELENOF"
_SELENOF_ORIGINAL = "SEP15"
_SELENOF_EXPECTED_MONTH = 9
_SELENOF_EXPECTED_YY = 15

# The two rows whose approved_symbol is blank in the live workbook, so they
# cannot be resolved from in-row data at all. Each is pinned to an external
# citation checked 2026-10-13, not guessed from the serial:
#   - ENSG00000173077: HGNC's live record for DELEC1 ("deleted in esophageal
#     cancer 1") lists prev_symbol=["DEC1"] --
#     rest.genenames.org/fetch/symbol/DELEC1.
#   - ENSG00000270466: a GRCh37-only, now-retired duplicate Ensembl gene ID
#     for the septin-1 locus (Ensembl's current, GRCh38 REST API returns "ID
#     not found" for it). The GRCh37 archive reports display_name "SEPT1" --
#     grch37.rest.ensembl.org/lookup/id/ENSG00000270466 -- matching the
#     *other*, current SEPT1 row (ENSG00000180096, approved_symbol SEPTIN1
#     today).
# Each override's (month, yy) is asserted below against the live serial as a
# staleness guard: if a future workbook revision ever changes this row's
# serial, the mismatch raises instead of silently keeping a now-wrong name.
_GENE_ID_OVERRIDES = {
    "ENSG00000173077": ("DEC1", 12, 1),
    "ENSG00000270466": ("SEPT1", 9, 1),
}


def _restore_gene_name(
    *, is_string: bool, raw_value, gene_id: str, approved_symbol, cell_ref: str
) -> str:
    """Return the restored ``gnomad_gene_name`` for one row's column-A cell.

    ``is_string``/``raw_value`` come from :func:`_cell_text` for column A.
    Ordinary cells (``is_string`` True) pass through unchanged. A numeric/
    date-typed cell is resolved via the gene-id override table, the SELENOF
    special case, or one of ``_FAMILY_PATTERNS`` applied to this row's own
    ``approved_symbol`` -- in that order. Anything else raises ``ValueError``
    naming the cell, the gene id, the decoded date, and the approved_symbol
    seen, rather than emitting an unverified guess.
    """
    if is_string:
        if not raw_value:
            raise ValueError(f"{cell_ref} ({gene_id}): empty gnomad_gene_name cell")
        return raw_value

    if raw_value is None:
        raise ValueError(
            f"{cell_ref} ({gene_id}): numeric gnomad_gene_name cell has no stored value"
        )

    serial = float(raw_value)
    decoded = _EXCEL_EPOCH + timedelta(days=serial)
    yy = decoded.year % 100

    override = _GENE_ID_OVERRIDES.get(gene_id)
    if override is not None:
        name, expected_month, expected_yy = override
        if decoded.month != expected_month or decoded.day != 1 or yy != expected_yy:
            raise ValueError(
                f"{cell_ref} ({gene_id}): pinned override {name!r} expects a serial "
                f"decoding to month={expected_month}, yy={expected_yy:02d}, but this "
                f"workbook's serial={serial!r} decodes to {decoded.date().isoformat()} "
                "-- the override may be stale; refusing to apply it blindly."
            )
        return name

    if approved_symbol == _SELENOF_APPROVED_SYMBOL:
        if (
            decoded.month != _SELENOF_EXPECTED_MONTH
            or decoded.day != 1
            or yy != _SELENOF_EXPECTED_YY
        ):
            raise ValueError(
                f"{cell_ref} ({gene_id}): approved_symbol is SELENOF but the serial "
                f"decodes to {decoded.date().isoformat()}, not the expected "
                f"month={_SELENOF_EXPECTED_MONTH}, yy={_SELENOF_EXPECTED_YY:02d}"
            )
        return _SELENOF_ORIGINAL

    if approved_symbol:
        for expected_month, pattern, template in _FAMILY_PATTERNS:
            m = pattern.match(approved_symbol)
            if not m:
                continue
            number = int(m.group(1))
            if decoded.month != expected_month or decoded.day != 1 or yy != number:
                raise ValueError(
                    f"{cell_ref} ({gene_id}): approved_symbol {approved_symbol!r} implies "
                    f"number {number}, but the serial decodes to "
                    f"{decoded.date().isoformat()} (expected month={expected_month}, "
                    f"yy={number:02d}) -- refusing to emit a possibly-wrong symbol."
                )
            return template.format(number)

    raise ValueError(
        f"{cell_ref} ({gene_id}): numeric/date-typed gnomad_gene_name cell "
        f"(serial={raw_value!r}, decodes to {decoded.date().isoformat()}) matches none "
        f"of the known Excel-date-mangled gene families (approved_symbol="
        f"{approved_symbol!r}); refusing to guess a symbol. See this module's "
        "docstring for the five families handled today."
    )


# ---------------------------------------------------------------------------
# XLSX parsing (stdlib only: zipfile + xml.etree)
# ---------------------------------------------------------------------------


def _find_sheet_xml(zf: zipfile.ZipFile, sheet_name: str) -> bytes:
    """Resolve *sheet_name* (a workbook tab name) to its worksheet XML bytes."""
    workbook = ET.fromstring(zf.read("xl/workbook.xml"))
    sheets_el = workbook.find("m:sheets", _NS)
    rid = None
    if sheets_el is not None:
        for sheet in sheets_el.findall("m:sheet", _NS):
            if sheet.get("name") == sheet_name:
                rid = sheet.get(_R_ID_ATTR)
                break
    if rid is None:
        raise ValueError(f"xl/workbook.xml has no sheet named {sheet_name!r}")

    rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    target = None
    for rel in rels:
        if rel.get("Id") == rid:
            target = rel.get("Target")
            break
    if target is None:
        raise ValueError(f"xl/_rels/workbook.xml.rels has no relationship {rid!r}")

    return zf.read(f"xl/{target}")


def _load_shared_strings(zf: zipfile.ZipFile) -> list:
    """Return the shared-string table as a list indexed by shared-string id."""
    try:
        data = zf.read("xl/sharedStrings.xml")
    except KeyError:
        return []
    root = ET.fromstring(data)
    out = []
    for si in root.findall("m:si", _NS):
        texts = si.findall(".//m:t", _NS)
        out.append("".join(t.text or "" for t in texts))
    return out


def _cell_text(cell, shared_strings: list):
    """Return ``(is_string, text)`` for one ``<c>`` element, or ``(True, None)``.

    Handles the three cell-content encodings OOXML uses: ``t="s"`` (shared
    string, looked up by index), ``t="inlineStr"`` (text inline under
    ``<is>``), ``t="str"`` (cached formula string, inline under ``<v>``), and
    numeric (``t`` absent or ``"n"``, value inline under ``<v>``) -- the last
    is what a date-serial cell looks like, which is exactly the case this
    module needs to detect (``is_string=False``).
    """
    if cell is None:
        return True, None

    t = cell.get("t")
    if t == "s":
        v = cell.find("m:v", _NS)
        if v is None or v.text is None:
            return True, None
        return True, shared_strings[int(v.text)]
    if t == "inlineStr":
        is_el = cell.find("m:is", _NS)
        if is_el is None:
            return True, None
        texts = is_el.findall(".//m:t", _NS)
        return True, "".join(x.text or "" for x in texts)
    if t == "str":
        v = cell.find("m:v", _NS)
        return True, (v.text if v is not None else None)

    # Numeric: t is absent or "n".
    v = cell.find("m:v", _NS)
    return False, (v.text if v is not None else None)


def _row_cells(row_el) -> dict:
    """Map column letter -> ``<c>`` element for one ``<row>``."""
    cells = {}
    for c in row_el.findall("m:c", _NS):
        ref = c.get("r")
        if not ref:
            continue
        cells[_COL_RE.match(ref).group(0)] = c
    return cells


def extract_gevir_table_2(
    xlsx_path: str, *, expected_data_rows: int | None = _EXPECTED_DATA_ROWS
) -> list:
    """Parse sheet ``table_2`` of the GeVIR supplementary workbook.

    Returns a list of ``(gnomad_gene_name, gene_id, canonical_transcript,
    gevir_pct, virlof_pct)`` tuples (the last two as raw, unrounded ``<v>``
    text), in row order, ready to write as the output TSV's data rows.

    Validates the row-2 header text for every column read (see
    ``_EXPECTED_HEADERS``) before extracting any data, and (when
    *expected_data_rows* is not ``None``) that rows 3..``2+expected_data_rows``
    are present contiguously with no gaps -- both fail loudly on a layout
    change rather than silently misreading a shifted sheet. Pass
    ``expected_data_rows=None`` to skip the row-count check (used by tests
    against a deliberately tiny synthetic workbook).
    """
    with zipfile.ZipFile(xlsx_path) as zf:
        sheet_xml = _find_sheet_xml(zf, _SHEET_NAME)
        shared_strings = _load_shared_strings(zf)

    root = ET.fromstring(sheet_xml)
    sheet_data = root.find("m:sheetData", _NS)
    if sheet_data is None:
        raise ValueError(f"sheet {_SHEET_NAME!r} has no <sheetData>")

    rows_by_num = {}
    for row_el in sheet_data.findall("m:row", _NS):
        r = row_el.get("r")
        if r is not None:
            rows_by_num[int(r)] = row_el

    header_row = rows_by_num.get(_HEADER_ROW_NUM)
    if header_row is None:
        raise ValueError(
            f"sheet {_SHEET_NAME!r} has no row {_HEADER_ROW_NUM} (expected header row)"
        )

    header_cells = _row_cells(header_row)
    for col, expected in _EXPECTED_HEADERS.items():
        is_str, text = _cell_text(header_cells.get(col), shared_strings)
        if not is_str or text != expected:
            raise ValueError(
                f"unexpected header at {col}{_HEADER_ROW_NUM}: expected {expected!r}, "
                f"got {text!r} -- sheet {_SHEET_NAME!r}'s layout may have changed"
            )

    data_row_numbers = sorted(n for n in rows_by_num if n >= _FIRST_DATA_ROW_NUM)
    if expected_data_rows is not None:
        expected_numbers = list(
            range(_FIRST_DATA_ROW_NUM, _FIRST_DATA_ROW_NUM + expected_data_rows)
        )
        if data_row_numbers != expected_numbers:
            raise ValueError(
                f"expected {expected_data_rows} contiguous data row(s) starting at row "
                f"{_FIRST_DATA_ROW_NUM}; found {len(data_row_numbers)} data row(s) "
                f"(first={data_row_numbers[0] if data_row_numbers else None}, "
                f"last={data_row_numbers[-1] if data_row_numbers else None})"
            )

    out = []
    for n in data_row_numbers:
        row_el = rows_by_num[n]
        cells = _row_cells(row_el)

        gene_id_is_str, gene_id = _cell_text(cells.get("B"), shared_strings)
        canon_is_str, canonical = _cell_text(cells.get("C"), shared_strings)
        _, approved_symbol = _cell_text(cells.get("D"), shared_strings)
        gevir_is_str, gevir_raw = _cell_text(cells.get("F"), shared_strings)
        virlof_is_str, virlof_raw = _cell_text(cells.get("H"), shared_strings)

        if not gene_id_is_str or not gene_id:
            raise ValueError(f"row {n}: gene_id (column B) is not a usable string cell")
        if not canon_is_str or not canonical:
            raise ValueError(
                f"row {n}: canonical_transcript (column C) is not a usable string cell"
            )
        if gevir_is_str or gevir_raw is None:
            raise ValueError(f"row {n}: gevir_percentile (column F) is not numeric")
        if virlof_is_str or virlof_raw is None:
            raise ValueError(f"row {n}: virlof_percentile (column H) is not numeric")

        name_is_str, name_raw = _cell_text(cells.get("A"), shared_strings)
        name = _restore_gene_name(
            is_string=name_is_str,
            raw_value=name_raw,
            gene_id=gene_id,
            approved_symbol=approved_symbol,
            cell_ref=f"A{n}",
        )

        out.append((name, gene_id, canonical, gevir_raw, virlof_raw))

    return out


# ---------------------------------------------------------------------------
# Download + write
# ---------------------------------------------------------------------------


def download_gevir_metrics(output_path: str, overwrite: bool = False) -> str:
    """Download the GeVIR supplementary workbook and write the derived TSV.

    Downloads ``GEVIR_SUPPLEMENTARY_URL`` to a temporary file, extracts sheet
    ``table_2`` (see :func:`extract_gevir_table_2`), and writes the result as
    a BGZF-compressed TSV at *output_path*. The final write is atomic: a temp
    file in *output_path*'s own directory, published via ``os.replace`` only
    once the full TSV has been written successfully.

    Parameters
    ----------
    output_path : str
        Destination ``.tsv.bgz`` path.
    overwrite : bool
        Whether to overwrite an existing file at *output_path*.

    Returns
    -------
    str
        *output_path*, unchanged, for chaining.

    Raises
    ------
    FileExistsError
        If *output_path* exists and *overwrite* is False.
    ValueError
        If the workbook's layout does not match what this module expects
        (see :func:`extract_gevir_table_2`).
    """
    output = Path(output_path)
    if output.exists() and not overwrite:
        raise FileExistsError(f"File already exists: {output_path}")
    output.parent.mkdir(parents=True, exist_ok=True)

    xlsx_filename = GEVIR_SUPPLEMENTARY_URL.rsplit("/", 1)[-1]
    with tempfile.TemporaryDirectory(dir=str(output.parent)) as tmp_dir:
        logger.info(
            "Downloading GeVIR supplementary workbook from %s", GEVIR_SUPPLEMENTARY_URL
        )
        xlsx_path = download_file(GEVIR_SUPPLEMENTARY_URL, tmp_dir, xlsx_filename)

        logger.info("Extracting sheet 'table_2' and restoring gene names...")
        # _EXPECTED_DATA_ROWS looked up by name (not as a bound default) so tests can
        # monkeypatch it against a deliberately tiny synthetic workbook.
        rows = extract_gevir_table_2(xlsx_path, expected_data_rows=_EXPECTED_DATA_ROWS)
        logger.info("Extracted %d data rows", len(rows))

        fd, tmp_out = tempfile.mkstemp(
            dir=str(output.parent), prefix=f".{output.name}."
        )
        os.close(fd)
        try:
            with BgzfWriter(tmp_out) as writer:
                writer.write("\t".join(_OUTPUT_HEADER) + "\n")
                for row in rows:
                    writer.write("\t".join(row) + "\n")
            os.replace(tmp_out, output_path)
        finally:
            if os.path.exists(tmp_out):
                os.remove(tmp_out)

    logger.info("Wrote %s", output_path)
    return str(output)


def download_dataset(raw_dir: str, overwrite: bool = False, **kwargs) -> str:
    """Lifecycle entry point for the plugin loader.

    Per the ``DatasetSpec`` contract (``hvantk.core.plugin.api``), a lifecycle
    ``download_fn`` accepts ``raw_dir=<path>`` and writes the raw upstream
    file(s) under that directory. gevir declares no ``lifecycle.parse``, so
    the single file written here (``GEVIR_RAW_FILENAME``) is handed straight
    to the builder via ``_resolve_gevir_path`` (``hvantk/skills/gevir/
    builder.py``), mirroring ``hgnc``.

    Parameters
    ----------
    raw_dir : str
        Directory under which the derived TSV is placed.
    overwrite : bool, optional
        Whether to overwrite an existing file (default: False).
    **kwargs
        Reserved for future lifecycle keyword arguments; ignored today.

    Returns
    -------
    str
        Path to the written raw file.
    """
    raw = Path(raw_dir)
    raw.mkdir(parents=True, exist_ok=True)
    target = raw / GEVIR_RAW_FILENAME
    return download_gevir_metrics(str(target), overwrite=overwrite)


@click.command("gevir-download", short_help="Download GeVIR gene metrics")
@click.option(
    "--output",
    "output_path",
    type=str,
    required=True,
    help="Output path for the derived, BGZF-compressed TSV.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing file if present.",
)
@click.pass_context
def download_cmd(ctx, output_path, overwrite):
    """Download the GeVIR Supplementary Table 2 workbook and derive the metrics TSV.

    Fetches the paper's Supplementary Table 2 (an Excel workbook) and extracts
    sheet ``table_2`` into the bgzipped TSV ``build_gevir_metrics`` reads:
    ``gnomad_gene_name``, ``gene_id``, ``canonical_transcript``, ``gevir_pct``,
    ``virlof_pct`` for 19,361 protein-coding genes.

    For more information: see ``hvantk/skills/gevir/SKILL.md``.

    Examples:

        # Download and derive the GeVIR metrics TSV

        hvantk download gevir --output data/gevir/gevir_metrics_pmid31873297.tsv.bgz

        # Overwrite an existing file

        hvantk download gevir --output data/gevir/gevir_metrics_pmid31873297.tsv.bgz --overwrite
    """
    click.echo("Downloading GeVIR supplementary workbook...")
    click.echo(f"Source: {GEVIR_SUPPLEMENTARY_URL}")

    try:
        result_path = download_gevir_metrics(output_path, overwrite=overwrite)
        click.echo(f"Downloaded to: {result_path}")
    except FileExistsError as e:
        click.echo(f"Error: {e}", err=True)
        click.echo("Use --overwrite to replace existing file.")
        ctx.exit(1)
    except (ValueError, RuntimeError) as e:
        click.echo(f"Download failed: {e}", err=True)
        ctx.exit(1)


if __name__ == "__main__":
    download_cmd()
