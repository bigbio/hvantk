"""Unit test for the gevir:metrics downloader (issue #386).

Mocked (no network): builds tiny synthetic ``.xlsx`` workbooks -- sheet
``table_2``, a row-2 header, and a handful of data rows -- in-memory with
``zipfile``, monkeypatches the network fetch to copy a local file instead of
downloading, and runs each through the real ``download_gevir_metrics``.

One parametrized test function (no new test functions, per review): the
"happy path" case builds a single workbook covering every restoration
family ``cli.py`` handles (MARCHF, SEPTIN-family-not-exercised-here/MARCHF
again via a shared serial, MTARC, SELENOF, and the blank-approved_symbol
DEC1 override) plus full-precision passthrough for both a plain-decimal and
a scientific-notation percentile, and asserts the exact output TSV. The
remaining cases each build a one-row workbook with a single defect and
assert ``pytest.raises`` -- these pin the specific mutants a reviewer found
would otherwise stay green: swapping the MTARC template for ``MARCH{}``,
swapping ``SEP15`` for ``SEPT15``, emptying ``_GENE_ID_OVERRIDES``, deleting
the header-text check, deleting the serial/approved_symbol cross-check, and
accepting an error/boolean-typed cell as numeric.

All gene symbols, Ensembl/Ensembl-transcript ids, and Excel date serials
below may be real (they are identifiers, not the paper's data -- some, like
row 805's serial 38777 or row 14875's serial 36951, are the live workbook's
own values for those exact rows, used here so the restoration math is
checked against a real case). Every ``gevir_pct``/``virlof_pct`` value is
fabricated and does not match the real, published GeVIR numbers.
"""

from __future__ import annotations

import gzip
import os
import stat
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import pytest

from hvantk.core.plugin.drift_runner import _current_umask
from hvantk.skills.gevir import cli as gevir_cli

_MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"

# Row-2 header text the real workbook uses for every column the parser reads
# (hvantk/skills/gevir/cli.py::_EXPECTED_HEADERS). Tests override individual
# entries to exercise the header-validation failure path.
_DEFAULT_HEADERS = {
    "A": "gnomad_gene_name",
    "B": "gene_id",
    "C": "canonical_transcript",
    "D": "approved_symbol",
    "F": "gevir_percentile",
    "H": "virlof_percentile",
}


class _WorkbookBuilder:
    """Accumulates a shared-string table while XML row fragments are built."""

    def __init__(self) -> None:
        self._strings: list[str] = []
        self._index: dict[str, int] = {}

    def _sid(self, text: str) -> int:
        if text not in self._index:
            self._index[text] = len(self._strings)
            self._strings.append(text)
        return self._index[text]

    def string_cell(self, ref: str, text: str) -> str:
        return f'<c r="{ref}" t="s"><v>{self._sid(text)}</v></c>'

    def numeric_cell(self, ref: str, raw: str, cell_type: str | None = None) -> str:
        t_attr = f' t="{cell_type}"' if cell_type else ""
        return f'<c r="{ref}"{t_attr}><v>{raw}</v></c>'

    def sst_xml(self) -> bytes:
        items = "".join(f"<si><t>{s}</t></si>" for s in self._strings)
        xml = (
            f'<sst xmlns="{_MAIN_NS}" count="{len(self._strings)}" '
            f'uniqueCount="{len(self._strings)}">{items}</sst>'
        )
        return xml.encode("utf-8")


def _header_row_xml(builder: _WorkbookBuilder, headers: dict) -> str:
    cells = "".join(
        builder.string_cell(f"{col}2", text) for col, text in headers.items()
    )
    return f'<row r="2">{cells}</row>'


def _data_row_xml(builder: _WorkbookBuilder, row_num: int, row: dict) -> str:
    cells = []
    if row.get("name_kind", "string") == "serial":
        cells.append(builder.numeric_cell(f"A{row_num}", row["name"]))
    else:
        cells.append(builder.string_cell(f"A{row_num}", row["name"]))
    cells.append(builder.string_cell(f"B{row_num}", row["gene_id"]))
    cells.append(builder.string_cell(f"C{row_num}", row["canonical"]))
    if row.get("approved_symbol") is not None:
        cells.append(builder.string_cell(f"D{row_num}", row["approved_symbol"]))
    cells.append(
        builder.numeric_cell(f"F{row_num}", row["gevir"], row.get("gevir_cell_type"))
    )
    cells.append(
        builder.numeric_cell(f"H{row_num}", row["virlof"], row.get("virlof_cell_type"))
    )
    return f'<row r="{row_num}">' + "".join(cells) + "</row>"


def _workbook_xml() -> bytes:
    xml = (
        f'<workbook xmlns="{_MAIN_NS}" xmlns:r="{_REL_NS}">'
        '<sheets><sheet name="table_2" sheetId="1" r:id="rId1"/></sheets>'
        "</workbook>"
    )
    return xml.encode("utf-8")


def _workbook_rels_xml() -> bytes:
    xml = (
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        'Target="worksheets/sheet1.xml"/>'
        "</Relationships>"
    )
    return xml.encode("utf-8")


def _build_synthetic_workbook(
    path: Path, rows: list[dict], *, header_overrides: dict | None = None
) -> None:
    """Write a minimal, valid ``.xlsx`` at *path* with one row per *rows* entry.

    Only the parts ``hvantk/skills/gevir/cli.py`` actually reads are present
    (no ``[Content_Types].xml``, ``styles.xml`` or theme). Every XML fragment
    round-trips through ``ET.fromstring`` before being written, so a malformed
    fixture fails at build time rather than inside the code under test.
    """
    headers = dict(_DEFAULT_HEADERS)
    if header_overrides:
        headers.update(header_overrides)

    builder = _WorkbookBuilder()
    row_xml = [_header_row_xml(builder, headers)]
    row_xml.extend(
        _data_row_xml(builder, n, row) for n, row in enumerate(rows, start=3)
    )
    sheet_xml = (
        f'<worksheet xmlns="{_MAIN_NS}"><sheetData>{"".join(row_xml)}'
        "</sheetData></worksheet>"
    ).encode("utf-8")

    parts = {
        "xl/workbook.xml": _workbook_xml(),
        "xl/_rels/workbook.xml.rels": _workbook_rels_xml(),
        "xl/sharedStrings.xml": builder.sst_xml(),
        "xl/worksheets/sheet1.xml": sheet_xml,
    }
    for data in parts.values():
        ET.fromstring(data)  # fail fast on a malformed fixture, not deep in the parser

    with zipfile.ZipFile(path, "w") as zf:
        for name, data in parts.items():
            zf.writestr(name, data)


# ---------------------------------------------------------------------------
# Row fixtures. gene_id/canonical/serial values are real (identifiers only);
# every gevir/virlof value is fabricated.
# ---------------------------------------------------------------------------

_ORDINARY_ROW = {
    "name": "TBX15",
    "gene_id": "ENSG00000092607",
    "canonical": "ENST00000207157",
    "gevir": "12.3456789012",
    "virlof": "87.6543210987",
}
_ORDINARY_LINE = "TBX15\tENSG00000092607\tENST00000207157\t12.3456789012\t87.6543210987"

# Serial 38777 is the real live workbook's value for this exact row
# (2006-03-01 on the 1899-12-30 epoch); approved_symbol MARCHF6 must restore
# it to MARCH6 (the MARCHF-family pattern).
_MARCHF6_ROW = {
    "name": "38777",
    "name_kind": "serial",
    "gene_id": "ENSG00000145495",
    "canonical": "ENST00000274140",
    "approved_symbol": "MARCHF6",
    "gevir": "1.1111111111",
    "virlof": "2.2222222222",
}
_MARCHF6_LINE = "MARCH6\tENSG00000145495\tENST00000274140\t1.1111111111\t2.2222222222"

# Serial 36951 is the real live workbook's value shared by TWO distinct rows
# (2001-03-01): one approved_symbol MTARC1 (-> MARC1, not MARCH1 -- MARC and
# MARCH parse to the same date and can only be told apart via approved_symbol),
# one approved_symbol MARCHF1 (-> MARCH1). Both in the same workbook pins that
# the serial alone never decides the family.
_MTARC1_ROW = {
    "name": "36951",
    "name_kind": "serial",
    "gene_id": "ENSG00000186205",
    "canonical": "ENST00000366910",
    "approved_symbol": "MTARC1",
    "gevir": "3.3333333333",
    "virlof": "4.4444444444",
}
_MTARC1_LINE = "MARC1\tENSG00000186205\tENST00000366910\t3.3333333333\t4.4444444444"

_MARCHF1_ROW = {
    "name": "36951",
    "name_kind": "serial",
    "gene_id": "ENSG00000145416",
    "canonical": "ENST00000503008",
    "approved_symbol": "MARCHF1",
    "gevir": "5.5555555555",
    "virlof": "6.6666666666",
}
_MARCHF1_LINE = "MARCH1\tENSG00000145416\tENST00000503008\t5.5555555555\t6.6666666666"

# Serial 42248 is the real live workbook's value for SELENOF (2015-09-01),
# which restores to SEP15 -- not SEPT15. gevir uses scientific notation (a
# fabricated value, not the real one) so the full-precision passthrough test
# covers both decimal shapes the real workbook uses (see this plugin's
# fixture README).
_SELENOF_ROW = {
    "name": "42248",
    "name_kind": "serial",
    "gene_id": "ENSG00000183291",
    "canonical": "ENST00000331835",
    "approved_symbol": "SELENOF",
    "gevir": "4.13E-2",
    "virlof": "7.7777777777",
}
_SELENOF_LINE = "SEP15\tENSG00000183291\tENST00000331835\t4.13E-2\t7.7777777777"

# Serial 37226 is the real live workbook's value for the gene now named
# DELEC1, whose approved_symbol cell is blank (column D omitted entirely) --
# it can only be resolved via the gene-id-keyed override, restoring DEC1.
_DEC1_ROW = {
    "name": "37226",
    "name_kind": "serial",
    "gene_id": "ENSG00000173077",
    "canonical": "ENST00000374016",
    "approved_symbol": None,
    "gevir": "8.8888888888",
    "virlof": "9.9999999999",
}
_DEC1_LINE = "DEC1\tENSG00000173077\tENST00000374016\t8.8888888888\t9.9999999999"

_HAPPY_PATH_ROWS = [
    _ORDINARY_ROW,
    _MARCHF6_ROW,
    _MTARC1_ROW,
    _MARCHF1_ROW,
    _SELENOF_ROW,
    _DEC1_ROW,
]
_HAPPY_PATH_LINES = [
    _ORDINARY_LINE,
    _MARCHF6_LINE,
    _MTARC1_LINE,
    _MARCHF1_LINE,
    _SELENOF_LINE,
    _DEC1_LINE,
]

# Serial 39142 is the real live workbook's value for MARCHF7 (2007-03-01), not
# MARCHF6 (number 6). Pairing it with approved_symbol MARCHF6 is a deliberate
# mismatch: the cross-check between the decoded serial and the number
# approved_symbol implies must raise rather than silently emit MARCH6 anyway.
_MISMATCHED_SERIAL_ROW = {
    "name": "39142",
    "name_kind": "serial",
    "gene_id": "ENSG00000136536",
    "canonical": "ENST00000259050",
    "approved_symbol": "MARCHF6",
    "gevir": "1.0",
    "virlof": "2.0",
}

_ERROR_TYPED_CELL_ROW = {
    "name": "TBX15",
    "gene_id": "ENSG00000092607",
    "canonical": "ENST00000207157",
    "gevir": "1.0",
    "virlof": "2.0",
    "gevir_cell_type": "e",  # OOXML error cell (<v> would hold "#N/A" et al.)
}

_BOOLEAN_TYPED_CELL_ROW = {
    "name": "TBX15",
    "gene_id": "ENSG00000092607",
    "canonical": "ENST00000207157",
    "gevir": "1.0",
    "virlof": "1",
    "virlof_cell_type": "b",  # OOXML boolean cell (<v> is "0"/"1")
}


def _run_download(tmp_path, monkeypatch, rows, *, header_overrides=None):
    """Build a synthetic workbook from *rows*, mock the fetch, and download it.

    Returns the output path. Raises whatever ``download_gevir_metrics`` raises
    (callers wrap this in ``pytest.raises`` for the failure cases).
    """
    workbook_path = tmp_path / "source.xlsx"
    _build_synthetic_workbook(workbook_path, rows, header_overrides=header_overrides)

    def fake_download_file(url, out_dir, file_name):
        dest = Path(out_dir) / file_name
        dest.write_bytes(workbook_path.read_bytes())
        return str(dest)

    monkeypatch.setattr(gevir_cli, "download_file", fake_download_file)
    monkeypatch.setattr(gevir_cli, "_EXPECTED_DATA_ROWS", len(rows))

    output_path = tmp_path / "gevir_metrics_pmid31873297.tsv.bgz"
    gevir_cli.download_gevir_metrics(str(output_path))
    return output_path


@pytest.mark.parametrize(
    "rows, header_overrides, expect_error, expected_lines",
    [
        pytest.param(
            _HAPPY_PATH_ROWS,
            None,
            None,
            _HAPPY_PATH_LINES,
            id="happy-path-all-families-and-precision-shapes",
        ),
        pytest.param(
            [_MISMATCHED_SERIAL_ROW],
            None,
            "refusing",
            None,
            id="mismatched-serial-vs-approved-symbol-raises",
        ),
        pytest.param(
            [_ORDINARY_ROW],
            {"H": "virlof_percentile_renamed"},
            "unexpected header",
            None,
            id="renamed-header-raises",
        ),
        pytest.param(
            [_ERROR_TYPED_CELL_ROW],
            None,
            "unsupported type",
            None,
            id="error-typed-numeric-cell-raises",
        ),
        pytest.param(
            [_BOOLEAN_TYPED_CELL_ROW],
            None,
            "unsupported type",
            None,
            id="boolean-typed-numeric-cell-raises",
        ),
    ],
)
def test_download_gevir_metrics_restores_name_and_keeps_full_precision(
    tmp_path, monkeypatch, rows, header_overrides, expect_error, expected_lines
):
    if expect_error is not None:
        with pytest.raises(ValueError, match=expect_error):
            _run_download(
                tmp_path, monkeypatch, rows, header_overrides=header_overrides
            )
        return

    output_path = _run_download(tmp_path, monkeypatch, rows)

    assert output_path.exists()
    # No leftover temp file from the atomic-publish step.
    assert [p for p in tmp_path.glob(".*") if p.is_file()] == []

    # mkstemp always creates mode 0600; the published file must not keep that
    # -- it should land at whatever a normal file create gets (0666 masked by
    # the process umask), so a shared raw dir stays group-readable.
    mode = stat.S_IMODE(os.stat(output_path).st_mode)
    assert mode == (0o666 & ~_current_umask())

    text = gzip.decompress(output_path.read_bytes()).decode("utf-8")
    lines = text.splitlines()

    assert lines[0] == (
        "gnomad_gene_name\tgene_id\tcanonical_transcript\tgevir_pct\tvirlof_pct"
    )
    assert lines[1:] == expected_lines
