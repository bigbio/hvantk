"""Unit test for the gevir:metrics downloader (issue #386).

Mocked (no network): builds a tiny synthetic ``.xlsx`` workbook -- sheet
``table_2``, a valid row-2 header, and two data rows -- in-memory with
``zipfile``, monkeypatches the network fetch to copy that local file instead
of downloading, and runs it through the real ``download_gevir_metrics``.

Covers the two behaviors issue #386 asks for: full-precision passthrough for
an ordinary string-named gene (no rounding of ``gevir_pct``/``virlof_pct``),
and restoration of one Excel-date-mangled ``gnomad_gene_name`` cell (a
MARCHF-family row, decoded via the serial + cross-checked against the row's
own ``approved_symbol``) -- see ``hvantk/skills/gevir/cli.py``'s module
docstring for the restoration rule this exercises.
"""

from __future__ import annotations

import gzip
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

from hvantk.skills.gevir import cli as gevir_cli

_MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"

# Shared-string table, indexed by position (shared-string id = list index).
_SHARED_STRINGS = [
    "gnomad_gene_name",  # 0
    "gene_id",  # 1
    "canonical_transcript",  # 2
    "approved_symbol",  # 3
    "gevir_percentile",  # 4
    "virlof_percentile",  # 5
    "TBX15",  # 6: ordinary row's name
    "ENSG00000092607",  # 7: ordinary row's gene_id
    "ENST00000207157",  # 8: ordinary row's canonical_transcript
    "ENSG00000145495",  # 9: mangled row's gene_id
    "ENST00000274140",  # 10: mangled row's canonical_transcript
    "MARCHF6",  # 11: mangled row's approved_symbol
]


def _sst_xml() -> bytes:
    items = "".join(f"<si><t>{s}</t></si>" for s in _SHARED_STRINGS)
    xml = (
        f'<sst xmlns="{_MAIN_NS}" count="{len(_SHARED_STRINGS)}" '
        f'uniqueCount="{len(_SHARED_STRINGS)}">{items}</sst>'
    )
    return xml.encode("utf-8")


def _s(ref: str, idx: int) -> str:
    """A shared-string cell."""
    return f'<c r="{ref}" t="s"><v>{idx}</v></c>'


def _n(ref: str, value) -> str:
    """A numeric cell (no ``t`` attribute -- Excel's default, and what a
    date-serial cell looks like)."""
    return f'<c r="{ref}"><v>{value}</v></c>'


def _sheet1_xml() -> bytes:
    # Row 2: header. Columns E/G are never read by the parser, so they are
    # omitted entirely -- _row_cells only looks up what it needs.
    row2 = (
        '<row r="2">'
        + _s("A2", 0)
        + _s("B2", 1)
        + _s("C2", 2)
        + _s("D2", 3)
        + _s("F2", 4)
        + _s("H2", 5)
        + "</row>"
    )
    # Row 3: an ordinary gene. F/H carry long, non-round decimals so a silent
    # rounding/reformatting bug would show up immediately.
    row3 = (
        '<row r="3">'
        + _s("A3", 6)
        + _s("B3", 7)
        + _s("C3", 8)
        + _n("F3", "6.7455193429999998")
        + _n("H3", "11.1822736429")
        + "</row>"
    )
    # Row 4: gnomad_gene_name is a numeric Excel-date serial (38777 = the real
    # live workbook's value for this exact row: 2006-03-01 on the 1899-12-30
    # epoch), with approved_symbol MARCHF6 -- the MARCHF family pattern must
    # restore this to "MARCH6".
    row4 = (
        '<row r="4">'
        + _n("A4", "38777")
        + _s("B4", 9)
        + _s("C4", 10)
        + _s("D4", 11)
        + _n("F4", "4.1630081090999997")
        + _n("H4", "1.5804968751999999")
        + "</row>"
    )
    xml = (
        f'<worksheet xmlns="{_MAIN_NS}"><sheetData>{row2}{row3}{row4}'
        "</sheetData></worksheet>"
    )
    return xml.encode("utf-8")


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


def _build_synthetic_workbook(path: Path) -> None:
    """Write a minimal, valid 2-data-row ``.xlsx`` at *path*.

    Only the parts ``hvantk/skills/gevir/cli.py`` actually reads are present
    (no ``[Content_Types].xml``, ``styles.xml`` or theme -- the parser never
    opens them). Every XML fragment round-trips through ``ET.fromstring``
    before being written, so a malformed fixture fails at build time rather
    than inside the code under test.
    """
    parts = {
        "xl/workbook.xml": _workbook_xml(),
        "xl/_rels/workbook.xml.rels": _workbook_rels_xml(),
        "xl/sharedStrings.xml": _sst_xml(),
        "xl/worksheets/sheet1.xml": _sheet1_xml(),
    }
    for data in parts.values():
        ET.fromstring(data)  # fail fast on a malformed fixture, not deep in the parser

    with zipfile.ZipFile(path, "w") as zf:
        for name, data in parts.items():
            zf.writestr(name, data)


def test_download_gevir_metrics_restores_name_and_keeps_full_precision(
    tmp_path, monkeypatch
):
    workbook_path = tmp_path / "source" / "workbook.xlsx"
    workbook_path.parent.mkdir()
    _build_synthetic_workbook(workbook_path)

    # No network: download_file's real signature is (url, out_dir, file_name);
    # stand in with a copy of the local synthetic workbook instead of a fetch.
    def fake_download_file(url, out_dir, file_name):
        dest = Path(out_dir) / file_name
        dest.write_bytes(workbook_path.read_bytes())
        return str(dest)

    monkeypatch.setattr(gevir_cli, "download_file", fake_download_file)
    # The real default (19,361) does not fit this 2-row fixture; this module
    # global is looked up by name inside download_gevir_metrics precisely so
    # a test can override it like this.
    monkeypatch.setattr(gevir_cli, "_EXPECTED_DATA_ROWS", 2)

    output_path = tmp_path / "gevir_metrics_pmid31873297.tsv.bgz"
    result = gevir_cli.download_gevir_metrics(str(output_path))

    assert result == str(output_path)
    assert output_path.exists()
    # No leftover temp file from the atomic-publish step.
    assert list(tmp_path.glob(".*")) == []

    text = gzip.decompress(output_path.read_bytes()).decode("utf-8")
    lines = text.splitlines()

    assert lines[0] == (
        "gnomad_gene_name\tgene_id\tcanonical_transcript\tgevir_pct\tvirlof_pct"
    )
    # Full precision, unrounded -- passed through exactly as the workbook stored it.
    assert lines[1] == (
        "TBX15\tENSG00000092607\tENST00000207157\t6.7455193429999998\t11.1822736429"
    )
    # Restored from the date serial via the MARCHF-family pattern on approved_symbol.
    assert lines[2] == (
        "MARCH6\tENSG00000145495\tENST00000274140\t4.1630081090999997\t1.5804968751999999"
    )
    assert len(lines) == 3
