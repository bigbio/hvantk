"""Tests for the hvantk genesets prepare CLI command."""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from hvantk.tools.genesets.genesets_cli import genesets_prepare as prepare_geneset_cmd

TESTDATA = Path(__file__).parent / "testdata" / "prepare_geneset"


def _make_mock_catalog_cli():
    """Return a mock satisfying GeneCatalogStreamer for CLI tests.

    validate_symbols is wired to the canonical / alias maps so that
    validate_with_catalog produces the expected validation result.
    """
    from hvantk.core.streamers.gene_catalog import GeneCatalogStreamer

    canonical = {
        "BRCA1",
        "BRCA2",
        "TP53",
        "EGFR",
        "ERCC1",
        "MYH7",
        "TNNT2",
        "LMNA",
        "SCN5A",
        "TTN",
        "SCN1A",
        "SCN2A",
        "KCNQ2",
        "STXBP1",
        "GABRA1",
        "BRAF",
        "KRAS",
        "NRAS",
        "HRAS",
        "MAP2K1",
    }
    alias_to_canonical = {"FANCD1": "BRCA2", "ERCC11": "ERCC1"}

    mock_catalog = MagicMock(spec=GeneCatalogStreamer)
    # Wire validate_symbols to the concrete ABC implementation logic.
    mock_catalog.is_canonical.side_effect = lambda s: s in canonical
    mock_catalog.resolve_alias.side_effect = lambda s: alias_to_canonical.get(s)

    def _validate_symbols(symbols):
        return GeneCatalogStreamer.validate_symbols(mock_catalog, symbols)

    mock_catalog.validate_symbols.side_effect = _validate_symbols
    return mock_catalog


@pytest.fixture
def runner():
    return CliRunner()


def test_valid_tsv_produces_complete_json(runner, tmp_path):
    """End-to-end: valid TSV → JSON with correct structure, contents, and background."""
    out = str(tmp_path / "output.json")
    result = runner.invoke(
        prepare_geneset_cmd,
        ["-i", str(TESTDATA / "valid_panels.tsv"), "-o", out],
    )
    assert result.exit_code == 0, result.output
    data = json.loads(Path(out).read_text())

    # Structure
    assert data["n_gene_sets"] == 3
    assert data["metadata"]["hgnc_validated"] is False
    assert set(data["gene_sets"].keys()) == {
        "cardiac_panel",
        "epilepsy_panel",
        "ras_pathway",
    }

    # Gene set contents
    cardiac = data["gene_sets"]["cardiac_panel"]
    assert set(cardiac["genes"]) == {"MYH7", "TNNT2", "LMNA", "SCN5A", "TTN"}
    assert cardiac["n_genes"] == 5

    # Background is union of all genes
    bg = set(data["background_genes"])
    all_genes = set()
    for gs in data["gene_sets"].values():
        all_genes.update(gs["genes"])
    assert bg == all_genes


def test_cli_rejects_invalid_input(runner, tmp_path):
    # Missing file
    result = runner.invoke(
        prepare_geneset_cmd,
        ["-i", "/nonexistent/file.tsv", "-o", str(tmp_path / "out.json")],
    )
    assert result.exit_code != 0

    # Wrong column count
    bad_tsv = tmp_path / "bad.tsv"
    bad_tsv.write_text("panel_a\tBRCA1\npanel_a\tTP53\textra\n")
    result = runner.invoke(
        prepare_geneset_cmd,
        ["-i", str(bad_tsv), "-o", str(tmp_path / "out2.json")],
    )
    assert result.exit_code != 0
    assert "expected 2 tab-separated columns" in result.output


def test_overwrite_flag(runner, tmp_path):
    """Refuse overwrite without flag; allow with --overwrite."""
    out = str(tmp_path / "out.json")
    runner.invoke(
        prepare_geneset_cmd, ["-i", str(TESTDATA / "valid_panels.tsv"), "-o", out]
    )
    # Without --overwrite → fail
    result = runner.invoke(
        prepare_geneset_cmd, ["-i", str(TESTDATA / "valid_panels.tsv"), "-o", out]
    )
    assert result.exit_code != 0
    assert "already exists" in result.output
    # With --overwrite → succeed
    result = runner.invoke(
        prepare_geneset_cmd,
        ["-i", str(TESTDATA / "valid_panels.tsv"), "-o", out, "--overwrite"],
    )
    assert result.exit_code == 0


@patch(
    "hvantk.skills.hgnc.streamers.HGNCGeneCatalogStreamer.from_path",
    side_effect=lambda path: _make_mock_catalog_cli(),
)
def test_hgnc_validation_resolves_aliases(mock_from_path, runner, tmp_path):
    """HGNC validation flag resolves aliases in output."""
    out = str(tmp_path / "out.json")
    result = runner.invoke(
        prepare_geneset_cmd,
        [
            "-i",
            str(TESTDATA / "with_aliases.tsv"),
            "-o",
            out,
            "--hgnc",
            str(TESTDATA / "valid_panels.tsv"),
        ],
    )
    assert result.exit_code == 0, result.output
    data = json.loads(Path(out).read_text())
    assert data["metadata"]["hgnc_validated"] is True
    assert "BRCA2" in set(data["gene_sets"]["panel_a"]["genes"])
    assert "FANCD1" not in set(data["gene_sets"]["panel_a"]["genes"])


def test_export_gmt(runner, tmp_path):
    """GMT export produces valid GMT format."""
    out_json = str(tmp_path / "out.json")
    out_gmt = str(tmp_path / "out.gmt")
    result = runner.invoke(
        prepare_geneset_cmd,
        [
            "-i",
            str(TESTDATA / "valid_panels.tsv"),
            "-o",
            out_json,
            "--export-gmt",
            out_gmt,
        ],
    )
    assert result.exit_code == 0, result.output
    lines = Path(out_gmt).read_text().strip().split("\n")
    assert len(lines) == 3
    for line in lines:
        assert len(line.split("\t")) >= 3
