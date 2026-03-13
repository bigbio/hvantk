"""Tests for the hvantk prepare-geneset CLI command."""

import json
from pathlib import Path
from unittest.mock import patch

import pytest
from click.testing import CliRunner

from hvantk.commands.prepare_geneset_cli import prepare_geneset_cmd

TESTDATA = Path(__file__).parent / "testdata" / "prepare_geneset"


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def output_json(tmp_path):
    return str(tmp_path / "output.json")


# ---------------------------------------------------------------------------
# Basic conversion
# ---------------------------------------------------------------------------


class TestBasicConversion:
    def test_valid_tsv_produces_json(self, runner, output_json):
        result = runner.invoke(
            prepare_geneset_cmd,
            ["-i", str(TESTDATA / "valid_panels.tsv"), "-o", output_json],
        )
        assert result.exit_code == 0, result.output
        data = json.loads(Path(output_json).read_text())
        assert "gene_sets" in data
        assert "cardiac_panel" in data["gene_sets"]
        assert "epilepsy_panel" in data["gene_sets"]
        assert "ras_pathway" in data["gene_sets"]
        assert data["n_gene_sets"] == 3
        assert "background_genes" in data
        assert data["metadata"]["created_by"] == "hvantk prepare-geneset"
        assert data["metadata"]["hgnc_validated"] is False

    def test_gene_set_contents(self, runner, output_json):
        runner.invoke(
            prepare_geneset_cmd,
            ["-i", str(TESTDATA / "valid_panels.tsv"), "-o", output_json],
        )
        data = json.loads(Path(output_json).read_text())
        cardiac = data["gene_sets"]["cardiac_panel"]
        assert set(cardiac["genes"]) == {"MYH7", "TNNT2", "LMNA", "SCN5A", "TTN"}
        assert cardiac["n_genes"] == 5
        assert cardiac["source"] == "prepare-geneset"

    def test_background_is_union(self, runner, output_json):
        runner.invoke(
            prepare_geneset_cmd,
            ["-i", str(TESTDATA / "valid_panels.tsv"), "-o", output_json],
        )
        data = json.loads(Path(output_json).read_text())
        bg = set(data["background_genes"])
        # Background should be union of all genes.
        all_genes = set()
        for gs in data["gene_sets"].values():
            all_genes.update(gs["genes"])
        assert bg == all_genes


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------


class TestErrorHandling:
    def test_missing_input(self, runner, output_json):
        result = runner.invoke(
            prepare_geneset_cmd,
            ["-i", "/nonexistent/file.tsv", "-o", output_json],
        )
        assert result.exit_code != 0

    def test_wrong_columns(self, runner, output_json):
        result = runner.invoke(
            prepare_geneset_cmd,
            ["-i", str(TESTDATA / "wrong_columns.tsv"), "-o", output_json],
        )
        assert result.exit_code != 0
        assert "expected 2 tab-separated columns" in result.output

    def test_ensembl_ids(self, runner, output_json):
        result = runner.invoke(
            prepare_geneset_cmd,
            ["-i", str(TESTDATA / "with_ensembl_ids.tsv"), "-o", output_json],
        )
        assert result.exit_code != 0
        assert "Ensembl" in result.output

    def test_overwrite_refused(self, runner, tmp_path):
        out = str(tmp_path / "out.json")
        # Create the file first.
        runner.invoke(
            prepare_geneset_cmd,
            ["-i", str(TESTDATA / "valid_panels.tsv"), "-o", out],
        )
        # Attempt without --overwrite.
        result = runner.invoke(
            prepare_geneset_cmd,
            ["-i", str(TESTDATA / "valid_panels.tsv"), "-o", out],
        )
        assert result.exit_code != 0
        assert "already exists" in result.output

    def test_overwrite_allowed(self, runner, tmp_path):
        out = str(tmp_path / "out.json")
        runner.invoke(
            prepare_geneset_cmd,
            ["-i", str(TESTDATA / "valid_panels.tsv"), "-o", out],
        )
        result = runner.invoke(
            prepare_geneset_cmd,
            [
                "-i", str(TESTDATA / "valid_panels.tsv"),
                "-o", out,
                "--overwrite",
            ],
        )
        assert result.exit_code == 0, result.output


# ---------------------------------------------------------------------------
# HGNC validation
# ---------------------------------------------------------------------------


def _mock_hgnc_maps():
    canonical = {
        "BRCA1", "BRCA2", "TP53", "EGFR", "ERCC1",
        "MYH7", "TNNT2", "LMNA", "SCN5A", "TTN",
        "SCN1A", "SCN2A", "KCNQ2", "STXBP1", "GABRA1",
        "BRAF", "KRAS", "NRAS", "HRAS", "MAP2K1",
    }
    alias_to_canonical = {
        "FANCD1": "BRCA2",
        "ERCC11": "ERCC1",
    }
    canonical_to_aliases = {
        "BRCA2": ["FANCD1"],
        "ERCC1": ["ERCC11"],
    }
    return canonical, alias_to_canonical, canonical_to_aliases


class TestHgncValidation:
    @patch("hvantk.utils.gene_aliases._load_hgnc_symbol_maps")
    def test_with_hgnc(self, mock_load, runner, output_json):
        mock_load.return_value = _mock_hgnc_maps()
        result = runner.invoke(
            prepare_geneset_cmd,
            [
                "-i", str(TESTDATA / "valid_panels.tsv"),
                "-o", output_json,
                "--hgnc", str(TESTDATA / "valid_panels.tsv"),  # path just needs to exist
            ],
        )
        assert result.exit_code == 0, result.output
        data = json.loads(Path(output_json).read_text())
        assert data["metadata"]["hgnc_validated"] is True

    @patch("hvantk.utils.gene_aliases._load_hgnc_symbol_maps")
    def test_aliases_resolved_in_output(self, mock_load, runner, output_json):
        mock_load.return_value = _mock_hgnc_maps()
        result = runner.invoke(
            prepare_geneset_cmd,
            [
                "-i", str(TESTDATA / "with_aliases.tsv"),
                "-o", output_json,
                "--hgnc", str(TESTDATA / "valid_panels.tsv"),
            ],
        )
        assert result.exit_code == 0, result.output
        data = json.loads(Path(output_json).read_text())
        # FANCD1 should be resolved to BRCA2.
        panel_a_genes = set(data["gene_sets"]["panel_a"]["genes"])
        assert "BRCA2" in panel_a_genes
        assert "FANCD1" not in panel_a_genes
        # ERCC11 should be resolved to ERCC1.
        panel_b_genes = set(data["gene_sets"]["panel_b"]["genes"])
        assert "ERCC1" in panel_b_genes
        assert "ERCC11" not in panel_b_genes
        # Metadata should contain alias mappings.
        assert "aliases_resolved" in data["metadata"]


# ---------------------------------------------------------------------------
# Filtering and background
# ---------------------------------------------------------------------------


class TestFilteringAndBackground:
    def test_min_genes_filters(self, runner, output_json):
        result = runner.invoke(
            prepare_geneset_cmd,
            [
                "-i", str(TESTDATA / "with_duplicates.tsv"),
                "-o", output_json,
                "--min-genes", "3",
            ],
        )
        assert result.exit_code == 0, result.output
        data = json.loads(Path(output_json).read_text())
        # Both panels have 2 genes (after dedup), so both should be filtered.
        assert data["n_gene_sets"] == 0

    def test_background_file(self, runner, output_json):
        result = runner.invoke(
            prepare_geneset_cmd,
            [
                "-i", str(TESTDATA / "valid_panels.tsv"),
                "-o", output_json,
                "--background", str(TESTDATA / "background_genes.txt"),
            ],
        )
        assert result.exit_code == 0, result.output
        data = json.loads(Path(output_json).read_text())
        bg = set(data["background_genes"])
        # Background should include genes from the file.
        assert "PTEN" in bg
        assert "APC" in bg


# ---------------------------------------------------------------------------
# GMT export
# ---------------------------------------------------------------------------


class TestGmtExport:
    def test_export_gmt(self, runner, tmp_path):
        out_json = str(tmp_path / "out.json")
        out_gmt = str(tmp_path / "out.gmt")
        result = runner.invoke(
            prepare_geneset_cmd,
            [
                "-i", str(TESTDATA / "valid_panels.tsv"),
                "-o", out_json,
                "--export-gmt", out_gmt,
            ],
        )
        assert result.exit_code == 0, result.output
        assert Path(out_gmt).exists()
        lines = Path(out_gmt).read_text().strip().split("\n")
        assert len(lines) == 3  # 3 gene sets
        # GMT format: name<TAB>description<TAB>gene1<TAB>gene2...
        for line in lines:
            fields = line.split("\t")
            assert len(fields) >= 3
