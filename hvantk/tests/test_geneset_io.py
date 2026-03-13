"""Tests for hvantk.utils.geneset_io — gene set TSV parsing and validation."""

from pathlib import Path
from unittest.mock import patch

import pytest

from hvantk.utils.geneset_io import (
    ParseResult,
    ValidationResult,
    detect_id_type,
    parse_geneset_tsv,
    validate_gene_ids,
    validate_with_hgnc,
)

TESTDATA = Path(__file__).parent / "testdata" / "prepare_geneset"


# ---------------------------------------------------------------------------
# detect_id_type
# ---------------------------------------------------------------------------


class TestDetectIdType:
    def test_ensembl_gene(self):
        assert detect_id_type("ENSG00000141510") == "ensembl_gene"

    def test_ensembl_gene_with_version(self):
        assert detect_id_type("ENSG00000141510.12") == "ensembl_gene"

    def test_ensembl_transcript(self):
        assert detect_id_type("ENST00000269305") == "ensembl_transcript"

    def test_ensembl_transcript_with_version(self):
        assert detect_id_type("ENST00000269305.8") == "ensembl_transcript"

    def test_entrez(self):
        assert detect_id_type("7157") == "entrez"

    def test_mouse_symbol(self):
        assert detect_id_type("Trp53") == "mouse_symbol"
        assert detect_id_type("Scn1a") == "mouse_symbol"
        assert detect_id_type("Brca1") == "mouse_symbol"

    def test_whitespace(self):
        assert detect_id_type("BRCA 1") == "whitespace"
        assert detect_id_type("BRCA\t1") == "whitespace"

    def test_human_symbol(self):
        assert detect_id_type("TP53") == "symbol"
        assert detect_id_type("BRCA1") == "symbol"
        assert detect_id_type("SCN1A") == "symbol"
        assert detect_id_type("HLA-DRB1") == "symbol"
        assert detect_id_type("MYH7") == "symbol"

    def test_short_symbols_not_mouse(self):
        # Short symbols like "Ab" should not be classified as mouse.
        # The regex requires at least 3 chars total.
        assert detect_id_type("TP") == "symbol"


# ---------------------------------------------------------------------------
# validate_gene_ids
# ---------------------------------------------------------------------------


class TestValidateGeneIds:
    def test_all_valid(self):
        valid, problems = validate_gene_ids(["BRCA1", "TP53", "EGFR"])
        assert valid == ["BRCA1", "TP53", "EGFR"]
        assert problems == {}

    def test_mixed(self):
        valid, problems = validate_gene_ids(
            ["BRCA1", "ENSG00000141510", "7157", "TP53"]
        )
        assert valid == ["BRCA1", "TP53"]
        assert "ensembl_gene" in problems
        assert "entrez" in problems

    def test_empty_list(self):
        valid, problems = validate_gene_ids([])
        assert valid == []
        assert problems == {}


# ---------------------------------------------------------------------------
# parse_geneset_tsv
# ---------------------------------------------------------------------------


class TestParseGenesetTsv:
    def test_valid_panels(self):
        result = parse_geneset_tsv(TESTDATA / "valid_panels.tsv")
        assert isinstance(result, ParseResult)
        assert set(result.gene_sets.keys()) == {
            "cardiac_panel",
            "epilepsy_panel",
            "ras_pathway",
        }
        assert result.gene_sets["cardiac_panel"] == [
            "MYH7", "TNNT2", "LMNA", "SCN5A", "TTN",
        ]
        assert result.gene_sets["epilepsy_panel"] == [
            "SCN1A", "SCN2A", "KCNQ2", "STXBP1", "GABRA1",
        ]
        assert result.gene_sets["ras_pathway"] == [
            "BRAF", "KRAS", "NRAS", "HRAS", "MAP2K1",
        ]
        assert result.n_lines_parsed == 15
        assert result.n_duplicates == 0
        assert result.warnings == []

    def test_comments_and_blanks_skipped(self):
        result = parse_geneset_tsv(TESTDATA / "with_comments.tsv")
        assert "cardiac_panel" in result.gene_sets
        assert "epilepsy_panel" in result.gene_sets
        assert result.n_lines_skipped > 0
        assert result.gene_sets["cardiac_panel"] == ["MYH7", "TNNT2", "LMNA"]

    def test_wrong_column_count(self):
        with pytest.raises(ValueError, match="expected 2 tab-separated columns"):
            parse_geneset_tsv(TESTDATA / "wrong_columns.tsv")

    def test_ensembl_ids_rejected(self):
        with pytest.raises(ValueError, match="Ensembl gene IDs"):
            parse_geneset_tsv(TESTDATA / "with_ensembl_ids.tsv")

    def test_mouse_symbols_rejected(self):
        with pytest.raises(ValueError, match="mouse gene symbols"):
            parse_geneset_tsv(TESTDATA / "with_mouse_symbols.tsv")

    def test_entrez_ids_rejected(self):
        with pytest.raises(ValueError, match="Entrez Gene IDs"):
            parse_geneset_tsv(TESTDATA / "with_entrez_ids.tsv")

    def test_duplicates_deduped_with_warning(self):
        result = parse_geneset_tsv(TESTDATA / "with_duplicates.tsv")
        assert result.gene_sets["panel_a"] == ["BRCA1", "TP53"]
        assert result.gene_sets["panel_b"] == ["EGFR", "KRAS"]
        assert result.n_duplicates == 2
        assert any("Duplicate gene" in w for w in result.warnings)

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            parse_geneset_tsv(Path("/nonexistent/file.tsv"))

    def test_mixed_valid_and_few_mouse_like_passes(self):
        """A few mouse-like symbols among mostly valid should not error."""
        import tempfile

        content = (
            "panel\tBRCA1\n"
            "panel\tTP53\n"
            "panel\tEGFR\n"
            "panel\tKRAS\n"
            "panel\tBRAF\n"
            # One ambiguous mouse-like entry below threshold
            "panel\tOct4\n"
        )
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write(content)
            f.flush()
            result = parse_geneset_tsv(Path(f.name))
            assert "BRCA1" in result.gene_sets["panel"]
            assert "Oct4" in result.gene_sets["panel"]


# ---------------------------------------------------------------------------
# validate_with_hgnc
# ---------------------------------------------------------------------------


def _mock_hgnc_maps():
    """Return mock HGNC data for testing."""
    canonical = {"BRCA1", "BRCA2", "TP53", "EGFR", "ERCC1"}
    alias_to_canonical = {
        "FANCD1": "BRCA2",
        "ERCC11": "ERCC1",
        "RNF53": "BRCA1",
    }
    canonical_to_aliases = {
        "BRCA2": ["FANCD1"],
        "ERCC1": ["ERCC11"],
        "BRCA1": ["RNF53"],
    }
    return canonical, alias_to_canonical, canonical_to_aliases


class TestValidateWithHgnc:
    @patch("hvantk.utils.gene_aliases._load_hgnc_symbol_maps")
    def test_all_canonical(self, mock_load):
        mock_load.return_value = _mock_hgnc_maps()
        gene_sets = {"panel": ["BRCA1", "TP53", "EGFR"]}

        vr = validate_with_hgnc(gene_sets, "/fake/hgnc.tsv")
        assert vr.recognized == {"BRCA1", "TP53", "EGFR"}
        assert vr.aliases_resolved == {}
        assert vr.unrecognized == set()
        assert vr.gene_sets["panel"] == ["BRCA1", "TP53", "EGFR"]

    @patch("hvantk.utils.gene_aliases._load_hgnc_symbol_maps")
    def test_alias_resolved(self, mock_load):
        mock_load.return_value = _mock_hgnc_maps()
        gene_sets = {"panel": ["FANCD1", "TP53", "ERCC11"]}

        vr = validate_with_hgnc(gene_sets, "/fake/hgnc.tsv")
        assert vr.aliases_resolved == {"FANCD1": "BRCA2", "ERCC11": "ERCC1"}
        assert "BRCA2" in vr.recognized
        assert "ERCC1" in vr.recognized
        assert vr.gene_sets["panel"] == ["BRCA2", "TP53", "ERCC1"]

    @patch("hvantk.utils.gene_aliases._load_hgnc_symbol_maps")
    def test_unrecognized_symbol(self, mock_load):
        mock_load.return_value = _mock_hgnc_maps()
        gene_sets = {"panel": ["BRCA1", "FAKEGENE"]}

        vr = validate_with_hgnc(gene_sets, "/fake/hgnc.tsv")
        assert vr.unrecognized == {"FAKEGENE"}
        assert "FAKEGENE" in vr.gene_sets["panel"]

    @patch("hvantk.utils.gene_aliases._load_hgnc_symbol_maps")
    def test_alias_creates_duplicate(self, mock_load):
        """If user listed both FANCD1 and BRCA2, alias resolution deduplicates."""
        mock_load.return_value = _mock_hgnc_maps()
        gene_sets = {"panel": ["BRCA2", "FANCD1", "TP53"]}

        vr = validate_with_hgnc(gene_sets, "/fake/hgnc.tsv")
        # FANCD1 resolves to BRCA2 which is already in the set.
        assert vr.gene_sets["panel"] == ["BRCA2", "TP53"]

    @patch("hvantk.utils.gene_aliases._load_hgnc_symbol_maps")
    def test_alias_across_multiple_sets(self, mock_load):
        mock_load.return_value = _mock_hgnc_maps()
        gene_sets = {
            "panel_a": ["FANCD1", "TP53"],
            "panel_b": ["ERCC11", "EGFR"],
        }

        vr = validate_with_hgnc(gene_sets, "/fake/hgnc.tsv")
        assert vr.gene_sets["panel_a"] == ["BRCA2", "TP53"]
        assert vr.gene_sets["panel_b"] == ["ERCC1", "EGFR"]
