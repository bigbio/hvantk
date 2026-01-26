"""
Tests for PSROC CLI module.

This module tests the CLI entry point for the PSROC pipeline, including
argument parsing, validation, and error handling.
"""

import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

from click.testing import CliRunner

from hvantk.commands.psroc_cli import psroc_cmd, _parse_comma_separated


class TestParseCommaSeparated:
    """Test the comma-separated string parser."""

    def test_parse_single_value(self):
        """Test parsing single value."""
        result = _parse_comma_separated("CADD_phred")
        assert result == ["CADD_phred"]

    def test_parse_multiple_values(self):
        """Test parsing multiple comma-separated values."""
        result = _parse_comma_separated("CADD_phred,REVEL_score,MetaLR_score")
        assert result == ["CADD_phred", "REVEL_score", "MetaLR_score"]

    def test_parse_with_spaces(self):
        """Test parsing values with surrounding spaces."""
        result = _parse_comma_separated("CADD_phred , REVEL_score , MetaLR_score")
        assert result == ["CADD_phred", "REVEL_score", "MetaLR_score"]

    def test_parse_empty_string(self):
        """Test parsing empty string."""
        result = _parse_comma_separated("")
        assert result == []

    def test_parse_none(self):
        """Test parsing None returns empty list."""
        result = _parse_comma_separated(None)
        assert result == []

    def test_parse_with_empty_elements(self):
        """Test parsing string with empty elements."""
        result = _parse_comma_separated("CADD_phred,,REVEL_score")
        assert result == ["CADD_phred", "REVEL_score"]


class TestPSROCCLIHelp:
    """Test CLI help output."""

    def test_help_displays(self):
        """Test that help text displays correctly."""
        runner = CliRunner()
        result = runner.invoke(psroc_cmd, ["--help"])

        assert result.exit_code == 0
        assert "PSROC: Prediction Score ROC Analysis" in result.output
        assert "--clinvar-ht" in result.output
        assert "--dbnsfp-ht" in result.output
        assert "--scores" in result.output
        assert "--output-dir" in result.output

    def test_help_shows_examples(self):
        """Test that help includes usage examples."""
        runner = CliRunner()
        result = runner.invoke(psroc_cmd, ["--help"])

        assert "Examples:" in result.output
        assert "hvantk psroc" in result.output


class TestPSROCCLIValidation:
    """Test CLI argument validation."""

    def test_missing_required_args(self):
        """Test that missing required arguments shows error."""
        runner = CliRunner()
        result = runner.invoke(psroc_cmd, [])

        assert result.exit_code != 0
        assert "Missing option" in result.output or "Error" in result.output

    def test_missing_clinvar_ht(self):
        """Test error when --clinvar-ht not provided."""
        runner = CliRunner()
        result = runner.invoke(
            psroc_cmd,
            [
                "--genes",
                "BRCA1",
                "--scores",
                "CADD_phred",
                "--output-dir",
                "/fake/output",
            ],
        )

        assert result.exit_code != 0
        # Click shows "Missing option" for required options
        assert "Missing option" in result.output or "clinvar" in result.output.lower()

    def test_missing_scores(self):
        """Test error when --scores not provided."""
        runner = CliRunner()
        result = runner.invoke(
            psroc_cmd,
            [
                "--genes",
                "BRCA1",
                "--output-dir",
                "/fake/output",
            ],
        )

        assert result.exit_code != 0
        assert "Missing option" in result.output or "scores" in result.output.lower()

    def test_missing_output_dir(self):
        """Test error when --output-dir not provided."""
        runner = CliRunner()
        result = runner.invoke(
            psroc_cmd,
            [
                "--genes",
                "BRCA1",
                "--scores",
                "CADD_phred",
            ],
        )

        assert result.exit_code != 0
        assert "Missing option" in result.output or "output" in result.output.lower()

    def test_no_variant_source(self):
        """Test error when no variant source provided."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create fake table directories
            clinvar_path = Path(tmpdir) / "clinvar.ht"
            dbnsfp_path = Path(tmpdir) / "dbnsfp.ht"
            clinvar_path.mkdir()
            dbnsfp_path.mkdir()

            result = runner.invoke(
                psroc_cmd,
                [
                    "--clinvar-ht",
                    str(clinvar_path),
                    "--dbnsfp-ht",
                    str(dbnsfp_path),
                    "--scores",
                    "CADD_phred",
                    "--output-dir",
                    tmpdir,
                ],
            )

            assert result.exit_code != 0
            assert "genes" in result.output.lower() or "variants" in result.output.lower()

    def test_multiple_variant_sources(self):
        """Test error when multiple variant sources provided."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create fake table directories and variant file
            clinvar_path = Path(tmpdir) / "clinvar.ht"
            dbnsfp_path = Path(tmpdir) / "dbnsfp.ht"
            variants_file = Path(tmpdir) / "variants.txt"
            clinvar_path.mkdir()
            dbnsfp_path.mkdir()
            variants_file.write_text("chr1:12345:A:T\n")

            result = runner.invoke(
                psroc_cmd,
                [
                    "--genes",
                    "BRCA1",
                    "--variants",
                    str(variants_file),  # Both genes and variants
                    "--clinvar-ht",
                    str(clinvar_path),
                    "--dbnsfp-ht",
                    str(dbnsfp_path),
                    "--scores",
                    "CADD_phred",
                    "--output-dir",
                    tmpdir,
                ],
            )

            assert result.exit_code != 0
            assert "multiple" in result.output.lower() or "cannot" in result.output.lower()

    def test_invalid_reference_genome(self):
        """Test error with invalid reference genome."""
        runner = CliRunner()
        result = runner.invoke(
            psroc_cmd,
            [
                "--genes",
                "BRCA1",
                "--clinvar-ht",
                "/fake/clinvar.ht",
                "--dbnsfp-ht",
                "/fake/dbnsfp.ht",
                "--scores",
                "CADD_phred",
                "--output-dir",
                "/fake/output",
                "--reference-genome",
                "hg19",  # Invalid
            ],
        )

        assert result.exit_code != 0

    def test_invalid_threshold_method(self):
        """Test error with invalid threshold method."""
        runner = CliRunner()
        result = runner.invoke(
            psroc_cmd,
            [
                "--genes",
                "BRCA1",
                "--clinvar-ht",
                "/fake/clinvar.ht",
                "--dbnsfp-ht",
                "/fake/dbnsfp.ht",
                "--scores",
                "CADD_phred",
                "--output-dir",
                "/fake/output",
                "--threshold-method",
                "invalid",  # Invalid
            ],
        )

        assert result.exit_code != 0


class TestPSROCCLIDryRun:
    """Test CLI dry-run functionality."""

    @patch("hvantk.psroc.pipeline.PSROCPipeline")
    @patch("hvantk.psroc.pipeline.PSROCConfig")
    def test_dry_run_shows_plan(self, mock_config_cls, mock_pipeline_cls):
        """Test that --dry-run shows execution plan without running."""
        runner = CliRunner()

        # Set up mocks
        mock_config = MagicMock()
        mock_config.validate.return_value = []
        mock_config_cls.return_value = mock_config

        mock_pipeline = MagicMock()
        mock_pipeline_cls.return_value = mock_pipeline

        with tempfile.TemporaryDirectory() as tmpdir:
            clinvar_path = Path(tmpdir) / "clinvar.ht"
            dbnsfp_path = Path(tmpdir) / "dbnsfp.ht"
            clinvar_path.mkdir()
            dbnsfp_path.mkdir()

            runner.invoke(
                psroc_cmd,
                [
                    "--genes",
                    "BRCA1",
                    "--clinvar-ht",
                    str(clinvar_path),
                    "--dbnsfp-ht",
                    str(dbnsfp_path),
                    "--scores",
                    "CADD_phred",
                    "--output-dir",
                    tmpdir,
                    "--dry-run",
                ],
            )

            # Should call show_plan but not run
            mock_pipeline.show_plan.assert_called_once()
            mock_pipeline.run.assert_not_called()


class TestPSROCCLIOptions:
    """Test CLI option handling."""

    def test_genes_parsing(self):
        """Test that genes are correctly parsed from comma-separated string."""
        # Test the parsing function directly
        genes = _parse_comma_separated("BRCA1,BRCA2,TP53")
        assert genes == ["BRCA1", "BRCA2", "TP53"]

    def test_scores_parsing(self):
        """Test that scores are correctly parsed."""
        scores = _parse_comma_separated("CADD_phred,REVEL_score")
        assert scores == ["CADD_phred", "REVEL_score"]

    @patch("hvantk.psroc.pipeline.PSROCPipeline")
    @patch("hvantk.psroc.pipeline.PSROCConfig")
    def test_max_missingness_passed_to_config(self, mock_config_cls, mock_pipeline_cls):
        """Test that max_missingness is passed to config."""
        runner = CliRunner()

        mock_config = MagicMock()
        mock_config.validate.return_value = []
        mock_config_cls.return_value = mock_config

        mock_pipeline = MagicMock()
        mock_pipeline_cls.return_value = mock_pipeline

        with tempfile.TemporaryDirectory() as tmpdir:
            clinvar_path = Path(tmpdir) / "clinvar.ht"
            dbnsfp_path = Path(tmpdir) / "dbnsfp.ht"
            clinvar_path.mkdir()
            dbnsfp_path.mkdir()

            runner.invoke(
                psroc_cmd,
                [
                    "--genes",
                    "BRCA1",
                    "--clinvar-ht",
                    str(clinvar_path),
                    "--dbnsfp-ht",
                    str(dbnsfp_path),
                    "--scores",
                    "CADD_phred",
                    "--output-dir",
                    tmpdir,
                    "--max-missingness",
                    "0.5",
                    "--dry-run",
                ],
            )

            # Check that config was created with max_missingness=0.5
            call_kwargs = mock_config_cls.call_args[1]
            assert call_kwargs["max_missingness"] == 0.5

    @patch("hvantk.psroc.pipeline.PSROCPipeline")
    @patch("hvantk.psroc.pipeline.PSROCConfig")
    def test_no_plots_flag(self, mock_config_cls, mock_pipeline_cls):
        """Test that --no-plots sets generate_plots=False."""
        runner = CliRunner()

        mock_config = MagicMock()
        mock_config.validate.return_value = []
        mock_config_cls.return_value = mock_config

        mock_pipeline = MagicMock()
        mock_pipeline_cls.return_value = mock_pipeline

        with tempfile.TemporaryDirectory() as tmpdir:
            clinvar_path = Path(tmpdir) / "clinvar.ht"
            dbnsfp_path = Path(tmpdir) / "dbnsfp.ht"
            clinvar_path.mkdir()
            dbnsfp_path.mkdir()

            runner.invoke(
                psroc_cmd,
                [
                    "--genes",
                    "BRCA1",
                    "--clinvar-ht",
                    str(clinvar_path),
                    "--dbnsfp-ht",
                    str(dbnsfp_path),
                    "--scores",
                    "CADD_phred",
                    "--output-dir",
                    tmpdir,
                    "--no-plots",
                    "--dry-run",
                ],
            )

            call_kwargs = mock_config_cls.call_args[1]
            assert call_kwargs["generate_plots"] is False

    @patch("hvantk.psroc.pipeline.PSROCPipeline")
    @patch("hvantk.psroc.pipeline.PSROCConfig")
    def test_export_tsv_flag(self, mock_config_cls, mock_pipeline_cls):
        """Test that --export-tsv sets export_tsv=True."""
        runner = CliRunner()

        mock_config = MagicMock()
        mock_config.validate.return_value = []
        mock_config_cls.return_value = mock_config

        mock_pipeline = MagicMock()
        mock_pipeline_cls.return_value = mock_pipeline

        with tempfile.TemporaryDirectory() as tmpdir:
            clinvar_path = Path(tmpdir) / "clinvar.ht"
            dbnsfp_path = Path(tmpdir) / "dbnsfp.ht"
            clinvar_path.mkdir()
            dbnsfp_path.mkdir()

            runner.invoke(
                psroc_cmd,
                [
                    "--genes",
                    "BRCA1",
                    "--clinvar-ht",
                    str(clinvar_path),
                    "--dbnsfp-ht",
                    str(dbnsfp_path),
                    "--scores",
                    "CADD_phred",
                    "--output-dir",
                    tmpdir,
                    "--export-tsv",
                    "--dry-run",
                ],
            )

            call_kwargs = mock_config_cls.call_args[1]
            assert call_kwargs["export_tsv"] is True


class TestPSROCCLIVariantSources:
    """Test different variant source options."""

    @patch("hvantk.psroc.pipeline.PSROCPipeline")
    @patch("hvantk.psroc.pipeline.PSROCConfig")
    def test_genes_option(self, mock_config_cls, mock_pipeline_cls):
        """Test using --genes option."""
        runner = CliRunner()

        mock_config = MagicMock()
        mock_config.validate.return_value = []
        mock_config_cls.return_value = mock_config

        mock_pipeline = MagicMock()
        mock_pipeline_cls.return_value = mock_pipeline

        with tempfile.TemporaryDirectory() as tmpdir:
            clinvar_path = Path(tmpdir) / "clinvar.ht"
            dbnsfp_path = Path(tmpdir) / "dbnsfp.ht"
            clinvar_path.mkdir()
            dbnsfp_path.mkdir()

            runner.invoke(
                psroc_cmd,
                [
                    "--genes",
                    "BRCA1,BRCA2",
                    "--clinvar-ht",
                    str(clinvar_path),
                    "--dbnsfp-ht",
                    str(dbnsfp_path),
                    "--scores",
                    "CADD_phred",
                    "--output-dir",
                    tmpdir,
                    "--dry-run",
                ],
            )

            call_kwargs = mock_config_cls.call_args[1]
            assert call_kwargs["genes"] == ["BRCA1", "BRCA2"]
            assert call_kwargs["genes_file"] is None
            assert call_kwargs["variants_path"] is None

    @patch("hvantk.psroc.pipeline.PSROCPipeline")
    @patch("hvantk.psroc.pipeline.PSROCConfig")
    def test_genes_file_option(self, mock_config_cls, mock_pipeline_cls):
        """Test using --genes-file option."""
        runner = CliRunner()

        mock_config = MagicMock()
        mock_config.validate.return_value = []
        mock_config_cls.return_value = mock_config

        mock_pipeline = MagicMock()
        mock_pipeline_cls.return_value = mock_pipeline

        with tempfile.TemporaryDirectory() as tmpdir:
            clinvar_path = Path(tmpdir) / "clinvar.ht"
            dbnsfp_path = Path(tmpdir) / "dbnsfp.ht"
            genes_file = Path(tmpdir) / "genes.txt"
            clinvar_path.mkdir()
            dbnsfp_path.mkdir()
            genes_file.write_text("BRCA1\nBRCA2\n")

            runner.invoke(
                psroc_cmd,
                [
                    "--genes-file",
                    str(genes_file),
                    "--clinvar-ht",
                    str(clinvar_path),
                    "--dbnsfp-ht",
                    str(dbnsfp_path),
                    "--scores",
                    "CADD_phred",
                    "--output-dir",
                    tmpdir,
                    "--dry-run",
                ],
            )

            call_kwargs = mock_config_cls.call_args[1]
            assert call_kwargs["genes"] is None
            assert call_kwargs["genes_file"] == str(genes_file)
            assert call_kwargs["variants_path"] is None

    @patch("hvantk.psroc.pipeline.PSROCPipeline")
    @patch("hvantk.psroc.pipeline.PSROCConfig")
    def test_variants_option(self, mock_config_cls, mock_pipeline_cls):
        """Test using --variants option."""
        runner = CliRunner()

        mock_config = MagicMock()
        mock_config.validate.return_value = []
        mock_config_cls.return_value = mock_config

        mock_pipeline = MagicMock()
        mock_pipeline_cls.return_value = mock_pipeline

        with tempfile.TemporaryDirectory() as tmpdir:
            clinvar_path = Path(tmpdir) / "clinvar.ht"
            dbnsfp_path = Path(tmpdir) / "dbnsfp.ht"
            variants_file = Path(tmpdir) / "variants.txt"
            clinvar_path.mkdir()
            dbnsfp_path.mkdir()
            variants_file.write_text("chr1:12345:A:T\nchr2:67890:G:C\n")

            runner.invoke(
                psroc_cmd,
                [
                    "--variants",
                    str(variants_file),
                    "--clinvar-ht",
                    str(clinvar_path),
                    "--dbnsfp-ht",
                    str(dbnsfp_path),
                    "--scores",
                    "CADD_phred",
                    "--output-dir",
                    tmpdir,
                    "--dry-run",
                ],
            )

            call_kwargs = mock_config_cls.call_args[1]
            assert call_kwargs["genes"] is None
            assert call_kwargs["genes_file"] is None
            assert call_kwargs["variants_path"] == str(variants_file)


if __name__ == "__main__":
    print("Running PSROC CLI tests...")

    print("\n1. Testing comma-separated parser")
    test_parser = TestParseCommaSeparated()
    test_parser.test_parse_single_value()
    print("  ✓ Single value parsing works")
    test_parser.test_parse_multiple_values()
    print("  ✓ Multiple value parsing works")
    test_parser.test_parse_with_spaces()
    print("  ✓ Space handling works")

    print("\n2. Testing CLI help")
    test_help = TestPSROCCLIHelp()
    test_help.test_help_displays()
    print("  ✓ Help displays correctly")

    print("\n3. Testing CLI validation")
    test_validation = TestPSROCCLIValidation()
    test_validation.test_missing_required_args()
    print("  ✓ Missing required args detected")
    test_validation.test_no_variant_source()
    print("  ✓ No variant source detected")
    test_validation.test_multiple_variant_sources()
    print("  ✓ Multiple variant sources detected")

    print("\n✅ All tests passed!")
