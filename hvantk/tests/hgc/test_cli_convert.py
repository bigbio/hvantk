"""
Tests for HGC CLI convert commands.
"""

from unittest.mock import patch
from click.testing import CliRunner

from hvantk.commands.hgc.convert_cli import vds2mt, mt2vcf


def test_vds2mt_cli_basic():
    """Test vds2mt command with basic options."""
    runner = CliRunner()
    with patch("hvantk.commands.hgc.convert_cli.convert_vds_to_mt") as mock_convert:
        with patch("hvantk.commands.hgc.convert_cli.validate_input_files") as mock_validate:
            mock_validate.return_value = (True, [])

            result = runner.invoke(
                vds2mt,
                ["--input", "/data.vds", "--output", "/out.mt"],
            )

            assert result.exit_code == 0
            mock_convert.assert_called_once()
            assert "Successfully converted" in result.output


def test_vds2mt_cli_validation_failure():
    """Test vds2mt command with validation failure."""
    runner = CliRunner()
    with patch("hvantk.commands.hgc.convert_cli.validate_input_files") as mock_validate:
        mock_validate.return_value = (False, ["VDS not found"])

        result = runner.invoke(
            vds2mt,
            ["--input", "/missing.vds", "--output", "/out.mt"],
        )

        assert result.exit_code == 1
        assert "Input file validation failed" in result.output
        assert "VDS not found" in result.output


def test_vds2mt_cli_dry_run():
    """Test vds2mt command with dry-run."""
    runner = CliRunner()
    with patch("hvantk.commands.hgc.convert_cli.validate_input_files") as mock_validate:
        mock_validate.return_value = (True, [])

        result = runner.invoke(
            vds2mt,
            [
                "--input", "/data.vds",
                "--output", "/out.mt",
                "--dry-run"
            ],
        )

        assert result.exit_code == 0
        assert "Dry run mode" in result.output
        assert "Adjust genotypes: True" in result.output


def test_vds2mt_cli_skip_validation():
    """Test vds2mt command with --skip-validation."""
    runner = CliRunner()
    with patch("hvantk.commands.hgc.convert_cli.convert_vds_to_mt") as mock_convert:
        with patch("hvantk.commands.hgc.convert_cli.validate_input_files") as mock_validate:
            mock_validate.return_value = (True, [])

            result = runner.invoke(
                vds2mt,
                [
                    "--input", "/data.vds",
                    "--output", "/out.mt",
                    "--skip-validation"
                ],
            )

            assert result.exit_code == 0
            call_kwargs = mock_convert.call_args[1]
            assert call_kwargs["skip_validation"] is True


def test_mt2vcf_cli_basic():
    """Test mt2vcf command with basic options."""
    runner = CliRunner()
    with patch("hvantk.commands.hgc.convert_cli.convert_mt_to_multi_sample_vcf") as mock_convert:
        with patch("hvantk.commands.hgc.convert_cli.validate_input_files") as mock_validate:
            mock_validate.return_value = (True, [])

            result = runner.invoke(
                mt2vcf,
                ["--input", "/data.mt", "--output", "/out.vcf.bgz"],
            )

            assert result.exit_code == 0
            mock_convert.assert_called_once()
            assert "Successfully converted" in result.output


def test_mt2vcf_cli_validation_failure():
    """Test mt2vcf command with validation failure."""
    runner = CliRunner()
    with patch("hvantk.commands.hgc.convert_cli.validate_input_files") as mock_validate:
        mock_validate.return_value = (False, ["MT not found"])

        result = runner.invoke(
            mt2vcf,
            ["--input", "/missing.mt", "--output", "/out.vcf.bgz"],
        )

        assert result.exit_code == 1
        assert "Input file validation failed" in result.output
        assert "MT not found" in result.output


def test_mt2vcf_cli_dry_run():
    """Test mt2vcf command with dry-run."""
    runner = CliRunner()
    with patch("hvantk.commands.hgc.convert_cli.validate_input_files") as mock_validate:
        mock_validate.return_value = (True, [])

        result = runner.invoke(
            mt2vcf,
            [
                "--input", "/data.mt",
                "--output", "/out.vcf.bgz",
                "--dry-run"
            ],
        )

        assert result.exit_code == 0
        assert "Dry run mode" in result.output
        assert "Filter adjusted genotypes: True" in result.output
        assert "Minimum AC: 1" in result.output


def test_mt2vcf_cli_min_ac():
    """Test mt2vcf command with custom min-ac."""
    runner = CliRunner()
    with patch("hvantk.commands.hgc.convert_cli.convert_mt_to_multi_sample_vcf") as mock_convert:
        with patch("hvantk.commands.hgc.convert_cli.validate_input_files") as mock_validate:
            mock_validate.return_value = (True, [])

            result = runner.invoke(
                mt2vcf,
                [
                    "--input", "/data.mt",
                    "--output", "/out.vcf.bgz",
                    "--min-ac", "5"
                ],
            )

            assert result.exit_code == 0
            call_kwargs = mock_convert.call_args[1]
            assert call_kwargs["min_ac"] == 5


def test_mt2vcf_cli_no_filter_adj():
    """Test mt2vcf command with --no-filter-adj."""
    runner = CliRunner()
    with patch("hvantk.commands.hgc.convert_cli.convert_mt_to_multi_sample_vcf") as mock_convert:
        with patch("hvantk.commands.hgc.convert_cli.validate_input_files") as mock_validate:
            mock_validate.return_value = (True, [])

            result = runner.invoke(
                mt2vcf,
                [
                    "--input", "/data.mt",
                    "--output", "/out.vcf.bgz",
                    "--no-filter-adj"
                ],
            )

            assert result.exit_code == 0
            call_kwargs = mock_convert.call_args[1]
            assert call_kwargs["filter_adj_genotypes"] is False
