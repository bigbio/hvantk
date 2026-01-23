"""
Tests for HGC CLI combine commands.
"""

from unittest.mock import patch, MagicMock
from click.testing import CliRunner
import pytest

from hvantk.commands.hgc.combine_cli import gvcf_combine, vds_combine


def test_gvcf_combine_cli_basic():
    """Test gvcf-combine command with basic options."""
    runner = CliRunner()
    with patch("hvantk.commands.hgc.combine_cli.combine_gvcfs") as mock_combine:
        with patch("hvantk.commands.hgc.combine_cli.validate_output_path") as mock_validate:
            mock_validate.return_value = True

            result = runner.invoke(
                gvcf_combine,
                ["--gvcf-dir", "/path/to/gvcfs", "--output", "/out.vds"],
            )

            assert result.exit_code == 0
            mock_combine.assert_called_once()
            assert "Successfully combined" in result.output


def test_gvcf_combine_cli_missing_inputs():
    """Test gvcf-combine command with missing inputs."""
    runner = CliRunner()
    result = runner.invoke(
        gvcf_combine,
        ["--output", "/out.vds"],
    )

    assert result.exit_code == 1
    assert "Either --gvcf-dir or --vds-paths must be provided" in result.output


def test_gvcf_combine_cli_dry_run():
    """Test gvcf-combine command with dry-run."""
    runner = CliRunner()
    with patch("hvantk.commands.hgc.combine_cli.validate_output_path") as mock_validate:
        mock_validate.return_value = True

        result = runner.invoke(
            gvcf_combine,
            [
                "--gvcf-dir", "/path/to/gvcfs",
                "--output", "/out.vds",
                "--dry-run"
            ],
        )

        assert result.exit_code == 0
        assert "Dry run mode" in result.output
        assert "GVCF directory: /path/to/gvcfs" in result.output


def test_gvcf_combine_cli_with_vds_paths():
    """Test gvcf-combine command with VDS paths."""
    runner = CliRunner()
    with patch("hvantk.commands.hgc.combine_cli.combine_gvcfs") as mock_combine:
        with patch("hvantk.commands.hgc.combine_cli.validate_output_path") as mock_validate:
            mock_validate.return_value = True

            result = runner.invoke(
                gvcf_combine,
                [
                    "--gvcf-dir", "/gvcfs",
                    "--vds-paths", "/vds1.vds",
                    "--vds-paths", "/vds2.vds",
                    "--output", "/out.vds"
                ],
            )

            assert result.exit_code == 0
            mock_combine.assert_called_once()
            # Check vdses parameter
            call_kwargs = mock_combine.call_args[1]
            assert call_kwargs["vdses"] == ["/vds1.vds", "/vds2.vds"]


def test_vds_combine_cli_basic():
    """Test vds-combine command with basic options."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        import os
        os.makedirs("input_dir")

        with patch("hvantk.commands.hgc.combine_cli.combine_vdses") as mock_combine:
            result = runner.invoke(
                vds_combine,
                ["--input-dir", "input_dir", "--output", "/out.vds"],
            )

            assert result.exit_code == 0
            mock_combine.assert_called_once()
            assert "Successfully combined" in result.output


def test_vds_combine_cli_missing_input_dir():
    """Test vds-combine command with missing input directory."""
    runner = CliRunner()
    result = runner.invoke(
        vds_combine,
        ["--input-dir", "/nonexistent", "--output", "/out.vds"],
    )

    assert result.exit_code == 1
    assert "Input directory not found" in result.output


def test_vds_combine_cli_dry_run():
    """Test vds-combine command with dry-run."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        import os
        os.makedirs("input_dir")

        result = runner.invoke(
            vds_combine,
            [
                "--input-dir", "input_dir",
                "--output", "/out.vds",
                "--dry-run"
            ],
        )

        assert result.exit_code == 0
        assert "Dry run mode" in result.output
        assert "Validate: True" in result.output


def test_vds_combine_cli_no_validate():
    """Test vds-combine command with --no-validate."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        import os
        os.makedirs("input_dir")

        with patch("hvantk.commands.hgc.combine_cli.combine_vdses") as mock_combine:
            result = runner.invoke(
                vds_combine,
                [
                    "--input-dir", "input_dir",
                    "--output", "/out.vds",
                    "--no-validate"
                ],
            )

            assert result.exit_code == 0
            call_kwargs = mock_combine.call_args[1]
            assert call_kwargs["validate"] is False
