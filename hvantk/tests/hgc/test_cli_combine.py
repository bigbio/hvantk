"""
Tests for HGC CLI combine commands.
"""

from unittest.mock import patch
from click.testing import CliRunner

from hvantk.tools.hgc.combine_cli import gvcf_combine, vds_combine


def test_gvcf_combine_cli_basic():
    """Test gvcf-combine command with basic options."""
    runner = CliRunner()
    with patch("hvantk.tools.hgc.combine_cli.combine_gvcfs") as mock_combine:
        with patch(
            "hvantk.tools.hgc.combine_cli.validate_output_path"
        ) as mock_validate:
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
    with patch("hvantk.tools.hgc.combine_cli.validate_output_path") as mock_validate:
        mock_validate.return_value = True

        result = runner.invoke(
            gvcf_combine,
            ["--gvcf-dir", "/path/to/gvcfs", "--output", "/out.vds", "--dry-run"],
        )

        assert result.exit_code == 0
        assert "Dry run mode" in result.output
        assert "GVCF directory: /path/to/gvcfs" in result.output


def test_gvcf_combine_cli_with_vds_paths():
    """Test gvcf-combine command with VDS paths."""
    runner = CliRunner()
    with patch("hvantk.tools.hgc.combine_cli.combine_gvcfs") as mock_combine:
        with patch(
            "hvantk.tools.hgc.combine_cli.validate_output_path"
        ) as mock_validate:
            mock_validate.return_value = True

            result = runner.invoke(
                gvcf_combine,
                [
                    "--gvcf-dir",
                    "/gvcfs",
                    "--vds-paths",
                    "/vds1.vds",
                    "--vds-paths",
                    "/vds2.vds",
                    "--output",
                    "/out.vds",
                ],
            )

            assert result.exit_code == 0
            mock_combine.assert_called_once()
            # Check vdses parameter
            call_kwargs = mock_combine.call_args[1]
            assert call_kwargs["vdses"] == ["/vds1.vds", "/vds2.vds"]


def test_gvcf_combine_cli_default_forwards_empty_kwargs():
    """Default invocation must forward no combiner kwargs (Hail's genome default applies).

    Regression guard: previously ``kwargs={}`` was hardcoded, so the partitioning was
    unreachable. Exposing the knobs must not change the default behaviour.
    """
    runner = CliRunner()
    with patch("hvantk.tools.hgc.combine_cli.combine_gvcfs") as mock_combine:
        with patch(
            "hvantk.tools.hgc.combine_cli.validate_output_path"
        ) as mock_validate:
            mock_validate.return_value = True

            result = runner.invoke(
                gvcf_combine,
                ["--gvcf-dir", "/gvcfs", "--output", "/out.vds"],
            )

            assert result.exit_code == 0
            assert mock_combine.call_args[1]["kwargs"] == {}


def test_gvcf_combine_cli_import_interval_size():
    """--import-interval-size is forwarded to the combiner (the parallelism knob)."""
    runner = CliRunner()
    with patch("hvantk.tools.hgc.combine_cli.combine_gvcfs") as mock_combine:
        with patch(
            "hvantk.tools.hgc.combine_cli.validate_output_path"
        ) as mock_validate:
            mock_validate.return_value = True

            result = runner.invoke(
                gvcf_combine,
                [
                    "--gvcf-dir",
                    "/gvcfs",
                    "--output",
                    "/out.vds",
                    "--import-interval-size",
                    "600000",
                ],
            )

            assert result.exit_code == 0
            assert mock_combine.call_args[1]["kwargs"] == {
                "import_interval_size": 600000
            }


def test_gvcf_combine_cli_tree_merge_options():
    """--gvcf-batch-size and --branch-factor are forwarded to the combiner."""
    runner = CliRunner()
    with patch("hvantk.tools.hgc.combine_cli.combine_gvcfs") as mock_combine:
        with patch(
            "hvantk.tools.hgc.combine_cli.validate_output_path"
        ) as mock_validate:
            mock_validate.return_value = True

            result = runner.invoke(
                gvcf_combine,
                [
                    "--gvcf-dir",
                    "/gvcfs",
                    "--output",
                    "/out.vds",
                    "--gvcf-batch-size",
                    "100",
                    "--branch-factor",
                    "50",
                ],
            )

            assert result.exit_code == 0
            assert mock_combine.call_args[1]["kwargs"] == {
                "gvcf_batch_size": 100,
                "branch_factor": 50,
            }


def test_gvcf_combine_cli_exome_default_intervals():
    """--use-exome-default-intervals is forwarded (previously unreachable from the CLI)."""
    runner = CliRunner()
    with patch("hvantk.tools.hgc.combine_cli.combine_gvcfs") as mock_combine:
        with patch(
            "hvantk.tools.hgc.combine_cli.validate_output_path"
        ) as mock_validate:
            mock_validate.return_value = True

            result = runner.invoke(
                gvcf_combine,
                [
                    "--gvcf-dir",
                    "/gvcfs",
                    "--output",
                    "/out.vds",
                    "--use-exome-default-intervals",
                ],
            )

            assert result.exit_code == 0
            assert mock_combine.call_args[1]["kwargs"] == {
                "use_exome_default_intervals": True
            }


def test_gvcf_combine_cli_partitioning_options_mutually_exclusive():
    """Colliding partitioning options must fail loudly.

    Hail itself only emits a warning and silently picks one, which is easy to miss in a
    long combiner log, so the CLI rejects the combination up front.
    """
    runner = CliRunner()
    with patch("hvantk.tools.hgc.combine_cli.validate_output_path") as mock_validate:
        mock_validate.return_value = True

        result = runner.invoke(
            gvcf_combine,
            [
                "--gvcf-dir",
                "/gvcfs",
                "--output",
                "/out.vds",
                "--import-interval-size",
                "600000",
                "--use-genome-default-intervals",
            ],
        )

        assert result.exit_code == 1
        assert "mutually exclusive" in result.output


def test_gvcf_combine_cli_dry_run_shows_combiner_options():
    """Dry-run surfaces the combiner options that would be used."""
    runner = CliRunner()
    with patch("hvantk.tools.hgc.combine_cli.validate_output_path") as mock_validate:
        mock_validate.return_value = True

        result = runner.invoke(
            gvcf_combine,
            [
                "--gvcf-dir",
                "/gvcfs",
                "--output",
                "/out.vds",
                "--import-interval-size",
                "300000",
                "--dry-run",
            ],
        )

        assert result.exit_code == 0
        assert "import_interval_size" in result.output
        assert "300000" in result.output


def test_vds_combine_cli_basic():
    """Test vds-combine command with basic options."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        import os

        os.makedirs("input_dir")

        with patch("hvantk.tools.hgc.combine_cli.combine_vdses") as mock_combine:
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
            ["--input-dir", "input_dir", "--output", "/out.vds", "--dry-run"],
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

        with patch("hvantk.tools.hgc.combine_cli.combine_vdses") as mock_combine:
            result = runner.invoke(
                vds_combine,
                ["--input-dir", "input_dir", "--output", "/out.vds", "--no-validate"],
            )

            assert result.exit_code == 0
            call_kwargs = mock_combine.call_args[1]
            assert call_kwargs["validate"] is False
