"""
Tests for HGC CLI QC commands.
"""

from unittest.mock import patch, MagicMock
from click.testing import CliRunner

from hvantk.tools.hgc.qc_cli import (
    compute_qc,
    filter_qc,
    qc_summary,
    plot_qc,
    qc_report,
)

# ===== compute_qc Tests =====


def test_compute_qc_cli_basic():
    """Test compute-qc command with basic options."""
    runner = CliRunner()
    with patch("hvantk.tools.hgc.qc_cli.compute_full_qc") as mock_compute:
        with patch("hvantk.tools.hgc.qc_cli.validate_input_files") as mock_validate:
            with patch("hvantk.tools.hgc.qc_cli.save_qc_metrics") as mock_save:
                with patch("hail.init"):
                    with patch("hail.read_matrix_table") as mock_read:
                        mock_validate.return_value = (True, [])
                        mock_mt = MagicMock()
                        mock_read.return_value = mock_mt
                        mock_qc = MagicMock()
                        mock_compute.return_value = mock_qc
                        mock_save.return_value = {"sample_qc": "/out/sample_qc.csv"}

                        result = runner.invoke(
                            compute_qc,
                            ["--input", "/data.mt", "--output-dir", "/out"],
                        )

                        assert result.exit_code == 0
                        assert "Successfully computed" in result.output
                        mock_compute.assert_called_once()


def test_compute_qc_cli_validation_failure():
    """Test compute-qc command with validation failure."""
    runner = CliRunner()
    with patch("hvantk.tools.hgc.qc_cli.validate_input_files") as mock_validate:
        mock_validate.return_value = (False, ["MT not found"])

        result = runner.invoke(
            compute_qc,
            ["--input", "/missing.mt", "--output-dir", "/out"],
        )

        assert result.exit_code == 1
        assert "Input file validation failed" in result.output


def test_compute_qc_cli_dry_run():
    """Test compute-qc command with dry-run."""
    runner = CliRunner()
    with patch("hvantk.tools.hgc.qc_cli.validate_input_files") as mock_validate:
        mock_validate.return_value = (True, [])

        result = runner.invoke(
            compute_qc,
            ["--input", "/data.mt", "--output-dir", "/out", "--dry-run"],
        )

        assert result.exit_code == 0
        assert "Dry run mode" in result.output
        assert "Sample QC: True" in result.output


def test_compute_qc_cli_no_qc_selected():
    """Test compute-qc command with no QC metrics selected."""
    runner = CliRunner()
    with patch("hvantk.tools.hgc.qc_cli.validate_input_files") as mock_validate:
        mock_validate.return_value = (True, [])

        result = runner.invoke(
            compute_qc,
            [
                "--input",
                "/data.mt",
                "--output-dir",
                "/out",
                "--no-sample-qc",
                "--no-variant-qc",
            ],
        )

        assert result.exit_code == 1
        assert (
            "At least one of --sample-qc or --variant-qc must be enabled"
            in result.output
        )


# ===== filter_qc Tests =====


def test_filter_qc_cli_basic():
    """Test filter-qc command with basic options."""
    runner = CliRunner()
    with patch("hvantk.tools.hgc.qc_cli.filter_samples_by_qc") as mock_filter_s:
        with patch("hvantk.tools.hgc.qc_cli.filter_variants_by_qc") as mock_filter_v:
            with patch(
                "hvantk.tools.hgc.qc_cli.validate_input_files"
            ) as mock_validate:
                with patch("hail.init"):
                    with patch("hail.read_matrix_table") as mock_read:
                        mock_validate.return_value = (True, [])
                        mock_mt = MagicMock()
                        mock_mt.count_cols.return_value = 100
                        mock_mt.count_rows.return_value = 1000
                        mock_read.return_value = mock_mt
                        mock_filter_s.return_value = mock_mt
                        mock_filter_v.return_value = mock_mt

                        result = runner.invoke(
                            filter_qc,
                            ["--input", "/data.mt", "--output", "/out.mt"],
                        )

                        assert result.exit_code == 0
                        assert "Successfully applied QC filters" in result.output


def test_filter_qc_cli_dry_run():
    """Test filter-qc command with dry-run."""
    runner = CliRunner()
    with patch("hvantk.tools.hgc.qc_cli.validate_input_files") as mock_validate:
        mock_validate.return_value = (True, [])

        result = runner.invoke(
            filter_qc,
            ["--input", "/data.mt", "--output", "/out.mt", "--dry-run"],
        )

        assert result.exit_code == 0
        assert "Dry run mode" in result.output
        assert "Min call rate: 0.85" in result.output


def test_filter_qc_cli_custom_thresholds():
    """Test filter-qc command with custom thresholds."""
    runner = CliRunner()
    with patch("hvantk.tools.hgc.qc_cli.filter_samples_by_qc") as mock_filter_s:
        with patch("hvantk.tools.hgc.qc_cli.filter_variants_by_qc") as mock_filter_v:
            with patch(
                "hvantk.tools.hgc.qc_cli.validate_input_files"
            ) as mock_validate:
                with patch("hail.init"):
                    with patch("hail.read_matrix_table") as mock_read:
                        mock_validate.return_value = (True, [])
                        mock_mt = MagicMock()
                        mock_mt.count_cols.return_value = 100
                        mock_mt.count_rows.return_value = 1000
                        mock_read.return_value = mock_mt
                        mock_filter_s.return_value = mock_mt
                        mock_filter_v.return_value = mock_mt

                        result = runner.invoke(
                            filter_qc,
                            [
                                "--input",
                                "/data.mt",
                                "--output",
                                "/out.mt",
                                "--min-ac",
                                "5",
                                "--min-af",
                                "0.01",
                            ],
                        )

                        assert result.exit_code == 0
                        # Check that custom thresholds were passed
                        call_kwargs = mock_filter_v.call_args[1]
                        assert call_kwargs["min_ac"] == 5
                        assert call_kwargs["min_af"] == 0.01


# ===== qc_summary Tests =====


def test_qc_summary_cli_basic():
    """Test qc-summary command with basic options."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        import os
        import pandas as pd

        # Create test directory and QC files
        os.makedirs("qc_dir")
        sample_df = pd.DataFrame({"call_rate": [0.9, 0.95], "mean_dp": [30, 35]})
        sample_df.to_csv("qc_dir/sample_qc.csv", index=False)

        with patch("hvantk.algorithms.hgc.qc.get_qc_summary_stats") as mock_stats:
            mock_stats.return_value = pd.DataFrame()

            result = runner.invoke(
                qc_summary,
                ["--qc-dir", "qc_dir"],
            )

            assert result.exit_code == 0
            assert "QC summary completed successfully" in result.output


def test_qc_summary_cli_no_files():
    """Test qc-summary command with no QC files found."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        import os

        os.makedirs("empty_dir")

        result = runner.invoke(
            qc_summary,
            ["--qc-dir", "empty_dir"],
        )

        assert result.exit_code == 1
        assert "No QC files found" in result.output


def test_qc_summary_cli_output_json():
    """Test qc-summary command with JSON output."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        import os
        import pandas as pd

        os.makedirs("qc_dir")
        sample_df = pd.DataFrame({"call_rate": [0.9, 0.95]})
        sample_df.to_csv("qc_dir/sample_qc.csv", index=False)

        with patch("hvantk.algorithms.hgc.qc.get_qc_summary_stats") as mock_stats:
            mock_stats.return_value = pd.DataFrame()

            result = runner.invoke(
                qc_summary,
                ["--qc-dir", "qc_dir", "--output", "summary.json", "--format", "json"],
            )

            assert result.exit_code == 0
            assert os.path.exists("summary.json")


# ===== plot_qc Tests =====


def test_plot_qc_cli_basic():
    """Test plot-qc command with basic options."""
    runner = CliRunner()
    with patch(
        "hvantk.tools.hgc.qc_cli.check_path_exists_and_readable"
    ) as mock_check:
        with patch("hail.init"):
            with patch("hail.read_matrix_table") as mock_read:
                mock_check.return_value = True
                mock_mt = MagicMock()
                mock_mt.col = {"sample_qc": MagicMock()}
                mock_read.return_value = mock_mt

                with patch("hvantk.algorithms.hgc.qc.QCMetrics") as mock_qc_metrics:
                    mock_qc_instance = MagicMock()
                    mock_qc_instance.has_sample_qc = True
                    mock_qc_instance.has_variant_qc = False

                    # Mock plot function to create the file
                    def mock_plot_overview(save_path=None, **kwargs):
                        if save_path:
                            from pathlib import Path

                            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
                            Path(save_path).write_text("mock plot")
                        return MagicMock()

                    mock_qc_instance.plot_sample_overview.side_effect = (
                        mock_plot_overview
                    )
                    mock_qc_metrics.return_value = mock_qc_instance

                    with runner.isolated_filesystem():
                        result = runner.invoke(
                            plot_qc,
                            ["--input", "/data.mt", "--output-dir", "plots"],
                        )

                        assert result.exit_code == 0
                        assert "QC plotting completed successfully" in result.output


def test_plot_qc_cli_dry_run():
    """Test plot-qc command with dry-run."""
    runner = CliRunner()
    with patch(
        "hvantk.tools.hgc.qc_cli.check_path_exists_and_readable"
    ) as mock_check:
        mock_check.return_value = True

        result = runner.invoke(
            plot_qc,
            ["--input", "/data.mt", "--output-dir", "plots", "--dry-run"],
        )

        assert result.exit_code == 0
        assert "Dry run mode" in result.output
        assert "Plot type: overview" in result.output


def test_plot_qc_cli_no_qc_annotations():
    """Test plot-qc command with MT lacking QC annotations."""
    runner = CliRunner()
    with patch(
        "hvantk.tools.hgc.qc_cli.check_path_exists_and_readable"
    ) as mock_check:
        with patch("hail.init"):
            with patch("hail.read_matrix_table") as mock_read:
                mock_check.return_value = True
                mock_mt = MagicMock()
                mock_mt.col = {}
                mock_mt.row = {}
                mock_read.return_value = mock_mt

                result = runner.invoke(
                    plot_qc,
                    ["--input", "/data.mt", "--output-dir", "plots"],
                )

                assert result.exit_code == 1
                assert "No QC annotations found" in result.output


# ===== qc_report Tests =====


def test_qc_report_cli_basic():
    """Test qc-report command with basic options."""
    runner = CliRunner()
    with patch(
        "hvantk.tools.hgc.qc_cli.check_path_exists_and_readable"
    ) as mock_check:
        with patch("hail.init"):
            with patch("hail.read_matrix_table") as mock_read:
                mock_check.return_value = True
                mock_mt = MagicMock()
                mock_mt.col = {"sample_qc": MagicMock()}
                mock_read.return_value = mock_mt

                with patch("hvantk.algorithms.hgc.qc.QCMetrics") as mock_qc_metrics:
                    from pathlib import Path

                    mock_qc_instance = MagicMock()
                    mock_qc_instance.has_sample_qc = True
                    mock_qc_instance.has_variant_qc = False
                    mock_qc_instance.generate_html_report.return_value = Path(
                        "report.html"
                    )
                    mock_qc_metrics.return_value = mock_qc_instance

                    with runner.isolated_filesystem():
                        # Create dummy report file for stat check
                        Path("report.html").write_text("<html></html>")

                        result = runner.invoke(
                            qc_report,
                            ["--input", "/data.mt", "--output", "report.html"],
                        )

                        assert result.exit_code == 0
                        assert "QC HTML report generated successfully" in result.output


def test_qc_report_cli_dry_run():
    """Test qc-report command with dry-run."""
    runner = CliRunner()
    with patch(
        "hvantk.tools.hgc.qc_cli.check_path_exists_and_readable"
    ) as mock_check:
        mock_check.return_value = True

        result = runner.invoke(
            qc_report,
            ["--input", "/data.mt", "--output", "report.html", "--dry-run"],
        )

        assert result.exit_code == 0
        assert "Dry run mode" in result.output
        assert "Report title: Quality Control Report" in result.output


def test_qc_report_cli_custom_title():
    """Test qc-report command with custom title."""
    runner = CliRunner()
    with patch(
        "hvantk.tools.hgc.qc_cli.check_path_exists_and_readable"
    ) as mock_check:
        with patch("hail.init"):
            with patch("hail.read_matrix_table") as mock_read:
                mock_check.return_value = True
                mock_mt = MagicMock()
                mock_mt.col = {"sample_qc": MagicMock()}
                mock_read.return_value = mock_mt

                with patch("hvantk.algorithms.hgc.qc.QCMetrics") as mock_qc_metrics:
                    from pathlib import Path

                    mock_qc_instance = MagicMock()
                    mock_qc_instance.has_sample_qc = True
                    mock_qc_instance.has_variant_qc = False
                    mock_qc_instance.generate_html_report.return_value = Path(
                        "custom_report.html"
                    )
                    mock_qc_metrics.return_value = mock_qc_instance

                    with runner.isolated_filesystem():
                        Path("custom_report.html").write_text("<html></html>")

                        result = runner.invoke(
                            qc_report,
                            [
                                "--input",
                                "/data.mt",
                                "--output",
                                "custom_report.html",
                                "--title",
                                "My Custom QC Report",
                            ],
                        )

                        assert result.exit_code == 0
                        # Check that custom title was passed
                        call_kwargs = mock_qc_instance.generate_html_report.call_args[1]
                        assert call_kwargs["title"] == "My Custom QC Report"
