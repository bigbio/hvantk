"""
Tests for ancestry inference CLI module.

This module tests the CLI entry point for the ancestry inference pipeline,
including argument parsing, validation, and error handling.
"""

import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest
from click.testing import CliRunner

from hvantk.commands.ancestry_cli import ancestry_inference_cmd


class TestAncestryCLIHelp:
    """Test CLI help output."""

    def test_help_displays(self):
        """Test that help text displays correctly."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert result.exit_code == 0
        assert "ancestry" in result.output.lower()
        assert "--query-mt" in result.output
        assert "--reference-mt" in result.output
        assert "--ancestry-col" in result.output
        assert "--output-ht" in result.output

    def test_help_shows_workflow(self):
        """Test that help shows workflow description."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "Workflow:" in result.output
        assert "PCA" in result.output
        assert "Random Forest" in result.output

    def test_help_shows_examples(self):
        """Test that help includes usage examples."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "Examples:" in result.output
        assert "hvantk ancestry-inference" in result.output

    def test_help_shows_options(self):
        """Test that help shows all major options."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        # Variant filtering options
        assert "--min-af" in result.output
        assert "--max-af" in result.output
        assert "--min-call-rate" in result.output

        # LD pruning options
        assert "--ld-r2" in result.output
        assert "--ld-window" in result.output
        assert "--skip-ld-pruning" in result.output

        # PCA options
        assert "--n-pcs" in result.output
        assert "--n-pcs-classify" in result.output

        # Classification options
        assert "--n-estimators" in result.output
        assert "--min-prob" in result.output

        # Output options
        assert "--generate-report" in result.output
        assert "--export-tsv" in result.output
        assert "--save-model" in result.output


class TestAncestryCLIValidation:
    """Test CLI argument validation."""

    def test_missing_required_args(self):
        """Test that missing required arguments shows error."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, [])

        assert result.exit_code != 0
        assert "Missing option" in result.output or "Error" in result.output

    def test_missing_query_mt(self):
        """Test error when --query-mt not provided."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            ref_path = Path(tmpdir) / "reference.mt"
            ref_path.mkdir()

            result = runner.invoke(
                ancestry_inference_cmd,
                [
                    "-r",
                    str(ref_path),
                    "-o",
                    str(Path(tmpdir) / "output.ht"),
                ],
            )

            assert result.exit_code != 0
            assert "Missing option" in result.output or "query" in result.output.lower()

    def test_missing_reference_mt(self):
        """Test error when --reference-mt not provided."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            query_path = Path(tmpdir) / "query.mt"
            query_path.mkdir()

            result = runner.invoke(
                ancestry_inference_cmd,
                [
                    "-q",
                    str(query_path),
                    "-o",
                    str(Path(tmpdir) / "output.ht"),
                ],
            )

            assert result.exit_code != 0
            assert "Missing option" in result.output or "reference" in result.output.lower()

    def test_missing_output_ht(self):
        """Test error when --output-ht not provided."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            query_path = Path(tmpdir) / "query.mt"
            ref_path = Path(tmpdir) / "reference.mt"
            query_path.mkdir()
            ref_path.mkdir()

            result = runner.invoke(
                ancestry_inference_cmd,
                [
                    "-q",
                    str(query_path),
                    "-r",
                    str(ref_path),
                ],
            )

            assert result.exit_code != 0
            assert "Missing option" in result.output or "output" in result.output.lower()

    def test_nonexistent_query_mt(self):
        """Test error when query MT doesn't exist."""
        runner = CliRunner()
        result = runner.invoke(
            ancestry_inference_cmd,
            [
                "-q",
                "/nonexistent/query.mt",
                "-r",
                "/fake/reference.mt",
                "-o",
                "/fake/output.ht",
            ],
        )

        assert result.exit_code != 0
        # Click validates that paths exist

    def test_nonexistent_reference_mt(self):
        """Test error when reference MT doesn't exist."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            query_path = Path(tmpdir) / "query.mt"
            query_path.mkdir()

            result = runner.invoke(
                ancestry_inference_cmd,
                [
                    "-q",
                    str(query_path),
                    "-r",
                    "/nonexistent/reference.mt",
                    "-o",
                    "/fake/output.ht",
                ],
            )

            assert result.exit_code != 0

    def test_invalid_log_level(self):
        """Test error with invalid log level."""
        runner = CliRunner()
        result = runner.invoke(
            ancestry_inference_cmd,
            [
                "-q",
                "/fake/query.mt",
                "-r",
                "/fake/reference.mt",
                "-o",
                "/fake/output.ht",
                "--log-level",
                "INVALID",
            ],
        )

        assert result.exit_code != 0


class TestAncestryCLIParameterValidation:
    """Test CLI parameter value validation.

    These tests require mocking Hail initialization and MatrixTable loading,
    which are done lazily inside the CLI command.
    """

    @patch("hvantk.core.hail_context.init_hail")
    @patch("hail.read_matrix_table")
    def test_invalid_min_af_negative(self, mock_read_mt, mock_init_hail):
        """Test error when min_af is negative."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            query_path = Path(tmpdir) / "query.mt"
            ref_path = Path(tmpdir) / "reference.mt"
            query_path.mkdir()
            ref_path.mkdir()

            result = runner.invoke(
                ancestry_inference_cmd,
                [
                    "-q",
                    str(query_path),
                    "-r",
                    str(ref_path),
                    "-o",
                    str(Path(tmpdir) / "output.ht"),
                    "--min-af",
                    "-0.1",
                ],
            )

            assert result.exit_code != 0
            assert "min-af" in result.output.lower() or "0" in result.output

    @patch("hvantk.core.hail_context.init_hail")
    @patch("hail.read_matrix_table")
    def test_invalid_min_af_too_high(self, mock_read_mt, mock_init_hail):
        """Test error when min_af is >= 0.5."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            query_path = Path(tmpdir) / "query.mt"
            ref_path = Path(tmpdir) / "reference.mt"
            query_path.mkdir()
            ref_path.mkdir()

            result = runner.invoke(
                ancestry_inference_cmd,
                [
                    "-q",
                    str(query_path),
                    "-r",
                    str(ref_path),
                    "-o",
                    str(Path(tmpdir) / "output.ht"),
                    "--min-af",
                    "0.6",
                ],
            )

            assert result.exit_code != 0

    @patch("hvantk.core.hail_context.init_hail")
    @patch("hail.read_matrix_table")
    def test_invalid_min_prob_out_of_range(self, mock_read_mt, mock_init_hail):
        """Test error when min_prob is out of range [0, 1]."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            query_path = Path(tmpdir) / "query.mt"
            ref_path = Path(tmpdir) / "reference.mt"
            query_path.mkdir()
            ref_path.mkdir()

            result = runner.invoke(
                ancestry_inference_cmd,
                [
                    "-q",
                    str(query_path),
                    "-r",
                    str(ref_path),
                    "-o",
                    str(Path(tmpdir) / "output.ht"),
                    "--min-prob",
                    "1.5",
                ],
            )

            assert result.exit_code != 0

    @patch("hvantk.core.hail_context.init_hail")
    @patch("hail.read_matrix_table")
    def test_n_pcs_classify_greater_than_n_pcs_warning(
        self, mock_read_mt, mock_init_hail
    ):
        """Test warning when n_pcs_classify > n_pcs."""
        runner = CliRunner()

        # Set up mocks to avoid actual pipeline execution
        mock_mt = MagicMock()
        mock_read_mt.return_value = mock_mt

        with tempfile.TemporaryDirectory() as tmpdir:
            query_path = Path(tmpdir) / "query.mt"
            ref_path = Path(tmpdir) / "reference.mt"
            query_path.mkdir()
            ref_path.mkdir()

            # This will fail later in pipeline but we just want to check the warning
            result = runner.invoke(
                ancestry_inference_cmd,
                [
                    "-q",
                    str(query_path),
                    "-r",
                    str(ref_path),
                    "-o",
                    str(Path(tmpdir) / "output.ht"),
                    "--n-pcs",
                    "5",
                    "--n-pcs-classify",
                    "10",
                ],
            )

            # Should show warning about n_pcs_classify
            assert "n-pcs-classify" in result.output.lower() or "5" in result.output


class TestAncestryCLIDefaults:
    """Test CLI default values."""

    def test_default_ancestry_col(self):
        """Test default value for --ancestry-col."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "[default: ancestry]" in result.output

    def test_default_min_af(self):
        """Test default value for --min-af."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "[default: 0.01]" in result.output

    def test_default_min_prob(self):
        """Test default value for --min-prob."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "[default: 0.75]" in result.output

    def test_default_n_pcs(self):
        """Test default value for --n-pcs."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "[default: 20]" in result.output

    def test_default_n_pcs_classify(self):
        """Test default value for --n-pcs-classify."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "[default: 10]" in result.output

    def test_default_ld_r2(self):
        """Test default value for --ld-r2."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "[default: 0.2]" in result.output

    def test_default_n_cv_folds(self):
        """Test default value for --n-cv-folds."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "[default: 5]" in result.output


class TestAncestryCLIFlags:
    """Test CLI flag options."""

    def test_skip_ld_pruning_flag_in_help(self):
        """Test that --skip-ld-pruning is documented."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "--skip-ld-pruning" in result.output

    def test_skip_validation_flag_in_help(self):
        """Test that --skip-validation is documented."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "--skip-validation" in result.output

    def test_generate_report_flag_in_help(self):
        """Test that --generate-report is documented."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "--generate-report" in result.output

    def test_export_tsv_flag_in_help(self):
        """Test that --export-tsv is documented."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "--export-tsv" in result.output

    def test_save_model_flag_in_help(self):
        """Test that --save-model is documented."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "--save-model" in result.output

    def test_save_loadings_flag_in_help(self):
        """Test that --save-loadings is documented."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "--save-loadings" in result.output

    def test_overwrite_flag_in_help(self):
        """Test that --overwrite is documented."""
        runner = CliRunner()
        result = runner.invoke(ancestry_inference_cmd, ["--help"])

        assert "--overwrite" in result.output


@pytest.mark.hail
class TestAncestryCLIIntegration:
    """Integration tests for ancestry CLI (require Hail)."""

    def test_basic_pipeline_execution(
        self, hail_session, synthetic_query_mt, synthetic_reference_mt
    ):
        """Test basic CLI execution with synthetic data."""
        import tempfile

        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            # Write synthetic MTs to disk
            query_path = str(Path(tmpdir) / "query.mt")
            ref_path = str(Path(tmpdir) / "reference.mt")
            output_path = str(Path(tmpdir) / "output.ht")

            synthetic_query_mt.write(query_path)
            synthetic_reference_mt.write(ref_path)

            result = runner.invoke(
                ancestry_inference_cmd,
                [
                    "-q",
                    query_path,
                    "-r",
                    ref_path,
                    "--ancestry-col",
                    "ancestry",
                    "-o",
                    output_path,
                    "--skip-ld-pruning",
                    "--skip-validation",
                    "--n-pcs",
                    "5",
                    "--n-pcs-classify",
                    "5",
                    "--min-shared-variants",
                    "100",
                ],
            )

            # Should succeed
            assert result.exit_code == 0, f"Failed with: {result.output}"
            assert "ANCESTRY INFERENCE COMPLETE" in result.output
            assert Path(output_path).exists()

    def test_pipeline_with_tsv_export(
        self, hail_session, synthetic_query_mt, synthetic_reference_mt
    ):
        """Test CLI with TSV export."""
        import tempfile

        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            query_path = str(Path(tmpdir) / "query.mt")
            ref_path = str(Path(tmpdir) / "reference.mt")
            output_path = str(Path(tmpdir) / "output.ht")
            output_dir = str(Path(tmpdir) / "results")

            synthetic_query_mt.write(query_path)
            synthetic_reference_mt.write(ref_path)

            result = runner.invoke(
                ancestry_inference_cmd,
                [
                    "-q",
                    query_path,
                    "-r",
                    ref_path,
                    "--ancestry-col",
                    "ancestry",
                    "-o",
                    output_path,
                    "--output-dir",
                    output_dir,
                    "--skip-ld-pruning",
                    "--skip-validation",
                    "--n-pcs",
                    "5",
                    "--n-pcs-classify",
                    "5",
                    "--min-shared-variants",
                    "100",
                    "--export-tsv",
                ],
            )

            assert result.exit_code == 0, f"Failed with: {result.output}"
            assert Path(output_dir, "predictions.tsv").exists()

    def test_pipeline_summary_output(
        self, hail_session, synthetic_query_mt, synthetic_reference_mt
    ):
        """Test CLI summary output content."""
        import tempfile

        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            query_path = str(Path(tmpdir) / "query.mt")
            ref_path = str(Path(tmpdir) / "reference.mt")
            output_path = str(Path(tmpdir) / "output.ht")

            synthetic_query_mt.write(query_path)
            synthetic_reference_mt.write(ref_path)

            result = runner.invoke(
                ancestry_inference_cmd,
                [
                    "-q",
                    query_path,
                    "-r",
                    ref_path,
                    "--ancestry-col",
                    "ancestry",
                    "-o",
                    output_path,
                    "--skip-ld-pruning",
                    "--skip-validation",
                    "--n-pcs",
                    "5",
                    "--n-pcs-classify",
                    "5",
                    "--min-shared-variants",
                    "100",
                ],
            )

            assert result.exit_code == 0

            # Check summary output
            assert "Query samples:" in result.output
            assert "Assigned:" in result.output
            assert "Unassigned:" in result.output
            assert "Ancestry distribution:" in result.output


if __name__ == "__main__":
    print("Running ancestry CLI tests...")

    print("\n1. Testing CLI help")
    test_help = TestAncestryCLIHelp()
    test_help.test_help_displays()
    print("  OK Help displays correctly")
    test_help.test_help_shows_workflow()
    print("  OK Workflow documented")
    test_help.test_help_shows_examples()
    print("  OK Examples included")
    test_help.test_help_shows_options()
    print("  OK All options shown")

    print("\n2. Testing CLI validation")
    test_validation = TestAncestryCLIValidation()
    test_validation.test_missing_required_args()
    print("  OK Missing required args detected")
    test_validation.test_missing_query_mt()
    print("  OK Missing query MT detected")
    test_validation.test_missing_reference_mt()
    print("  OK Missing reference MT detected")
    test_validation.test_missing_output_ht()
    print("  OK Missing output HT detected")

    print("\n3. Testing CLI defaults")
    test_defaults = TestAncestryCLIDefaults()
    test_defaults.test_default_ancestry_col()
    test_defaults.test_default_min_af()
    test_defaults.test_default_min_prob()
    test_defaults.test_default_n_pcs()
    print("  OK Defaults documented correctly")

    print("\n4. Testing CLI flags")
    test_flags = TestAncestryCLIFlags()
    test_flags.test_skip_ld_pruning_flag_in_help()
    test_flags.test_export_tsv_flag_in_help()
    test_flags.test_generate_report_flag_in_help()
    print("  OK All flags documented")

    print("\nAll basic tests passed!")
    print("Note: Integration tests require Hail and should be run with pytest -m hail")
