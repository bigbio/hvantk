"""
Tests for HGC CLI pipeline command.
"""

from unittest.mock import patch, MagicMock
from click.testing import CliRunner
import pytest

from hvantk.commands.hgc.pipeline_cli import pipeline


def test_pipeline_cli_basic():
    """Test pipeline command with basic options."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        import os
        os.makedirs("input_dir")
        os.makedirs("output_dir")

        with patch("hvantk.hgc.pipeline.PipelineConfig") as mock_config:
            with patch("hvantk.hgc.pipeline.PipelineRunner") as mock_runner:
                # Mock configuration
                mock_cfg_instance = MagicMock()
                mock_cfg_instance.validate.return_value = []
                mock_config.return_value = mock_cfg_instance

                # Mock runner
                mock_runner_instance = MagicMock()
                mock_state = MagicMock()
                mock_state.errors = []
                mock_state.outputs = {"vds": "/out/cohort.vds"}
                mock_state.completed_stages = ["combine_gvcfs"]
                mock_state.start_time = "00:00:00"
                mock_state.end_time = "00:10:00"
                mock_runner_instance.run.return_value = mock_state
                mock_runner.return_value = mock_runner_instance

                result = runner.invoke(
                    pipeline,
                    [
                        "--input-dir", "input_dir",
                        "--output-dir", "output_dir"
                    ],
                )

                assert result.exit_code == 0
                assert "Pipeline completed successfully" in result.output


def test_pipeline_cli_validation_failure():
    """Test pipeline command with validation failure."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        import os
        os.makedirs("input_dir")
        os.makedirs("output_dir")

        with patch("hvantk.hgc.pipeline.PipelineConfig") as mock_config:
            mock_cfg_instance = MagicMock()
            mock_cfg_instance.validate.return_value = ["Missing required parameter"]
            mock_config.return_value = mock_cfg_instance

            result = runner.invoke(
                pipeline,
                [
                    "--input-dir", "input_dir",
                    "--output-dir", "output_dir"
                ],
            )

            assert result.exit_code == 1
            assert "Configuration validation failed" in result.output
            assert "Missing required parameter" in result.output


def test_pipeline_cli_dry_run():
    """Test pipeline command with dry-run."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        import os
        os.makedirs("input_dir")
        os.makedirs("output_dir")

        with patch("hvantk.hgc.pipeline.PipelineConfig") as mock_config:
            with patch("hvantk.hgc.pipeline.PipelineRunner") as mock_runner:
                mock_cfg_instance = MagicMock()
                mock_cfg_instance.validate.return_value = []
                mock_config.return_value = mock_cfg_instance

                mock_runner_instance = MagicMock()
                mock_runner.return_value = mock_runner_instance

                result = runner.invoke(
                    pipeline,
                    [
                        "--input-dir", "input_dir",
                        "--output-dir", "output_dir",
                        "--dry-run"
                    ],
                )

                assert result.exit_code == 0
                mock_runner_instance.show_plan.assert_called_once()
                mock_runner_instance.run.assert_not_called()


def test_pipeline_cli_with_errors():
    """Test pipeline command when pipeline returns errors."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        import os
        os.makedirs("input_dir")
        os.makedirs("output_dir")

        with patch("hvantk.hgc.pipeline.PipelineConfig") as mock_config:
            with patch("hvantk.hgc.pipeline.PipelineRunner") as mock_runner:
                mock_cfg_instance = MagicMock()
                mock_cfg_instance.validate.return_value = []
                mock_config.return_value = mock_cfg_instance

                mock_runner_instance = MagicMock()
                mock_state = MagicMock()
                mock_state.errors = ["GVCF combination failed", "Unknown error"]
                mock_runner_instance.run.return_value = mock_state
                mock_runner.return_value = mock_runner_instance

                result = runner.invoke(
                    pipeline,
                    [
                        "--input-dir", "input_dir",
                        "--output-dir", "output_dir"
                    ],
                )

                assert result.exit_code == 1
                assert "Pipeline completed with errors" in result.output
                assert "GVCF combination failed" in result.output


def test_pipeline_cli_custom_options():
    """Test pipeline command with custom options."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        import os
        os.makedirs("input_dir")
        os.makedirs("output_dir")

        with patch("hvantk.hgc.pipeline.PipelineConfig") as mock_config:
            with patch("hvantk.hgc.pipeline.PipelineRunner") as mock_runner:
                mock_cfg_instance = MagicMock()
                mock_cfg_instance.validate.return_value = []
                mock_config.return_value = mock_cfg_instance

                mock_runner_instance = MagicMock()
                mock_state = MagicMock()
                mock_state.errors = []
                mock_state.outputs = {}
                mock_state.completed_stages = []
                mock_state.start_time = "00:00:00"
                mock_state.end_time = "00:10:00"
                mock_runner_instance.run.return_value = mock_state
                mock_runner.return_value = mock_runner_instance

                result = runner.invoke(
                    pipeline,
                    [
                        "--input-dir", "input_dir",
                        "--output-dir", "output_dir",
                        "--skip-compute-sample-qc",
                        "--reference-genome", "GRCh37",
                        "--output-prefix", "my_cohort"
                    ],
                )

                assert result.exit_code == 0
                # Check that custom options were passed to config
                call_kwargs = mock_config.call_args[1]
                assert call_kwargs["skip_compute_sample_qc"] is True
                assert call_kwargs["reference_genome"] == "GRCh37"
                assert call_kwargs["output_prefix"] == "my_cohort"


def test_pipeline_cli_skip_stages():
    """Test pipeline command with multiple skip flags."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        import os
        os.makedirs("input_dir")
        os.makedirs("output_dir")

        with patch("hvantk.hgc.pipeline.PipelineConfig") as mock_config:
            with patch("hvantk.hgc.pipeline.PipelineRunner") as mock_runner:
                mock_cfg_instance = MagicMock()
                mock_cfg_instance.validate.return_value = []
                mock_config.return_value = mock_cfg_instance

                mock_runner_instance = MagicMock()
                mock_state = MagicMock()
                mock_state.errors = []
                mock_state.outputs = {}
                mock_state.completed_stages = []
                mock_state.start_time = "00:00:00"
                mock_state.end_time = "00:10:00"
                mock_runner_instance.run.return_value = mock_state
                mock_runner.return_value = mock_runner_instance

                result = runner.invoke(
                    pipeline,
                    [
                        "--input-dir", "input_dir",
                        "--output-dir", "output_dir",
                        "--skip-combine-gvcfs",
                        "--skip-vds-to-mt",
                        "--vds-path", "/existing.vds",
                        "--mt-path", "/existing.mt"
                    ],
                )

                assert result.exit_code == 0
                call_kwargs = mock_config.call_args[1]
                assert call_kwargs["skip_combine_gvcfs"] is True
                assert call_kwargs["skip_vds_to_mt"] is True
                assert call_kwargs["vds_path"] == "/existing.vds"
                assert call_kwargs["mt_path"] == "/existing.mt"
