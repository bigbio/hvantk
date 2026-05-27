"""
Test script for HGC Pipeline functionality

This script tests the HGC pipeline module to ensure all components work correctly.
"""

import tempfile
import shutil
from pathlib import Path
import pytest

from hvantk.algorithms.hgc.pipeline import (
    PipelineConfig,
    PipelineState,
    PipelineRunner,
    PipelineStage,
)


class TestPipelineConfig:
    """Test PipelineConfig validation and creation."""

    def test_config_creation(self):
        """Test creating a basic configuration."""
        config = PipelineConfig(input_dir="/test/input", output_dir="/test/output")
        assert config.input_dir == "/test/input"
        assert config.output_dir == "/test/output"
        assert config.reference_genome == "GRCh38"
        assert config.overwrite is False

    def test_config_validation_missing_input(self):
        """Test validation catches missing input directory."""
        config = PipelineConfig(
            input_dir="/nonexistent/path", output_dir="/test/output"
        )
        errors = config.validate()
        assert len(errors) > 0
        assert any("does not exist" in err for err in errors)

    def test_config_validation_skip_without_path(self):
        """Test validation catches skip flags without required paths."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = PipelineConfig(
                input_dir=tmpdir,
                output_dir=tmpdir,
                skip_combine_gvcfs=True,
                # Missing vds_path
            )
            errors = config.validate()
            assert len(errors) > 0
            assert any("vds-path required" in err.lower() for err in errors)

    def test_config_validation_qc_thresholds(self):
        """Test validation of QC threshold values."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = PipelineConfig(
                input_dir=tmpdir,
                output_dir=tmpdir,
                min_sample_call_rate=1.5,  # Invalid: > 1.0
            )
            errors = config.validate()
            assert len(errors) > 0
            assert any("must be between 0 and 1" in err for err in errors)


class TestPipelineState:
    """Test PipelineState save/load functionality."""

    def test_state_creation(self):
        """Test creating a pipeline state."""
        config = PipelineConfig(input_dir="/test/input", output_dir="/test/output")
        state = PipelineState(config=config)
        assert state.config == config
        assert state.current_stage is None
        assert len(state.completed_stages) == 0
        assert len(state.errors) == 0

    def test_state_mark_complete(self):
        """Test marking stages as complete."""
        config = PipelineConfig(input_dir="/test/input", output_dir="/test/output")
        state = PipelineState(config=config)

        state.mark_stage_complete(
            PipelineStage.COMBINE_GVCFS, output_path="/test/output.vds"
        )

        assert PipelineStage.COMBINE_GVCFS.value in state.completed_stages
        assert state.outputs[PipelineStage.COMBINE_GVCFS.value] == "/test/output.vds"
        assert state.is_stage_complete(PipelineStage.COMBINE_GVCFS)
        assert not state.is_stage_complete(PipelineStage.VDS_TO_MT)

    def test_state_save_load(self):
        """Test saving and loading pipeline state."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = PipelineConfig(input_dir=tmpdir, output_dir=tmpdir)
            state = PipelineState(config=config)
            state.mark_stage_complete(
                PipelineStage.COMBINE_GVCFS, output_path="/test/output.vds"
            )
            state.start_time = "2026-01-14T10:00:00"

            # Save state
            state_path = Path(tmpdir) / "state.json"
            state.save(state_path)

            # Load state
            loaded_state = PipelineState.load(state_path)

            assert loaded_state.config.input_dir == tmpdir
            assert loaded_state.config.output_dir == tmpdir
            assert PipelineStage.COMBINE_GVCFS.value in loaded_state.completed_stages
            assert loaded_state.start_time == "2026-01-14T10:00:00"


@pytest.mark.hail
class TestPipelineRunner:
    """Test PipelineRunner orchestration (without actually running Hail)."""

    def test_runner_initialization(self):
        """Test runner initialization and path setup."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = PipelineConfig(
                input_dir=tmpdir, output_dir=tmpdir, output_prefix="test_cohort"
            )
            runner = PipelineRunner(config)

            assert runner.config == config
            assert "vds" in runner.paths
            assert "mt" in runner.paths
            assert "pvcf" in runner.paths
            assert "qc_report" in runner.paths
            assert runner.paths["vds"].endswith("test_cohort.vds")

    def test_runner_show_plan(self, capsys):
        """Test that show_plan displays execution plan."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = PipelineConfig(input_dir=tmpdir, output_dir=tmpdir)
            runner = PipelineRunner(config)
            runner.show_plan()

            captured = capsys.readouterr()
            assert "HGC PIPELINE EXECUTION PLAN" in captured.out
            assert "Combine gVCFs" in captured.out
            assert "Convert VDS" in captured.out
            assert "Compute Sample QC" in captured.out
            assert "Compute Variant QC" in captured.out
            assert "Export cohort VCF" in captured.out

    def test_runner_show_plan_with_skips(self, capsys):
        """Test show_plan with skip flags."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create dummy vds and mt paths
            vds_path = Path(tmpdir) / "test.vds"
            vds_path.mkdir()

            config = PipelineConfig(
                input_dir=tmpdir,
                output_dir=tmpdir,
                skip_combine_gvcfs=True,
                skip_compute_sample_qc=True,
                vds_path=str(vds_path),
            )
            runner = PipelineRunner(config)
            runner.show_plan()

            captured = capsys.readouterr()
            assert "[SKIP]" in captured.out
            assert "using existing VDS" in captured.out


def test_pipeline_stages_enum():
    """Test PipelineStage enum values."""
    assert PipelineStage.COMBINE_GVCFS.value == "combine_gvcfs"
    assert PipelineStage.VDS_TO_MT.value == "vds_to_mt"
    assert PipelineStage.COMPUTE_SAMPLE_QC.value == "compute_sample_qc"
    assert PipelineStage.COMPUTE_VARIANT_QC.value == "compute_variant_qc"
    assert PipelineStage.EXPORT_PVCF.value == "export_pvcf"


if __name__ == "__main__":
    print("Running HGC Pipeline tests...")
    print("\nTest 1: Configuration validation")
    test = TestPipelineConfig()
    test.test_config_creation()
    print("  ✓ Config creation works")

    test.test_config_validation_missing_input()
    print("  ✓ Validates missing input directory")

    test.test_config_validation_qc_thresholds()
    print("  ✓ Validates QC thresholds")

    print("\nTest 2: State management")
    test_state = TestPipelineState()
    test_state.test_state_creation()
    print("  ✓ State creation works")

    test_state.test_state_mark_complete()
    print("  ✓ Stage completion tracking works")

    test_state.test_state_save_load()
    print("  ✓ State save/load works")

    print("\nTest 3: Pipeline runner")
    test_runner = TestPipelineRunner()
    test_runner.test_runner_initialization()
    print("  ✓ Runner initialization works")

    print("\nTest 4: Enums")
    test_pipeline_stages_enum()
    print("  ✓ Pipeline stages enum works")

    print("\n✅ All tests passed!")
