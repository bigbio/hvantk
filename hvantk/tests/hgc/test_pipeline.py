"""
Test script for HGC Pipeline functionality

This script tests the HGC pipeline module to ensure all components work correctly.
"""

import json
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

    def test_state_load_ignores_removed_config_keys(self, tmp_path):
        """A state file from 0.6.0 still loads: it carries the removed keep_intermediates."""
        config = PipelineConfig(input_dir=str(tmp_path), output_dir=str(tmp_path))
        state_path = tmp_path / "state.json"
        PipelineState(config=config).save(state_path)
        saved = json.loads(state_path.read_text())
        saved["config"]["keep_intermediates"] = True
        state_path.write_text(json.dumps(saved))

        loaded_state = PipelineState.load(state_path)

        assert loaded_state.config.output_dir == str(tmp_path)


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


def test_densify_stage_reads_the_vds_the_combine_stage_wrote(monkeypatch, tmp_path):
    """A fresh run (no ``--vds-path``) must hand the combine stage's VDS to the densify
    stage. The stage runner recorded the output under the stage name only, while
    ``_run_vds_to_mt`` reads ``outputs["vds"]``, so ``convert_vds_to_mt`` received
    ``vds_path=None`` and the pipeline died at its second stage on every real input.
    """
    from unittest.mock import MagicMock

    from hvantk.algorithms.hgc import pipeline as pipeline_module

    monkeypatch.setattr(PipelineRunner, "_initialize_hail", lambda self: None)
    combine = MagicMock(name="combine_gvcfs")
    densify = MagicMock(name="convert_vds_to_mt")
    monkeypatch.setattr(pipeline_module, "combine_gvcfs", combine)
    monkeypatch.setattr(pipeline_module, "convert_vds_to_mt", densify)

    config = PipelineConfig(
        input_dir=str(tmp_path),
        output_dir=str(tmp_path / "out"),
        output_prefix="cohort",
    )
    runner = PipelineRunner(config)
    runner._run_stage(PipelineStage.COMBINE_GVCFS)
    runner._run_stage(PipelineStage.VDS_TO_MT)

    assert runner.state.outputs["vds"] == runner.paths["vds"]
    assert (
        runner.state.outputs[PipelineStage.COMBINE_GVCFS.value] == runner.paths["vds"]
    )
    assert densify.call_args.kwargs["vds_path"] == runner.paths["vds"]
    assert runner.state.outputs["mt"] == runner.paths["mt"]


def test_run_writes_the_pipeline_log(monkeypatch, tmp_path, caplog):
    """run() writes logs/pipeline_<timestamp>.log with the messages of every hgc module,
    and detaches the file handler when it ends."""
    import logging

    caplog.set_level(logging.INFO, logger="hvantk.algorithms.hgc")
    # No Hail start-up; the stand-in logs as the qc module does.
    monkeypatch.setattr(
        PipelineRunner,
        "_initialize_hail",
        lambda self: logging.getLogger("hvantk.algorithms.hgc.qc").info("a qc message"),
    )
    config = PipelineConfig(
        input_dir=str(tmp_path),
        output_dir=str(tmp_path / "out"),
        skip_combine_gvcfs=True,
        vds_path=str(tmp_path / "cohort.vds"),
        skip_vds_to_mt=True,
        mt_path=str(tmp_path / "cohort.mt"),
        skip_compute_sample_qc=True,
        skip_compute_variant_qc=True,
        skip_export_pvcf=True,
    )
    runner = PipelineRunner(config)
    runner.run()

    log_text = Path(runner.paths["log"]).read_text()
    assert "Starting HGC pipeline execution" in log_text
    assert "a qc message" in log_text
    assert not any(
        isinstance(h, logging.FileHandler)
        for h in logging.getLogger("hvantk.algorithms.hgc").handlers
    )


@pytest.mark.hail
def test_resumed_run_filters_with_the_earlier_runs_qc_tables(tmp_path):
    """--apply-qc-filters on a resume that skips both QC stages filters with the tables
    an earlier run left in the output directory, never exporting an unfiltered VCF, and
    stops before any stage when they are not there."""
    import gzip

    import hail as hl

    mt = hl.balding_nichols_model(
        n_populations=2, n_samples=20, n_variants=50, n_partitions=2
    )
    mt = mt.key_cols_by(s=hl.str(mt.sample_idx))
    # Sample "0" is never called and the first five variants are called in half the
    # samples: both fail the 0.85 call-rate thresholds.
    mt = mt.annotate_entries(
        GT=hl.or_missing(
            (mt.s != "0") & ~((mt.locus.position <= 5) & (mt.sample_idx % 2 == 0)),
            mt.GT,
        )
    )
    mt_path = str(tmp_path / "input.mt")
    mt.write(mt_path)
    dummy_vds = tmp_path / "unused.vds"
    dummy_vds.mkdir()

    def runner(**overrides):
        return PipelineRunner(
            PipelineConfig(
                input_dir=str(tmp_path),
                output_dir=str(tmp_path / "output"),
                skip_combine_gvcfs=True,
                vds_path=str(dummy_vds),
                skip_vds_to_mt=True,
                mt_path=mt_path,
                **overrides,
            )
        )

    resume = dict(
        skip_compute_sample_qc=True, skip_compute_variant_qc=True, apply_qc_filters=True
    )
    with pytest.raises(FileNotFoundError, match="cohort_sample_qc.ht"):
        runner(**resume).run()

    runner(skip_export_pvcf=True).run()
    resumed = runner(**resume)
    resumed.run()

    with gzip.open(resumed.paths["pvcf"], "rt") as fh:
        lines = [line for line in fh if not line.startswith("##")]
    samples = lines[0].rstrip("\n").split("\t")[9:]
    positions = {int(line.split("\t")[1]) for line in lines[1:]}
    assert len(samples) == 19 and "0" not in samples
    assert positions and not positions & {1, 2, 3, 4, 5}


@pytest.mark.hail
def test_pipeline_qc_filter_and_report_use_the_qc_tables_without_an_annotated_mt(
    tmp_path,
):
    """``--apply-qc-filters`` and ``--generate-qc-report`` work from the QC tables the
    compute stages wrote, and no QC-annotated MatrixTable copy is written (see
    ``PipelineRunner._read_qc_table``)."""
    import hail as hl

    from hvantk.algorithms.hgc.pipeline import PipelineConfig, PipelineRunner

    mt = hl.balding_nichols_model(
        n_populations=2, n_samples=10, n_variants=50, n_partitions=2
    )
    mt = mt.annotate_entries(
        GQ=hl.int32(hl.rand_unif(10, 40)), DP=hl.int32(hl.rand_unif(5, 30))
    )
    # export_vcf (stage 5) requires a single string column key; balding_nichols_model
    # keys columns by the integer sample_idx.
    mt = mt.key_cols_by(s=hl.str(mt.sample_idx))
    mt_path = str(tmp_path / "input.mt")
    mt.write(mt_path)

    # --skip-combine-gvcfs needs an existing --vds-path even though stage 1 never
    # runs: PipelineConfig.validate() only checks that the path exists.
    dummy_vds = tmp_path / "unused.vds"
    dummy_vds.mkdir()

    output_dir = tmp_path / "output"
    config = PipelineConfig(
        input_dir=str(tmp_path),
        output_dir=str(output_dir),
        skip_combine_gvcfs=True,
        vds_path=str(dummy_vds),
        skip_vds_to_mt=True,
        mt_path=mt_path,
        apply_qc_filters=True,
        generate_qc_report=True,
    )
    assert config.validate() == []

    runner = PipelineRunner(config)
    state = runner.run()

    assert state.errors == []
    # In the output directory, not next to --mt-path.
    assert (output_dir / "cohort_filtered.mt").exists()
    report_path = Path(runner.paths["qc_report"])
    assert report_path.exists()
    assert report_path.stat().st_size > 0

    # Each QC text file is delimited the way its extension says.
    qc_dir = output_dir / "qc"
    with open(qc_dir / "cohort_sample_qc.csv") as fh:
        sample_header = fh.readline()
    assert "," in sample_header
    assert "\t" not in sample_header

    with open(qc_dir / "cohort_variant_qc.tsv") as fh:
        variant_header = fh.readline()
    assert "\t" in variant_header

    # Not even with both QC consumers on.
    assert not (qc_dir / "cohort_with_qc.mt").exists()
