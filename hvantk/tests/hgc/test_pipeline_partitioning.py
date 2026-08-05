"""The HGC pipeline must be able to tune gVCF-combiner partitioning.

`hvantk hgc pipeline` is the documented, recommended entry point, but it previously
passed ``kwargs={}`` to ``combine_gvcfs`` unconditionally -- so a pipeline user was stuck
with Hail's 1.2 Mb genome default and had no way to raise the partition count above their
core count. These tests pin the wiring so that cannot silently regress.
"""

import re

import pytest

from hvantk.algorithms.hgc.pipeline import PipelineConfig


def _cfg(tmp_path, **overrides):
    return PipelineConfig(
        input_dir=str(tmp_path), output_dir=str(tmp_path / "out"), **overrides
    )


def test_combiner_kwargs_empty_by_default(tmp_path):
    """No tuning options -> empty kwargs -> Hail's genome default (legacy behaviour)."""
    assert _cfg(tmp_path).combiner_kwargs() == {}


def test_combiner_kwargs_carries_import_interval_size(tmp_path):
    cfg = _cfg(tmp_path, import_interval_size=600_000)
    assert cfg.combiner_kwargs() == {"import_interval_size": 600_000}


def test_combiner_kwargs_carries_tree_merge_options(tmp_path):
    cfg = _cfg(tmp_path, gvcf_batch_size=100, branch_factor=50)
    assert cfg.combiner_kwargs() == {"gvcf_batch_size": 100, "branch_factor": 50}


def test_combiner_kwargs_carries_exome_default(tmp_path):
    cfg = _cfg(tmp_path, use_exome_default_intervals=True)
    assert cfg.combiner_kwargs() == {"use_exome_default_intervals": True}


@pytest.mark.parametrize(
    "overrides, expected",
    [
        ({"import_interval_size": 0}, "at least 1 bp"),
        ({"import_interval_size": -5}, "at least 1 bp"),
        ({"gvcf_batch_size": 0}, "at least 1"),
        ({"branch_factor": 1}, "at least 2"),
        (
            {"import_interval_size": 600_000, "use_exome_default_intervals": True},
            "mutually exclusive",
        ),
    ],
)
def test_invalid_combiner_tuning_is_rejected_by_validate(tmp_path, overrides, expected):
    """Bad values must be caught by config validation, not deep inside Hail.

    Hail does reject these, but only after init + gVCF validation, and the error is then
    wrapped in generic "check your Spark version / GVCF files" guidance that points the
    user at the wrong cause entirely.
    """
    errors = _cfg(tmp_path, **overrides).validate()
    assert any(expected in e for e in errors), errors


def test_valid_combiner_tuning_passes_validate(tmp_path):
    errors = _cfg(
        tmp_path, import_interval_size=600_000, gvcf_batch_size=50, branch_factor=100
    ).validate()
    combiner_errors = [
        e for e in errors if "interval" in e or "batch" in e or "branch" in e
    ]
    assert combiner_errors == []


# --- #208: n_partitions reached no pipeline stage ----------------------------------
#
# The field was accepted from the CLI and echoed back in the run plan while being read by
# nothing, so the run plan affirmatively told the user a setting had taken effect when it
# had not. These pin the wiring; the end-to-end partition-count assertion lives in
# test_hgc_core_hail.py, which needs Hail.


def _runner_without_hail(pipeline_mod, monkeypatch, cfg):
    """Build a PipelineRunner without booting a JVM.

    PipelineRunner.__init__ calls _initialize_hail(). None of these tests need Hail --
    show_plan() only prints strings and _run_vds_to_mt's converter is stubbed -- but
    without this the file boots Spark inside the DEFAULT `pytest -q` run, which
    pytest.ini deliberately configures to deselect the `hail` marker, and fails outright
    on any machine with no JVM. Round-3 review caught it: the pre-PR file produced zero
    Spark banners, this one produced three.
    """
    monkeypatch.setattr(pipeline_mod.PipelineRunner, "_initialize_hail", lambda self: None)
    return pipeline_mod.PipelineRunner(cfg)


def test_pipeline_forwards_n_partitions_to_the_converter(tmp_path, monkeypatch):
    """The regression: PipelineConfig.n_partitions must REACH convert_vds_to_mt."""
    from hvantk.algorithms.hgc import pipeline as pipeline_mod

    seen = {}

    def _spy(**kwargs):
        seen.update(kwargs)

    monkeypatch.setattr(pipeline_mod, "convert_vds_to_mt", _spy)

    runner = _runner_without_hail(pipeline_mod, monkeypatch, _cfg(tmp_path, n_partitions=64))
    runner.state.outputs["vds"] = str(tmp_path / "cohort.vds")
    runner._run_vds_to_mt()

    assert seen["n_partitions"] == 64


def test_pipeline_passes_none_when_unset(tmp_path, monkeypatch):
    """Unset must stay None rather than becoming a number, so the VDS layout is kept."""
    from hvantk.algorithms.hgc import pipeline as pipeline_mod

    seen = {}
    monkeypatch.setattr(pipeline_mod, "convert_vds_to_mt", lambda **kw: seen.update(kw))

    runner = _runner_without_hail(pipeline_mod, monkeypatch, _cfg(tmp_path))
    runner.state.outputs["vds"] = str(tmp_path / "cohort.vds")
    runner._run_vds_to_mt()

    assert seen["n_partitions"] is None


# --- adversarial-review findings on #261 --------------------------------------------


def test_invalid_n_partitions_is_rejected_by_validate(tmp_path):
    """Bad values must be caught by config validation, not deep inside stage 2.

    n_partitions was the only tuning knob with no validate() clause, so `--n-partitions 0`
    passed the gate, ran the multi-hour gVCF combine, and only then hit the guard inside
    convert_vds_to_mt. Its three siblings are all checked here for exactly this reason.
    """
    for bad in (0, -5):
        errors = _cfg(tmp_path, n_partitions=bad).validate()
        assert any("n_partitions" in e for e in errors), (bad, errors)


def test_valid_n_partitions_passes_validate(tmp_path):
    assert [e for e in _cfg(tmp_path, n_partitions=64).validate() if "n_partitions" in e] == []
    assert [e for e in _cfg(tmp_path).validate() if "n_partitions" in e] == []


def test_run_plan_does_not_render_zero_as_auto():
    """`or` printed 0 -- a value that aborts the run -- as "auto"."""
    from hvantk.algorithms.hgc.pipeline import _shown_partitions

    assert _shown_partitions(None) == "auto (VDS layout)"
    assert _shown_partitions(0) == "0"
    assert _shown_partitions(64) == "64"


def test_pipeline_help_names_flags_that_exist():
    """The --n-partitions help pointed at a `--combiner-*` family that does not exist."""
    from hvantk.tools.hgc.pipeline_cli import pipeline as pipeline_cmd

    opt = next(p for p in pipeline_cmd.params if "--n-partitions" in getattr(p, "opts", []))
    declared = {o for p in pipeline_cmd.params for o in getattr(p, "opts", [])}
    referenced = re.findall(r"--[a-z][a-z0-9-]+", opt.help)

    missing = [f for f in referenced if f not in declared]
    assert not missing, f"help references flags this command does not define: {missing}"


def test_vds2mt_dry_run_shows_zero_rather_than_auto(tmp_path):
    """The CALL SITE, not just the helper.

    Round-2 review: the helper test alone left both call sites unprotected -- reverting
    convert_cli.py to `n_partitions or 'auto (VDS layout)'` kept the suite green. This is
    the live half: `vds2mt` has no config.validate() gate, so 0 reaches the dry-run
    printer and the user is told a plan is fine for an invocation that aborts.
    """
    from unittest.mock import patch
    from click.testing import CliRunner

    from hvantk.tools.hgc.convert_cli import vds2mt

    with patch("hvantk.tools.hgc.convert_cli.validate_input_files", return_value=(True, [])), \
         patch("hvantk.tools.hgc.convert_cli.validate_output_path", return_value=True):
        result = CliRunner().invoke(
            vds2mt,
            ["-i", str(tmp_path / "in.vds"), "-o", str(tmp_path / "out.mt"),
             "--n-partitions", "0", "--dry-run"],
        )

    assert result.exit_code == 0, result.output
    # The Partitions LINE only -- result.output embeds tmp_path, so a blanket
    # `assert "auto" not in result.output` depends on pytest's tmp-dir naming rather
    # than on the code under test.
    line = next(ln for ln in result.output.splitlines() if "Partitions:" in ln)
    assert line.split(":", 1)[1].strip() == "0", line


def test_run_plan_call_site_shows_zero_rather_than_auto(tmp_path, monkeypatch):
    """The other call site: PipelineRunner.show_plan."""
    import contextlib
    import io

    from hvantk.algorithms.hgc import pipeline as pipeline_mod

    runner = _runner_without_hail(
        pipeline_mod, monkeypatch, _cfg(tmp_path, n_partitions=0)
    )
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        runner.show_plan()

    # Assert on the Partitions LINE, not the whole capture. A blanket
    # `assert "auto" not in out` scans output that embeds tmp_path, whose directory name
    # pytest derives from this function's own name truncated to 30 chars -- it passed
    # only because the cut landed one character before the trailing "auto". Renaming the
    # test would have broken it for reasons having nothing to do with the code.
    line = next(ln for ln in buf.getvalue().splitlines() if "Partitions (MT):" in ln)
    assert line.split(":", 1)[1].strip() == "0", line


def test_run_plan_does_not_advertise_partitions_when_the_stage_is_skipped(
    tmp_path, monkeypatch
):
    """--skip-vds-to-mt means convert_vds_to_mt never runs, so the existing MatrixTable
    keeps its layout. Printing the requested number there would repeat #208's mistake in
    a new place: affirming a setting no stage will read."""
    import contextlib
    import io

    from hvantk.algorithms.hgc import pipeline as pipeline_mod

    cfg = _cfg(tmp_path, n_partitions=64, skip_vds_to_mt=True, mt_path=str(tmp_path))
    runner = _runner_without_hail(pipeline_mod, monkeypatch, cfg)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        runner.show_plan()

    line = next(ln for ln in buf.getvalue().splitlines() if "Partitions (MT):" in ln)
    assert "64" not in line, line
    assert "skipped" in line.lower(), line
