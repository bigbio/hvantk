"""Tests for `hvantk reprocess <provider:dataset>` orchestrator."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest
from click.testing import CliRunner

from hvantk.tools.plugins.reprocess_cli import reprocess_cmd
from hvantk.core.plugin import loader as plugin_loader
from hvantk.core.plugin.api import DatasetSpec, TestPaths


def _make_test_paths() -> TestPaths:
    return TestPaths(
        command="pytest",
        fixture="/x/fixture",
        schema_snapshot="/x/schema.json",
        row_snapshot="/x/rows.json",
        drift_fingerprint="/x/fp.json",
    )


def _make_spec(
    name: str = "stub:default",
    *,
    download_fn=None,
    parse_fn=None,
    builder=None,
) -> DatasetSpec:
    return DatasetSpec(
        name=name,
        domain="genomics",
        backend="hail",
        builder=builder or MagicMock(),
        drift_probe=MagicMock(return_value={"probe_version": 1}),
        skill_path="/x/SKILL.md",
        test_paths=_make_test_paths(),
        download_fn=download_fn,
        parse_fn=parse_fn,
    )


@pytest.fixture(autouse=True)
def isolated_registry():
    """Each test starts with a fresh, empty registry singleton."""
    plugin_loader.reset_registry_for_tests()
    yield
    plugin_loader.reset_registry_for_tests()


def _install_registry(monkeypatch, spec: DatasetSpec) -> plugin_loader.PluginRegistry:
    reg = plugin_loader.PluginRegistry()
    reg._providers["stub"] = MagicMock(name=spec.name.split(":")[0])
    reg._datasets[spec.name] = spec
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)
    return reg


def test_reprocess_full_pipeline_calls_all_stages(tmp_path: Path, monkeypatch):
    download = MagicMock()
    parse = MagicMock()
    builder = MagicMock()
    spec = _make_spec(download_fn=download, parse_fn=parse, builder=builder)
    _install_registry(monkeypatch, spec)

    raw = tmp_path / "raw"
    intermediate = tmp_path / "mid.parquet"
    output = tmp_path / "out.parquet"

    runner = CliRunner()
    result = runner.invoke(
        reprocess_cmd,
        [
            "stub:default",
            "--raw-dir",
            str(raw),
            "--intermediate",
            str(intermediate),
            "--output",
            str(output),
            "--no-check-drift",
        ],
    )
    assert result.exit_code == 0, result.output

    download.assert_called_once_with(raw_dir=str(raw))
    parse.assert_called_once_with(raw_dir=str(raw), output_path=str(intermediate))
    builder.assert_called_once_with(str(intermediate), str(output))


def test_reprocess_skip_download(tmp_path: Path, monkeypatch):
    download = MagicMock()
    parse = MagicMock()
    builder = MagicMock()
    spec = _make_spec(download_fn=download, parse_fn=parse, builder=builder)
    _install_registry(monkeypatch, spec)

    raw = tmp_path / "raw"
    intermediate = tmp_path / "mid.parquet"
    output = tmp_path / "out.parquet"

    runner = CliRunner()
    result = runner.invoke(
        reprocess_cmd,
        [
            "stub:default",
            "--raw-dir",
            str(raw),
            "--intermediate",
            str(intermediate),
            "--output",
            str(output),
            "--skip-download",
            "--no-check-drift",
        ],
    )
    assert result.exit_code == 0, result.output

    download.assert_not_called()
    parse.assert_called_once_with(raw_dir=str(raw), output_path=str(intermediate))
    builder.assert_called_once_with(str(intermediate), str(output))


def test_reprocess_no_lifecycle_declared_requires_skip_download(
    tmp_path: Path, monkeypatch
):
    builder = MagicMock()
    # No download_fn / parse_fn declared.
    spec = _make_spec(builder=builder)
    _install_registry(monkeypatch, spec)

    raw = tmp_path / "raw"
    output = tmp_path / "out.parquet"

    runner = CliRunner()
    result = runner.invoke(
        reprocess_cmd,
        [
            "stub:default",
            "--raw-dir",
            str(raw),
            "--output",
            str(output),
            "--no-check-drift",
        ],
    )
    # UsageError -> Click exit code 2.
    assert result.exit_code != 0
    assert "no lifecycle.download declared" in result.output
    builder.assert_not_called()


def test_reprocess_skip_download_no_parse_passes_raw_dir_to_builder(
    tmp_path: Path, monkeypatch
):
    """When no parse_fn is declared, the builder consumes raw_dir directly."""
    builder = MagicMock()
    spec = _make_spec(builder=builder)
    _install_registry(monkeypatch, spec)

    raw = tmp_path / "raw"
    output = tmp_path / "out.parquet"

    runner = CliRunner()
    result = runner.invoke(
        reprocess_cmd,
        [
            "stub:default",
            "--raw-dir",
            str(raw),
            "--output",
            str(output),
            "--skip-download",
            "--no-check-drift",
        ],
    )
    assert result.exit_code == 0, result.output
    builder.assert_called_once_with(str(raw), str(output))


def test_reprocess_unknown_dataset_errors(tmp_path: Path, monkeypatch):
    reg = plugin_loader.PluginRegistry()
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)

    runner = CliRunner()
    result = runner.invoke(
        reprocess_cmd,
        [
            "missing:dataset",
            "--raw-dir",
            str(tmp_path / "raw"),
            "--output",
            str(tmp_path / "out.parquet"),
            "--no-check-drift",
        ],
    )
    assert result.exit_code != 0
    assert "unknown dataset" in result.output.lower()


def test_reprocess_parse_requires_intermediate_path(tmp_path: Path, monkeypatch):
    """When a plugin declares lifecycle.parse, --intermediate must be supplied."""
    download = MagicMock()
    parse = MagicMock()
    builder = MagicMock()
    spec = _make_spec(download_fn=download, parse_fn=parse, builder=builder)
    _install_registry(monkeypatch, spec)

    runner = CliRunner()
    result = runner.invoke(
        reprocess_cmd,
        [
            "stub:default",
            "--raw-dir",
            str(tmp_path / "raw"),
            "--output",
            str(tmp_path / "out.parquet"),
            "--no-check-drift",
        ],
    )
    assert result.exit_code != 0
    assert "--intermediate" in result.output
    parse.assert_not_called()
    builder.assert_not_called()


def test_reprocess_skip_parse_without_intermediate_falls_back_to_raw_dir(
    tmp_path: Path, monkeypatch
):
    """`--skip-parse` without `--intermediate` must feed raw_dir to the builder.

    Regression test for the bug where parsed_path stayed None and crashed the
    builder with a confusing TypeError.
    """
    download = MagicMock()
    parse = MagicMock()
    builder = MagicMock()
    spec = _make_spec(download_fn=download, parse_fn=parse, builder=builder)
    _install_registry(monkeypatch, spec)

    raw = tmp_path / "raw"
    output = tmp_path / "out.parquet"

    runner = CliRunner()
    result = runner.invoke(
        reprocess_cmd,
        [
            "stub:default",
            "--raw-dir",
            str(raw),
            "--output",
            str(output),
            "--skip-parse",
            "--no-check-drift",
        ],
    )
    assert result.exit_code == 0, result.output
    download.assert_called_once_with(raw_dir=str(raw))
    parse.assert_not_called()
    # Builder consumes the raw dir because parse was skipped and intermediate
    # was not provided.
    builder.assert_called_once_with(str(raw), str(output))


def test_reprocess_skip_parse_with_intermediate_uses_intermediate(
    tmp_path: Path, monkeypatch
):
    """`--skip-parse --intermediate <path>` feeds the existing intermediate to build."""
    download = MagicMock()
    parse = MagicMock()
    builder = MagicMock()
    spec = _make_spec(download_fn=download, parse_fn=parse, builder=builder)
    _install_registry(monkeypatch, spec)

    raw = tmp_path / "raw"
    intermediate = tmp_path / "mid.tsv"
    output = tmp_path / "out.parquet"

    runner = CliRunner()
    result = runner.invoke(
        reprocess_cmd,
        [
            "stub:default",
            "--raw-dir",
            str(raw),
            "--intermediate",
            str(intermediate),
            "--output",
            str(output),
            "--skip-parse",
            "--no-check-drift",
        ],
    )
    assert result.exit_code == 0, result.output
    parse.assert_not_called()
    builder.assert_called_once_with(str(intermediate), str(output))


def test_reprocess_plugin_arg_forwarded_to_download_and_parse(
    tmp_path: Path, monkeypatch
):
    """`--plugin-arg KEY=VALUE` is forwarded as kwargs to download_fn and parse_fn.

    Regression test for the silent default-to-all behavior in cptac:phospho's
    download_dataset (cancer_type=None expands to ALL cancer types).
    """
    download = MagicMock()
    parse = MagicMock()
    builder = MagicMock()
    spec = _make_spec(download_fn=download, parse_fn=parse, builder=builder)
    _install_registry(monkeypatch, spec)

    raw = tmp_path / "raw"
    intermediate = tmp_path / "mid.tsv"
    output = tmp_path / "out.parquet"

    runner = CliRunner()
    result = runner.invoke(
        reprocess_cmd,
        [
            "stub:default",
            "--raw-dir",
            str(raw),
            "--intermediate",
            str(intermediate),
            "--output",
            str(output),
            "--plugin-arg",
            "cancer_type=brca",
            "--plugin-arg",
            "overwrite=true",
            "--no-check-drift",
        ],
    )
    assert result.exit_code == 0, result.output
    # --plugin-arg values are coerced: 'true'/'false' -> bool, numeric strings
    # -> int/float. 'brca' stays str because it matches no other type.
    download.assert_called_once_with(
        raw_dir=str(raw), cancer_type="brca", overwrite=True
    )
    parse.assert_called_once_with(
        raw_dir=str(raw),
        output_path=str(intermediate),
        cancer_type="brca",
        overwrite=True,
    )


def test_reprocess_plugin_arg_bad_format_errors(tmp_path: Path, monkeypatch):
    """`--plugin-arg` without '=' is rejected as a usage error."""
    download = MagicMock()
    builder = MagicMock()
    spec = _make_spec(download_fn=download, builder=builder)
    _install_registry(monkeypatch, spec)

    runner = CliRunner()
    result = runner.invoke(
        reprocess_cmd,
        [
            "stub:default",
            "--raw-dir",
            str(tmp_path / "raw"),
            "--output",
            str(tmp_path / "out.parquet"),
            "--plugin-arg",
            "noequals",
            "--no-check-drift",
        ],
    )
    assert result.exit_code != 0
    assert "KEY=VALUE" in result.output
    download.assert_not_called()
    builder.assert_not_called()


def test_reprocess_phase_b_plugin_uses_run_builder_for_spec(
    tmp_path: Path, monkeypatch
):
    """A spec with artifact_type set should route through run_builder_for_spec
    (which uses BuildContext), not the legacy spec.builder(input, output) shape.

    Also guards that --plugin-arg values reach the Phase B builder with type
    coercion applied — regressing this would silently drop build-time params
    (reference_genome, p_threshold, tissue, ...) and break the deep-dive docs.
    """
    import pandas as pd
    from hvantk.core.models import AnnotationTable, BuildContext

    # Real Phase B-style builder accepting (parsed, ctx, **params)
    captured_ctx: dict = {}
    captured_params: dict = {}

    def phase_b_builder(parsed, ctx, **params):
        captured_ctx["ctx"] = ctx
        captured_params.update(params)
        df = pd.DataFrame({"x": [1, 2, 3]})
        return AnnotationTable.from_pandas(
            df, provenance=ctx.provenance(schema_id="test-rows-v1")
        )

    spec = _make_spec(builder=phase_b_builder)
    # Set the Phase B fields on the spec — _make_spec doesn't set them by default
    object.__setattr__(spec, "artifact_type", AnnotationTable)
    object.__setattr__(spec, "schema_id", "test-rows-v1")
    object.__setattr__(spec, "plugin_version", "0.1.0")
    object.__setattr__(spec, "drift_probe", lambda: {"source_version": "v1"})
    _install_registry(monkeypatch, spec)

    raw = tmp_path / "raw"
    raw.mkdir()
    output = tmp_path / "out.parquet"

    result = CliRunner().invoke(
        reprocess_cmd,
        [
            "stub:default",
            "--raw-dir",
            str(raw),
            "--output",
            str(output),
            "--skip-download",
            "--plugin-arg",
            "reference_genome=GRCh38",
            "--plugin-arg",
            "p_threshold=5e-8",
            "--plugin-arg",
            "overwrite=false",
            "--no-check-drift",
        ],
    )
    assert result.exit_code == 0, result.output

    # Confirm the builder got a real BuildContext (not a string)
    assert isinstance(captured_ctx["ctx"], BuildContext)
    assert captured_ctx["ctx"].source_fingerprint  # nonempty
    # Confirm --plugin-arg values reached the Phase B builder with coercion
    assert captured_params["reference_genome"] == "GRCh38"
    assert captured_params["p_threshold"] == 5e-8
    assert captured_params["overwrite"] is False
    # Confirm the artifact was saved (file exists + sidecar manifest exists)
    assert output.exists()
    assert output.with_name(output.name + ".provenance.json").exists()
