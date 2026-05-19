"""Tests for `hvantk reprocess <provider:dataset>` orchestrator."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest
from click.testing import CliRunner

from hvantk.tools.plugins.reprocess_cli import reprocess_cmd
from hvantk.core import plugin_loader
from hvantk.core.plugin_api import DatasetSpec, TestPaths


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
    download.assert_called_once_with(
        raw_dir=str(raw), cancer_type="brca", overwrite="true"
    )
    parse.assert_called_once_with(
        raw_dir=str(raw),
        output_path=str(intermediate),
        cancer_type="brca",
        overwrite="true",
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
