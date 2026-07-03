"""#198: reprocess must (a) fail fast when --output's extension can't hold the
plugin's backend (e.g. a pandas AnnotationTable written to .ht silently invokes
to_hail(), needing a JVM), and (b) not hard-error on a missing --intermediate when
the plugin declares a parse stage.
"""
from __future__ import annotations

import click
import pytest


def test_expected_extensions():
    from hvantk.tools.plugins.reprocess_cli import _expected_extensions

    assert _expected_extensions("pandas") == (".parquet",)
    assert _expected_extensions("anndata") == (".h5ad",)
    assert ".ht" in _expected_extensions("hail")
    assert _expected_extensions("unknown") is None


def test_check_output_extension_mismatch_raises():
    from hvantk.tools.plugins.reprocess_cli import _check_output_extension

    class _Spec:
        name = "peptideatlas:phospho"
        backend = "pandas"

    with pytest.raises(click.UsageError):
        _check_output_extension(_Spec(), "out.ht")
    # Matching extension must not raise.
    _check_output_extension(_Spec(), "out.parquet")


def test_check_output_extension_unknown_backend_is_noop():
    from hvantk.tools.plugins.reprocess_cli import _check_output_extension

    class _Spec:
        name = "x:y"
        backend = "something-else"

    # No expected extensions -> no opinion, no raise.
    _check_output_extension(_Spec(), "out.whatever")


def test_default_intermediate_path():
    from hvantk.tools.plugins.reprocess_cli import _default_intermediate

    p = _default_intermediate("peptideatlas:phospho", "/raw")
    assert p.endswith("peptideatlas_phospho.intermediate")
    assert p.startswith("/raw")


def test_reprocess_cmd_rejects_backend_extension_mismatch(tmp_path, monkeypatch):
    """The extension check fires inside reprocess_cmd, before any download."""
    from click.testing import CliRunner
    import hvantk.tools.plugins.reprocess_cli as rc

    class _Spec:
        name = "peptideatlas:phospho"
        backend = "pandas"

    monkeypatch.setattr(
        "hvantk.core.plugin.loader.get_registry",
        lambda: type("_Reg", (), {"get_dataset": lambda self, k: _Spec()})(),
    )

    result = CliRunner().invoke(
        rc.reprocess_cmd,
        [
            "peptideatlas:phospho",
            "--raw-dir", str(tmp_path / "raw"),
            "--output", str(tmp_path / "out.ht"),  # pandas backend -> .ht is wrong
        ],
    )
    assert result.exit_code != 0
    assert "backend" in result.output.lower()


def test_reprocess_cmd_defaults_missing_intermediate(tmp_path, monkeypatch):
    """A parse-declaring plugin with no --intermediate defaults it instead of erroring."""
    from click.testing import CliRunner
    import hvantk.tools.plugins.reprocess_cli as rc

    calls = {}

    def _fake_parse(*, raw_dir, output_path, **kw):
        calls["output_path"] = output_path
        with open(output_path, "w") as fh:
            fh.write("ok")

    class _Spec:
        name = "peptideatlas:phospho"
        backend = "pandas"
        download_fn = None
        parse_fn = staticmethod(_fake_parse)
        artifact_type = None  # legacy Phase-A build shape: builder(input, output)
        plugin_version = "0.1.0"

        def builder(self, parsed_path, output):
            calls["built_from"] = parsed_path

    monkeypatch.setattr(
        "hvantk.core.plugin.loader.get_registry",
        lambda: type("_Reg", (), {"get_dataset": lambda self, k: _Spec()})(),
    )
    monkeypatch.setattr(
        "hvantk.core.plugin.drift_runner.run_drift_check",
        lambda ds: type("_R", (), {"status": "clean", "diff": None})(),
    )

    result = CliRunner().invoke(
        rc.reprocess_cmd,
        [
            "peptideatlas:phospho",
            "--raw-dir", str(tmp_path / "raw"),
            "--output", str(tmp_path / "out.parquet"),  # matches pandas backend
            "--skip-download",
        ],
    )
    assert result.exit_code == 0, result.output
    assert calls["output_path"].endswith("peptideatlas_phospho.intermediate")
    assert calls["built_from"] == calls["output_path"]
