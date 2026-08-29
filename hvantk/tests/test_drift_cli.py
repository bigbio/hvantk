"""Tests for `hvantk drift ...` Click commands."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from hvantk.core.plugin import loader as plugin_loader
from hvantk.tools.plugins.drift_cli import drift_cmd


FIXTURE_ROOT = Path(__file__).parent / "testdata" / "raw" / "plugins"


@pytest.fixture(autouse=True)
def reset_registry(monkeypatch):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)


def test_drift_clean_exit_zero():
    runner = CliRunner()
    result = runner.invoke(drift_cmd, ["fake:default"])
    assert result.exit_code == 0
    assert "clean" in result.output.lower()


def test_drift_all_runs_every_dataset():
    runner = CliRunner()
    result = runner.invoke(drift_cmd, ["--all"])
    assert result.exit_code == 0
    assert "fake:default" in result.output


def test_drift_json_output_is_parseable():
    runner = CliRunner()
    result = runner.invoke(drift_cmd, ["--json", "fake:default"])
    assert result.exit_code == 0
    parsed = json.loads(result.output)
    # --json with a single dataset returns a list with one entry (matches --all behavior).
    assert isinstance(parsed, list)
    assert parsed[0]["dataset_name"] == "fake:default"
    assert parsed[0]["status"] == "clean"


def test_drift_stub_emits_warning_to_stderr_in_json_mode():
    """Regression (PR #186 review): a stub probe must surface a WARNING on
    stderr even in --json mode, so the scheduled drift workflow — which captures
    `drift --all --json` stdout to a file — stays visibly non-green for
    documentation-only sources. stdout must remain clean machine-readable JSON.
    """
    from hvantk.core.plugin.api import stub_fingerprint

    spec = plugin_loader.get_registry().get_dataset("fake:default")
    object.__setattr__(
        spec, "drift_probe", lambda: stub_fingerprint("doc-only; no probeable URL")
    )

    # Click >= 8.2 removed the `mix_stderr` kwarg and always captures stdout
    # and stderr separately; Click 8.1.x needs `mix_stderr=False` to do so.
    # Support both so the test runs on either version.
    try:
        runner = CliRunner(mix_stderr=False)
    except TypeError:
        runner = CliRunner()
    result = runner.invoke(drift_cmd, ["--json", "fake:default"])

    assert result.exit_code == 0
    # `result.stdout` is stdout-only on both Click 8.1.x (via mix_stderr=False)
    # and 8.2+ (always separated); `result.output` mixes stderr in on 8.2+.
    parsed = json.loads(result.stdout)
    assert parsed[0]["status"] == "stub"
    assert "WARNING" not in result.stdout
    # the WARNING is on stderr so it shows up in CI step logs
    assert "WARNING" in result.stderr
    assert "stub probe" in result.stderr
    assert "doc-only; no probeable URL" in result.stderr


def test_drift_regenerate_overwrites_fingerprint(tmp_path: Path, monkeypatch):
    # Point the fixture at a tmpdir-copy so we don't mutate the test asset.
    import shutil
    plugin_dir = tmp_path / "fake_plugin"
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", plugin_dir)
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(plugin_dir)
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)
    fp_path = plugin_dir / "tests" / "drift_fingerprint.json"
    old = json.loads(fp_path.read_text())
    # Mutate the expected file so a regenerate visibly changes it.
    fp_path.write_text(json.dumps({"probe_version": 1, "stale": True}))
    runner = CliRunner()
    result = runner.invoke(drift_cmd, ["--regenerate", "fake:default"])
    assert result.exit_code == 0
    new = json.loads(fp_path.read_text())
    assert "stale" not in new
    assert new["probe_version"] == old["probe_version"]


def test_unknown_dataset_returns_registry_error_exit_code():
    runner = CliRunner()
    result = runner.invoke(drift_cmd, ["does:not:exist"])
    assert result.exit_code == 3  # EXIT_REGISTRY_ERROR
    assert "unknown dataset" in result.output.lower() or "unknown dataset" in (result.stderr or "")


def test_regenerate_unknown_dataset_returns_registry_error_exit_code():
    runner = CliRunner()
    result = runner.invoke(drift_cmd, ["--regenerate", "does:not:exist"])
    assert result.exit_code == 3


# --- --ledger ------------------------------------------------------------------------
#
# A fingerprint bump accepted into a PR is also the signal that a built artifact may
# now be stale. `hvantk drift --ledger` reads the rebuild ledger (written by the drift
# bot; see .github/scripts/drift_to_pr.py) and lists what is still pending a rebuild.

def test_ledger_flag_lists_datasets_needing_rebuild(tmp_path, monkeypatch):
    """A dataset whose upstream moved after its last rebuild is stale. Never-rebuilt
    (rebuilt_at None) counts as stale."""
    from hvantk.tools.plugins import drift_cli

    ledger = tmp_path / "drift_ledger.json"
    ledger.write_text(json.dumps({
        "clinvar:variants": {"last_upstream_change": "2026-08-23T00:00:00+00:00",
                             "accepted_in": "PR #288", "signal": "routine",
                             "rebuilt_at": None},
        "hgnc:lookup": {"last_upstream_change": "2026-08-01T00:00:00+00:00",
                        "accepted_in": "PR #286", "signal": "routine",
                        "rebuilt_at": "2026-08-20T00:00:00+00:00"},
    }))
    monkeypatch.setattr(drift_cli, "LEDGER_PATH", ledger)

    result = CliRunner().invoke(drift_cli.drift_cmd, ["--ledger"])

    assert result.exit_code == 0
    assert "clinvar:variants" in result.output
    assert "hgnc:lookup" not in result.output


def test_ledger_flag_reports_nothing_pending_on_empty_ledger(tmp_path, monkeypatch):
    from hvantk.tools.plugins import drift_cli

    ledger = tmp_path / "drift_ledger.json"
    ledger.write_text("{}")
    monkeypatch.setattr(drift_cli, "LEDGER_PATH", ledger)

    result = CliRunner().invoke(drift_cli.drift_cmd, ["--ledger"])
    assert result.exit_code == 0
    assert "no datasets pending rebuild" in result.output


def test_load_ledger_returns_empty_dict_on_truthy_non_dict_json(tmp_path, monkeypatch):
    """`json.loads(text) or {}` only substitutes `{}` for FALSY JSON -- a populated
    list or a bare string is truthy and passes straight through, breaking the
    docstring's "a missing or corrupt file yields {}" promise."""
    from hvantk.tools.plugins import drift_cli

    ledger = tmp_path / "drift_ledger.json"
    monkeypatch.setattr(drift_cli, "LEDGER_PATH", ledger)

    ledger.write_text("[1, 2, 3]")
    assert drift_cli._load_ledger() == {}

    ledger.write_text('"x"')
    assert drift_cli._load_ledger() == {}


def test_ledger_flag_survives_a_truthy_non_dict_ledger_file(tmp_path, monkeypatch):
    """`drift --ledger` against a ledger file containing a JSON list must exit 0 with
    "no datasets pending rebuild" rather than crash calling `.items()` on a list."""
    from hvantk.tools.plugins import drift_cli

    ledger = tmp_path / "drift_ledger.json"
    ledger.write_text("[1, 2, 3]")
    monkeypatch.setattr(drift_cli, "LEDGER_PATH", ledger)

    result = CliRunner().invoke(drift_cli.drift_cmd, ["--ledger"])

    assert result.exit_code == 0, result.output
    assert "no datasets pending rebuild" in result.output
