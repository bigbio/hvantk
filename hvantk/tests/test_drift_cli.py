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


# --- _stale_datasets timestamp comparison ---------------------------------------------
#
# `rebuilt < entry.get("last_upstream_change", "")` compared ISO-8601 strings
# LEXICOGRAPHICALLY. That happens to agree with chronological order only when both
# timestamps share the same UTC offset and the same sub-second precision -- neither is
# guaranteed for values written across timezones/tools/library versions.


def test_stale_datasets_compares_timestamps_as_real_instants_not_strings():
    """Two counterexamples where the artifact WAS rebuilt after the upstream change,
    but naive string comparison says otherwise.

    "2026-08-23T09:00:00-05:00" (=14:00 UTC) sorts BEFORE
    "2026-08-23T10:00:00+00:00" lexicographically despite being the LATER instant.
    "2026-08-23T10:00:00.500000+00:00" is a later instant than
    "2026-08-23T10:00:00Z", but string comparison here happens to agree only by
    accident of digit count -- the offset case above proves it isn't reliable.
    """
    from hvantk.tools.plugins import drift_cli

    ledger = {
        "a:one": {
            "last_upstream_change": "2026-08-23T10:00:00+00:00",
            "rebuilt_at": "2026-08-23T09:00:00-05:00",
        },
        "b:two": {
            "last_upstream_change": "2026-08-23T10:00:00Z",
            "rebuilt_at": "2026-08-23T10:00:00.500000+00:00",
        },
    }
    assert drift_cli._stale_datasets(ledger) == []


def test_stale_datasets_reports_a_genuinely_stale_dataset():
    from hvantk.tools.plugins import drift_cli

    ledger = {
        "c:three": {
            "last_upstream_change": "2026-08-23T10:00:00+00:00",
            "rebuilt_at": "2026-08-20T00:00:00+00:00",
        },
    }
    assert [name for name, _ in drift_cli._stale_datasets(ledger)] == ["c:three"]


def test_stale_datasets_reports_an_unparseable_rebuilt_at_as_stale():
    """An unparseable timestamp must be reported, not hidden: a false "stale" merely
    prompts someone to look, but a false "fresh" would hide a genuinely stale artifact
    from the report."""
    from hvantk.tools.plugins import drift_cli

    ledger = {
        "d:four": {
            "last_upstream_change": "2026-08-23T10:00:00+00:00",
            "rebuilt_at": "not-a-timestamp",
        },
    }
    assert [name for name, _ in drift_cli._stale_datasets(ledger)] == ["d:four"]


# --- --mark-rebuilt --------------------------------------------------------------------
#
# The ledger records that a dataset drifted, but nothing ever cleared `rebuilt_at` --
# every drifted dataset stayed pending-rebuild forever and `--ledger` could never report
# "nothing pending". `--mark-rebuilt <dataset>` is the missing other half: it records
# that a human (or a follow-up pipeline run) actually rebuilt the artifact.


def test_mark_rebuilt_clears_a_dataset_from_the_stale_list(tmp_path, monkeypatch):
    from hvantk.tools.plugins import drift_cli

    ledger = tmp_path / "drift_ledger.json"
    ledger.write_text(json.dumps({
        "clinvar:variants": {"last_upstream_change": "2026-08-23T00:00:00+00:00",
                             "accepted_in": "PR #288", "signal": "routine",
                             "rebuilt_at": None},
    }))
    monkeypatch.setattr(drift_cli, "LEDGER_PATH", ledger)

    result = CliRunner().invoke(drift_cli.drift_cmd, ["--mark-rebuilt", "clinvar:variants"])
    assert result.exit_code == 0, result.output

    stale = CliRunner().invoke(drift_cli.drift_cmd, ["--ledger"])
    assert stale.exit_code == 0
    assert "no datasets pending rebuild" in stale.output
    assert "clinvar:variants" not in stale.output


def test_mark_rebuilt_unknown_dataset_errors_without_creating_an_entry(tmp_path, monkeypatch):
    """Marking a dataset that never drifted (so it has no ledger row) must fail loudly,
    not silently fabricate a row the drift bot never wrote."""
    from hvantk.tools.plugins import drift_cli

    ledger = tmp_path / "drift_ledger.json"
    ledger.write_text("{}")
    monkeypatch.setattr(drift_cli, "LEDGER_PATH", ledger)

    result = CliRunner().invoke(drift_cli.drift_cmd, ["--mark-rebuilt", "does:not:exist"])

    assert result.exit_code != 0
    assert "does:not:exist" in (result.output or "")
    assert json.loads(ledger.read_text()) == {}


def test_mark_rebuilt_preserves_every_other_entry_untouched(tmp_path, monkeypatch):
    from hvantk.tools.plugins import drift_cli

    ledger = tmp_path / "drift_ledger.json"
    other_entry = {"last_upstream_change": "2026-08-01T00:00:00+00:00",
                    "accepted_in": "PR #286", "signal": "routine",
                    "rebuilt_at": "2026-08-20T00:00:00+00:00"}
    ledger.write_text(json.dumps({
        "clinvar:variants": {"last_upstream_change": "2026-08-23T00:00:00+00:00",
                             "accepted_in": "PR #288", "signal": "routine",
                             "rebuilt_at": None},
        "hgnc:lookup": other_entry,
    }))
    monkeypatch.setattr(drift_cli, "LEDGER_PATH", ledger)

    result = CliRunner().invoke(drift_cli.drift_cmd, ["--mark-rebuilt", "clinvar:variants"])
    assert result.exit_code == 0, result.output

    on_disk = json.loads(ledger.read_text())
    assert on_disk["hgnc:lookup"] == other_entry
    assert on_disk["clinvar:variants"]["accepted_in"] == "PR #288"
    assert on_disk["clinvar:variants"]["signal"] == "routine"
    assert on_disk["clinvar:variants"]["rebuilt_at"] is not None


def test_mark_rebuilt_and_ledger_flag_are_mutually_exclusive(tmp_path, monkeypatch):
    from hvantk.tools.plugins import drift_cli

    ledger = tmp_path / "drift_ledger.json"
    ledger.write_text("{}")
    monkeypatch.setattr(drift_cli, "LEDGER_PATH", ledger)

    result = CliRunner().invoke(
        drift_cli.drift_cmd, ["--ledger", "--mark-rebuilt", "clinvar:variants"]
    )
    assert result.exit_code != 0
    # Must be OUR validation catching the combination, not e.g. an unrecognized-option
    # error from click -- the wording should name both flags.
    assert "--ledger" in result.output and "--mark-rebuilt" in result.output, result.output


# --- --ledger must not silently swallow co-occurring flags ----------------------------
#
# `--ledger` used to return before the --all/dataset validation, so `--ledger --all`,
# `--ledger somedataset`, `--ledger --regenerate`, and `--ledger --json` all exited 0
# and printed the same whole-ledger dump, silently discarding whichever other flag was
# passed -- exactly the kind of surprise the --all/dataset mutual-exclusion check below
# already guards against for the non-ledger path.


def test_ledger_flag_rejects_all_flag(tmp_path, monkeypatch):
    from hvantk.tools.plugins import drift_cli

    monkeypatch.setattr(drift_cli, "LEDGER_PATH", tmp_path / "drift_ledger.json")
    result = CliRunner().invoke(drift_cli.drift_cmd, ["--ledger", "--all"])
    assert result.exit_code != 0
    assert "--ledger" in result.output and "--all" in result.output, result.output


def test_ledger_flag_rejects_a_dataset_argument(tmp_path, monkeypatch):
    from hvantk.tools.plugins import drift_cli

    monkeypatch.setattr(drift_cli, "LEDGER_PATH", tmp_path / "drift_ledger.json")
    result = CliRunner().invoke(drift_cli.drift_cmd, ["--ledger", "clinvar:variants"])
    assert result.exit_code != 0
    assert "--ledger" in result.output, result.output


def test_ledger_flag_rejects_regenerate(tmp_path, monkeypatch):
    from hvantk.tools.plugins import drift_cli

    monkeypatch.setattr(drift_cli, "LEDGER_PATH", tmp_path / "drift_ledger.json")
    result = CliRunner().invoke(drift_cli.drift_cmd, ["--ledger", "--regenerate"])
    assert result.exit_code != 0
    assert "--ledger" in result.output and "--regenerate" in result.output, result.output


def test_ledger_flag_rejects_json(tmp_path, monkeypatch):
    from hvantk.tools.plugins import drift_cli

    monkeypatch.setattr(drift_cli, "LEDGER_PATH", tmp_path / "drift_ledger.json")
    result = CliRunner().invoke(drift_cli.drift_cmd, ["--ledger", "--json"])
    assert result.exit_code != 0
    assert "--ledger" in result.output and "--json" in result.output, result.output
