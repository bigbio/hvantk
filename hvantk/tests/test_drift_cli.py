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
    assert "unknown dataset" in result.output.lower() or "unknown dataset" in (
        result.stderr or ""
    )


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
    ledger.write_text(
        json.dumps(
            {
                "clinvar:variants": {
                    "last_upstream_change": "2026-08-23T00:00:00+00:00",
                    "accepted_in": "PR #288",
                    "signal": "routine",
                    "rebuilt_at": None,
                },
                "hgnc:lookup": {
                    "last_upstream_change": "2026-08-01T00:00:00+00:00",
                    "accepted_in": "PR #286",
                    "signal": "routine",
                    "rebuilt_at": "2026-08-20T00:00:00+00:00",
                },
            }
        )
    )
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
    ledger.write_text(
        json.dumps(
            {
                "clinvar:variants": {
                    "last_upstream_change": "2026-08-23T00:00:00+00:00",
                    "accepted_in": "PR #288",
                    "signal": "routine",
                    "rebuilt_at": None,
                },
            }
        )
    )
    monkeypatch.setattr(drift_cli, "LEDGER_PATH", ledger)

    result = CliRunner().invoke(
        drift_cli.drift_cmd, ["--mark-rebuilt", "clinvar:variants"]
    )
    assert result.exit_code == 0, result.output

    stale = CliRunner().invoke(drift_cli.drift_cmd, ["--ledger"])
    assert stale.exit_code == 0
    assert "no datasets pending rebuild" in stale.output
    assert "clinvar:variants" not in stale.output


def test_mark_rebuilt_unknown_dataset_errors_without_creating_an_entry(
    tmp_path, monkeypatch
):
    """Marking a dataset that never drifted (so it has no ledger row) must fail loudly,
    not silently fabricate a row the drift bot never wrote."""
    from hvantk.tools.plugins import drift_cli

    ledger = tmp_path / "drift_ledger.json"
    ledger.write_text("{}")
    monkeypatch.setattr(drift_cli, "LEDGER_PATH", ledger)

    result = CliRunner().invoke(
        drift_cli.drift_cmd, ["--mark-rebuilt", "does:not:exist"]
    )

    assert result.exit_code != 0
    assert "does:not:exist" in (result.output or "")
    assert json.loads(ledger.read_text()) == {}


def test_mark_rebuilt_preserves_every_other_entry_untouched(tmp_path, monkeypatch):
    from hvantk.tools.plugins import drift_cli

    ledger = tmp_path / "drift_ledger.json"
    other_entry = {
        "last_upstream_change": "2026-08-01T00:00:00+00:00",
        "accepted_in": "PR #286",
        "signal": "routine",
        "rebuilt_at": "2026-08-20T00:00:00+00:00",
    }
    ledger.write_text(
        json.dumps(
            {
                "clinvar:variants": {
                    "last_upstream_change": "2026-08-23T00:00:00+00:00",
                    "accepted_in": "PR #288",
                    "signal": "routine",
                    "rebuilt_at": None,
                },
                "hgnc:lookup": other_entry,
            }
        )
    )
    monkeypatch.setattr(drift_cli, "LEDGER_PATH", ledger)

    result = CliRunner().invoke(
        drift_cli.drift_cmd, ["--mark-rebuilt", "clinvar:variants"]
    )
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
    assert "--ledger" in result.output and "--mark-rebuilt" in result.output, (
        result.output
    )


# --- --ledger must not silently swallow co-occurring flags ----------------------------
#
# `--ledger` used to return before the --all/dataset validation, so `--ledger --all`,
# `--ledger somedataset`, `--ledger --regenerate`, and `--ledger --json` all exited 0
# and printed the same whole-ledger dump, silently discarding whichever other flag was
# passed -- exactly the kind of surprise drift_cmd's own `--all`/dataset
# mutual-exclusion check (further down in drift_cli.py) already guards against on the
# non-ledger path.


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
    assert "--ledger" in result.output and "--regenerate" in result.output, (
        result.output
    )


def test_ledger_flag_rejects_json(tmp_path, monkeypatch):
    from hvantk.tools.plugins import drift_cli

    monkeypatch.setattr(drift_cli, "LEDGER_PATH", tmp_path / "drift_ledger.json")
    result = CliRunner().invoke(drift_cli.drift_cmd, ["--ledger", "--json"])
    assert result.exit_code != 0
    assert "--ledger" in result.output and "--json" in result.output, result.output


# --- #351: a unit that never registered must become a row and an exit code -----------
#
# The mechanism #351 added had no test at all: `grep -c 'load_error\|probe_failed'` over
# this file returned 0 until these landed. `test_plugin_loader.py` pins the logging half,
# but nothing pinned the synthetic JSON row or EXIT_PROBE_FAILED -- so a refactor could
# delete the behaviour and leave every signal green, which is verbatim the failure #351
# documents. This module's sibling says it best: a check that cannot fail is decoration.


def _registry_with_load_error(monkeypatch, unit="brokenprov"):
    """A registry holding one healthy dataset and one recorded load failure."""
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    reg._record_load_error(unit, plugin_loader.PluginLoadError("boom"))
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)
    return reg


def test_load_error_becomes_a_probe_failed_row_in_json_mode(monkeypatch):
    """drift-health.yml greps the JSON for status == 'probe_failed' and reads
    dataset_name. A unit that never registered contributes no DriftResult, so without
    this it is absent from the report entirely and the workflow prints 'all probes
    healthy' over a silently smaller set."""
    _registry_with_load_error(monkeypatch)
    result = CliRunner().invoke(drift_cmd, ["--all", "--json"])

    # .stdout, not .output: the command also writes a WARNING to stderr, and CliRunner
    # folds stderr into .output. Machine-readable stdout staying clean is the point.
    rows = json.loads(result.stdout)
    failed = [r for r in rows if r["status"] == "probe_failed"]
    assert [r["dataset_name"] for r in failed] == ["brokenprov"], rows
    assert "plugin failed to load" in failed[0]["probe_error"]
    # Same shape as a genuinely failing probe, so no consumer learns a new key.
    assert set(failed[0]) == set(rows[0])


def test_load_error_sets_exit_probe_failed_even_when_every_probe_is_clean(monkeypatch):
    """The exit code is the other half: drift.yml and drift-health.yml both gate on
    `rc > 2`, so 2 reads as 'probe trouble' rather than 'the command broke'."""
    _registry_with_load_error(monkeypatch)
    result = CliRunner().invoke(drift_cmd, ["--all"])
    assert result.exit_code == 2  # EXIT_PROBE_FAILED, despite fake:default being clean


def test_load_error_is_named_by_provider_not_by_absolute_path(tmp_path, monkeypatch):
    """The recorded unit goes verbatim into the drift-health issue body, and
    drift_to_pr.py splits it on ':' to look up maintainers. An absolute path would leak
    the build machine's layout (/home/runner/work/...) and hand the lookup a directory.
    """
    broken = tmp_path / "skills" / "wonky"
    broken.mkdir(parents=True)
    (broken / "plugin.yaml").write_text("api_version: 2\nname: wonky\ndatasets: nope\n")

    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(broken)

    units = [u for u, _ in reg.load_errors()]
    assert units == ["wonky"], units
    assert not any(str(tmp_path) in u for u in units), units


def test_a_provider_directory_with_no_manifest_is_recorded_not_skipped(tmp_path):
    """A manifest that fails to LOAD was recorded; one that fails to EXIST was not, so
    the provider left the registry with no row, no warning and exit 0. The loader accepts
    only `plugin.yaml`, so a `.yml` typo or a packaging glob that stops shipping it lands
    here."""
    root = tmp_path / "skills"
    (root / "ghost").mkdir(parents=True)
    (root / "ghost" / "plugin.yml").write_text("name: ghost\n")  # note: .yml
    (root / "_conventions").mkdir()  # underscore dirs are not providers

    reg = plugin_loader.PluginRegistry()
    reg.load_from_skills_root(root)

    units = [u for u, _ in reg.load_errors()]
    assert units == ["ghost"], units
    assert "_conventions" not in units


def test_drift_all_refuses_to_report_clean_over_an_empty_registry(monkeypatch):
    """An empty sweep is the limiting case of the silently-smaller set: no targets means
    no results, so exit_codes stays {EXIT_CLEAN} and drift.yml (rc > 2) passes on a
    report of []. Green, having checked nothing."""
    monkeypatch.setattr(plugin_loader, "get_registry", plugin_loader.PluginRegistry)
    result = CliRunner().invoke(drift_cmd, ["--all"])
    assert result.exit_code == 2  # EXIT_PROBE_FAILED, not 0


def test_drift_all_rejects_a_domain_that_matches_nothing(monkeypatch):
    """--domain takes a free-form string and its help text names no valid value, so a
    typo silently meant 'check nothing' and exited 0."""
    result = CliRunner().invoke(drift_cmd, ["--all", "--domain", "genomicz"])
    assert result.exit_code != 0
    assert "matches no dataset" in result.output
    assert "genomics" in result.output  # names the real ones


# --- #361: --regenerate must not write a failing probe's output, and must not exit 1 -----


def test_regenerate_reports_a_failing_probe_as_probe_failed_exit_code(
    tmp_path: Path, monkeypatch
):
    import shutil

    from hvantk.core.plugin.api import DriftProbeError

    plugin_dir = tmp_path / "fake_plugin"
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", plugin_dir)
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(plugin_dir)
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)
    spec = reg.get_dataset("fake:default")

    def boom():
        raise DriftProbeError("upstream down")

    object.__setattr__(spec, "drift_probe", boom)
    fp_path = plugin_dir / "tests" / "drift_fingerprint.json"
    before = fp_path.read_text()

    result = CliRunner().invoke(drift_cmd, ["--regenerate", "fake:default"])

    assert result.exit_code == 2, (
        result.output
    )  # EXIT_PROBE_FAILED, not a traceback's 1
    assert "upstream down" in result.output
    assert fp_path.read_text() == before, (
        "a failed probe must not overwrite the baseline"
    )
    assert sorted(p.name for p in (plugin_dir / "tests").iterdir()) == [
        "drift_fingerprint.json"
    ], "no temp sibling left behind"


def test_regenerate_reports_a_non_driftprobe_exception_as_probe_failed_exit_code(
    tmp_path: Path, monkeypatch
):
    """A probe that raises something other than DriftProbeError (KeyError,
    AttributeError, zlib.error, ...) must not escape as a traceback -- exit 2 either
    way, with the underlying exception type named."""
    import shutil

    plugin_dir = tmp_path / "fake_plugin"
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", plugin_dir)
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(plugin_dir)
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)
    spec = reg.get_dataset("fake:default")

    def boom():
        raise KeyError("boom")

    object.__setattr__(spec, "drift_probe", boom)
    fp_path = plugin_dir / "tests" / "drift_fingerprint.json"
    before = fp_path.read_text()

    result = CliRunner().invoke(drift_cmd, ["--regenerate", "fake:default"])

    assert result.exit_code == 2, (
        result.output
    )  # EXIT_PROBE_FAILED, not a traceback's 1
    assert "KeyError" in result.output
    assert fp_path.read_text() == before, (
        "a failed probe must not overwrite the baseline"
    )
    assert sorted(p.name for p in (plugin_dir / "tests").iterdir()) == [
        "drift_fingerprint.json"
    ], "no temp sibling left behind"


def test_regenerate_reports_a_write_failure_as_probe_failed_exit_code(
    tmp_path: Path, monkeypatch
):
    """A probe that succeeds but a write_fingerprint that fails (disk full, permission
    denied, ...) must also exit 2, not let the OSError escape as a traceback's 1."""
    import shutil

    from hvantk.core.plugin import drift_runner

    plugin_dir = tmp_path / "fake_plugin"
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", plugin_dir)
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(plugin_dir)
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)
    fp_path = plugin_dir / "tests" / "drift_fingerprint.json"
    before = fp_path.read_text()

    def boom_write(path, fingerprint):
        raise OSError("disk full")

    monkeypatch.setattr(drift_runner, "write_fingerprint", boom_write)

    result = CliRunner().invoke(drift_cmd, ["--regenerate", "fake:default"])

    assert result.exit_code == 2, (
        result.output
    )  # EXIT_PROBE_FAILED, not a traceback's 1
    assert "disk full" in result.output
    assert fp_path.read_text() == before, (
        "a write failure must not corrupt the baseline"
    )


def test_regenerate_reports_a_fingerprint_serialization_failure_as_probe_failed_exit_code(
    tmp_path: Path, monkeypatch
):
    """A probe returning a dict with a tuple key fails `_coerce_fingerprint`'s
    JSON-serialisability round-trip as a DriftProbeError, uniformly with
    the plain check path -- so it no longer needs to fall through to
    write_fingerprint's own json.dumps and the CLI's generic `except Exception`
    backstop. Exit code and the no-corruption guarantee are unchanged; only the
    message is now the same shape a probe that raises DriftProbeError directly
    would produce.
    """
    import shutil

    plugin_dir = tmp_path / "fake_plugin"
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", plugin_dir)
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(plugin_dir)
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)
    spec = reg.get_dataset("fake:default")

    object.__setattr__(spec, "drift_probe", lambda: {("a", "b"): 1})
    fp_path = plugin_dir / "tests" / "drift_fingerprint.json"
    before = fp_path.read_text()

    result = CliRunner().invoke(drift_cmd, ["--regenerate", "fake:default"])

    assert result.exit_code == 2, (
        result.output
    )  # EXIT_PROBE_FAILED, not a traceback's 1
    assert "fingerprint NOT rewritten" in result.output
    assert "not JSON-serialisable" in result.output
    assert fp_path.read_text() == before, (
        "a serialization failure must not corrupt the baseline"
    )
    assert sorted(p.name for p in (plugin_dir / "tests").iterdir()) == [
        "drift_fingerprint.json"
    ], "no temp sibling left behind"


# --- #364: load errors are scoped to what the caller asked about --------------------
#
# `load_errors` was applied unconditionally, so `hvantk drift clinvar:variants --json`
# emitted probe_failed rows for units the caller never asked about and exited 2 even
# when clinvar was clean; `--all --domain X` filtered list_datasets but not load_errors;
# and `--regenerate` returned before the exit-code block, so it exited 0 with load errors
# present -- while a PluginLoadError from get_dataset escaped as a traceback, exit 1.

_BROKEN_MANIFEST = (
    "api_version: 2\n"
    "name: brokenprov\n"
    "version: 0.1.0\n"
    "datasets:\n"
    "  - name: thing\n"
    "    domain: proteomics\n"
    "    backend: hail\n"
    "    builder: {module: hvantk.nope.missing, function: build_thing}\n"
    "    drift_probe: {module: hvantk.nope.missing, function: fetch_fingerprint}\n"
    "    skill: SKILL.md\n"
    "    tests: {command: pytest, fixture: f, schema_snapshot: s,\n"
    "            row_snapshot: r, drift_fingerprint: d}\n"
)


def _registry_with_broken_dataset(tmp_path, monkeypatch):
    """fake:default (genomics, healthy) plus brokenprov:thing (proteomics, will not bind)."""
    plugin = tmp_path / "brokenprov"
    plugin.mkdir()
    (plugin / "plugin.yaml").write_text(_BROKEN_MANIFEST)
    (plugin / "SKILL.md").write_text("---\nname: x\ndescription: y\n---\n# x\n")
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    reg.load_from_directory(plugin)
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)
    return reg


def _stdout_rows(result):
    return json.loads(result.stdout)


def test_single_dataset_drift_ignores_unrelated_load_errors(monkeypatch):
    _registry_with_load_error(monkeypatch, unit="brokenprov")
    result = CliRunner().invoke(drift_cmd, ["--json", "fake:default"])
    assert result.exit_code == 0, result.output
    assert [r["dataset_name"] for r in _stdout_rows(result)] == ["fake:default"]
    assert "brokenprov" not in result.output


@pytest.mark.parametrize("unit", ["fake", "entry-point:fake"])
def test_single_dataset_drift_reports_its_own_providers_failure(monkeypatch, unit):
    _registry_with_load_error(monkeypatch, unit=unit)
    result = CliRunner().invoke(drift_cmd, ["--json", "fake:default"])
    assert result.exit_code == 2, result.output
    rows = _stdout_rows(result)
    assert {r["dataset_name"] for r in rows} == {"fake:default", unit}
    assert [r["status"] for r in rows if r["dataset_name"] == unit] == ["probe_failed"]


def test_domain_filter_drops_dataset_level_errors_from_other_domains(
    tmp_path, monkeypatch
):
    _registry_with_broken_dataset(tmp_path, monkeypatch)
    result = CliRunner().invoke(drift_cmd, ["--all", "--domain", "genomics", "--json"])
    assert result.exit_code == 0, result.output
    assert [r["dataset_name"] for r in _stdout_rows(result)] == ["fake:default"]
    assert "brokenprov" not in result.output


def test_all_without_domain_still_reports_every_load_error(tmp_path, monkeypatch):
    _registry_with_broken_dataset(tmp_path, monkeypatch)
    result = CliRunner().invoke(drift_cmd, ["--all", "--json"])
    assert result.exit_code == 2, result.output
    assert "brokenprov:thing" in {r["dataset_name"] for r in _stdout_rows(result)}


def test_domain_filter_keeps_provider_level_errors(monkeypatch):
    """A provider whose manifest did not load has no domain to filter on; dropping it
    would recreate the silently-smaller sweep #351 is about."""
    _registry_with_load_error(monkeypatch, unit="ghostprov")
    result = CliRunner().invoke(drift_cmd, ["--all", "--domain", "genomics", "--json"])
    assert result.exit_code == 2, result.output
    assert "ghostprov" in {r["dataset_name"] for r in _stdout_rows(result)}


def test_regenerate_ignores_unrelated_load_errors(tmp_path, monkeypatch):
    """drift_to_pr.py runs --regenerate unattended and discards the staged fingerprint on
    any non-zero exit, so an unrelated broken plugin must not fail every regenerate."""
    import shutil

    plugin_dir = tmp_path / "fake_plugin"
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", plugin_dir)
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(plugin_dir)
    reg._record_load_error("brokenprov", plugin_loader.PluginLoadError("boom"))
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)
    result = CliRunner().invoke(drift_cmd, ["--regenerate", "fake:default"])
    assert result.exit_code == 0, result.output


def test_regenerate_exits_probe_failed_when_the_datasets_own_provider_failed(
    tmp_path, monkeypatch
):
    import shutil

    plugin_dir = tmp_path / "fake_plugin"
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", plugin_dir)
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(plugin_dir)
    reg._record_load_error("entry-point:fake", plugin_loader.PluginLoadError("boom"))
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)
    result = CliRunner().invoke(drift_cmd, ["--regenerate", "fake:default"])
    assert result.exit_code == 2, result.output
    assert "entry-point:fake" in result.output


@pytest.mark.parametrize(
    "args", [["brokenprov:thing"], ["--regenerate", "brokenprov:thing"]]
)
def test_a_dataset_that_failed_to_bind_exits_probe_failed_not_a_traceback(
    tmp_path, monkeypatch, args
):
    """get_dataset re-raises the cached PluginLoadError, which is not a KeyError, so it
    escaped the `except KeyError` and a wrapper reading exit codes saw 1 == 'drifted'."""
    _registry_with_broken_dataset(tmp_path, monkeypatch)
    result = CliRunner().invoke(drift_cmd, args)
    assert result.exit_code == 2, result.output
    assert not isinstance(result.exception, plugin_loader.PluginLoadError), (
        result.exception
    )
    assert "failed to load" in result.output


def test_a_dataset_that_failed_to_bind_gets_a_json_row(tmp_path, monkeypatch):
    _registry_with_broken_dataset(tmp_path, monkeypatch)
    result = CliRunner().invoke(drift_cmd, ["--json", "brokenprov:thing"])
    assert result.exit_code == 2, result.output
    rows = _stdout_rows(result)
    assert [(r["dataset_name"], r["status"]) for r in rows] == [
        ("brokenprov:thing", "probe_failed")
    ]


# --- follow-up: provider-level load errors recorded under the DIRECTORY name must
# still be relevant to a dataset addressed by the manifest's hyphenated `name:` -------
#
# `load_from_skills_root` records `child.name` (the directory) when `plugin.yaml` is
# missing, and `_provider_id_hint` falls back to `plugin_dir.name` too when the
# manifest's own `name` cannot be read. 9 of the 23 in-tree providers have a directory
# name that differs from the declared `name:` only by `_` vs `-` (`gwas_catalog` dir,
# `gwas-catalog` name; `uniprot_ptm` dir, `uniprot-ptm` name; ...), so exact string
# equality in `_relevant_load_errors` silently dropped a directly-relevant failure.


def test_single_dataset_drift_matches_a_provider_level_error_across_underscore_hyphen(
    tmp_path, monkeypatch
):
    """A provider whose real name uses hyphens (as `gwas-catalog` does) can have a
    load failure recorded under its directory's underscored spelling. That unit must
    still be reported as relevant to `hvantk drift gwas-catalog:associations`."""
    import shutil

    plugin_dir = tmp_path / "gwas_catalog"
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", plugin_dir)
    manifest = plugin_dir / "plugin.yaml"
    manifest.write_text(
        manifest.read_text()
        .replace("name: fake", "name: gwas-catalog", 1)
        .replace("name: default", "name: associations", 1)
    )

    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(plugin_dir)
    reg._record_load_error("gwas_catalog", plugin_loader.PluginLoadError("boom"))
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)

    result = CliRunner().invoke(drift_cmd, ["--json", "gwas-catalog:associations"])

    assert result.exit_code == 2, result.output
    assert "gwas_catalog" in result.stderr
    rows = _stdout_rows(result)
    assert {r["dataset_name"] for r in rows} == {
        "gwas-catalog:associations",
        "gwas_catalog",
    }
    assert [r["status"] for r in rows if r["dataset_name"] == "gwas_catalog"] == [
        "probe_failed"
    ]


# --- follow-up: an unknown dataset should hint at `plugins errors` when there is a
# recorded load failure that might explain it -----------------------------------------


def test_unknown_dataset_hints_at_plugins_errors_when_something_failed_to_load(
    monkeypatch,
):
    _registry_with_load_error(monkeypatch, unit="brokenprov")
    result = CliRunner().invoke(drift_cmd, ["does:not:exist"])
    assert result.exit_code == 3, result.output
    assert "plugins errors" in result.stderr


def test_unknown_dataset_prints_no_hint_when_nothing_failed_to_load():
    result = CliRunner().invoke(drift_cmd, ["does:not:exist"])
    assert result.exit_code == 3, result.output
    assert "plugins errors" not in result.stderr


# --- minor: --regenerate --json writes prose to stdout, breaking the machine-readable
# contract; reject the combination up front instead ------------------------------------


def test_regenerate_and_json_are_mutually_exclusive():
    result = CliRunner().invoke(drift_cmd, ["--regenerate", "--json", "fake:default"])
    assert result.exit_code != 0
    assert "--regenerate" in result.output and "--json" in result.output, result.output


# --- item 1: a non-str-keyed fingerprint must not crash the CHECK path's echo --------
#
# `default=str` in `json.dumps(rows, indent=2, default=str)` rescues non-serialisable
# VALUES, not KEYS. A probe returning `{("a", "b"): 1}` (or a bytes/frozenset/Path key)
# passes `_coerce_fingerprint`'s Mapping check, so the drift RUN succeeds and produces a
# DriftResult carrying the bad-keyed dict in `observed` (and, once diffed, possibly in
# `diff`) -- and only blows up later, as an uncaught TypeError, when the CLI tries to
# `json.dumps` it for output. `--regenerate` already turned this into exit 2 via its own
# catch-all; the plain check path had no such backstop.


def _tuple_keyed_probe_registry(monkeypatch):
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(FIXTURE_ROOT / "fake_plugin")
    spec = reg.get_dataset("fake:default")
    object.__setattr__(spec, "drift_probe", lambda: {("a", "b"): 1})
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)
    return reg


def test_drift_json_single_dataset_exits_probe_failed_on_non_str_keyed_fingerprint(
    monkeypatch,
):
    _tuple_keyed_probe_registry(monkeypatch)
    result = CliRunner().invoke(drift_cmd, ["--json", "fake:default"])
    assert result.exit_code == 2, (
        result.output
    )  # EXIT_PROBE_FAILED, not a traceback's 1
    rows = json.loads(result.output)
    assert rows[0]["dataset_name"] == "fake:default"
    assert rows[0]["status"] == "probe_failed"


def test_drift_all_json_exits_probe_failed_on_non_str_keyed_fingerprint(monkeypatch):
    _tuple_keyed_probe_registry(monkeypatch)
    result = CliRunner().invoke(drift_cmd, ["--all", "--json"])
    assert result.exit_code == 2, (
        result.output
    )  # EXIT_PROBE_FAILED, not a traceback's 1
    rows = json.loads(result.output)
    matching = [r for r in rows if r["dataset_name"] == "fake:default"]
    assert matching and matching[0]["status"] == "probe_failed"


def test_drift_human_readable_exits_probe_failed_on_non_str_keyed_fingerprint(
    monkeypatch,
):
    _tuple_keyed_probe_registry(monkeypatch)
    result = CliRunner().invoke(drift_cmd, ["fake:default"])
    assert result.exit_code == 2, (
        result.output
    )  # EXIT_PROBE_FAILED, not a traceback's 1
    assert "probe_failed" in result.output


# --- a probe whose output is JSON-serialisable but not JSON-native must be clean right
# after --regenerate ---------------------------------------------------------------------
#
# `--regenerate` writes the JSON form (tuple -> list); the next check compared the RAW
# object against `json.loads(file)`, and `("a", "b") != ["a", "b"]`, so such a probe
# reported `drifted` forever -- the bot would regenerate it nightly and never converge.


def _tmp_fake_plugin_registry(tmp_path, monkeypatch):
    """The fake plugin on a tmp copy, so --regenerate cannot touch the committed asset."""
    import shutil

    plugin_dir = tmp_path / "fake_plugin"
    shutil.copytree(FIXTURE_ROOT / "fake_plugin", plugin_dir)
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(plugin_dir)
    monkeypatch.setattr(plugin_loader, "get_registry", lambda: reg)
    return reg, plugin_dir


def test_regenerate_then_check_is_clean_for_a_tuple_valued_fingerprint(
    tmp_path, monkeypatch
):
    reg, _ = _tmp_fake_plugin_registry(tmp_path, monkeypatch)
    spec = reg.get_dataset("fake:default")
    object.__setattr__(
        spec, "drift_probe", lambda: {"probe_version": 1, "headers": ("a", "b")}
    )

    regen = CliRunner().invoke(drift_cmd, ["--regenerate", "fake:default"])
    assert regen.exit_code == 0, regen.output

    check = CliRunner().invoke(drift_cmd, ["--json", "fake:default"])
    assert check.exit_code == 0, check.output
    assert [r["status"] for r in _stdout_rows(check)] == ["clean"]


# --- --regenerate must not commit a stub sentinel or a placeholder-shaped result as the
# baseline ----------------------------------------------------------------------------


@pytest.mark.parametrize("label", ["stub", "placeholder"])
def test_regenerate_exits_probe_failed_instead_of_writing_a_non_baseline(
    tmp_path, monkeypatch, label
):
    from hvantk.core.plugin.api import stub_fingerprint

    probe = {
        "stub": lambda: stub_fingerprint("doc-only; no probeable URL"),
        "placeholder": lambda: {"probe_version": 1, "checksums": {"a.tsv": ""}},
    }[label]
    reg, plugin_dir = _tmp_fake_plugin_registry(tmp_path, monkeypatch)
    spec = reg.get_dataset("fake:default")
    object.__setattr__(spec, "drift_probe", probe)
    fp_path = plugin_dir / "tests" / "drift_fingerprint.json"
    before = fp_path.read_text()

    result = CliRunner().invoke(drift_cmd, ["--regenerate", "fake:default"])

    assert result.exit_code == 2, (label, result.output)  # EXIT_PROBE_FAILED
    assert "fingerprint NOT rewritten" in result.output
    assert fp_path.read_text() == before, "the baseline must be untouched"
    assert sorted(p.name for p in (plugin_dir / "tests").iterdir()) == [
        "drift_fingerprint.json"
    ], "no temp sibling left behind"
