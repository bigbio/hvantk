"""The scheduled drift bot must not report success when it opened nothing.

``.github/scripts/drift_to_pr.py`` turns a ``hvantk drift --all --json`` report into one
draft PR per drifted dataset. It had no tests, and it caught every ``gh`` failure per
dataset and still returned 0. With the repository's "Allow GitHub Actions to create and
approve pull requests" setting off, the daily job pushed 11 branches, failed all 11
``gh pr create`` calls, and went green -- so nothing surfaced the fact that the drift
loop had stopped delivering.

Loaded by path because ``.github/scripts`` is not an importable package.
"""

from __future__ import annotations

import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parents[2] / ".github" / "scripts" / "drift_to_pr.py"


def _load_module():
    """Import the bot helper from its path in .github/scripts (not a package)."""
    spec = importlib.util.spec_from_file_location("drift_to_pr", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def drift_to_pr():
    """The loaded bot helper module, imported once per test module."""
    return _load_module()


def _write_report(tmp_path: Path, entries: list[dict]) -> Path:
    """Write a `hvantk drift --all --json`-shaped report and return its path."""
    path = tmp_path / "drift_report.json"
    path.write_text(json.dumps(entries))
    return path


DRIFTED = {
    "dataset_name": "clinvar:variants",
    "status": "drifted",
    "observed": {"source_version": "new"},
    "expected": {"source_version": "old"},
    "diff": {"added": {}, "removed": {}, "changed": {}},
    "probe_error": None,
}


def test_exits_nonzero_when_a_pr_cannot_be_opened(drift_to_pr, tmp_path, monkeypatch):
    """The regression: drift detected, branch pushed, `gh pr create` refused, exit 0.

    Reproduces the exact stderr the repository setting produced in run 30255526454.
    """
    report = _write_report(tmp_path, [DRIFTED])

    def _explode(cmd, **kwargs):
        """Succeed for every git call; fail only on `gh pr create`."""
        if cmd[:3] == ["gh", "pr", "create"]:
            raise subprocess.CalledProcessError(
                1,
                cmd,
                output="",
                stderr="GitHub Actions is not permitted to create or approve "
                "pull requests (createPullRequest)",
            )
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(drift_to_pr, "_run", _explode)
    monkeypatch.setattr(drift_to_pr, "pr_exists_for_branch", lambda b, dry_run: None)
    # The real emptiness check shells out to git; the branch under test is what
    # happens after a commit was made.
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *a, **k: subprocess.CompletedProcess(a[0] if a else [], 1),
    )

    rc = drift_to_pr.main(["--report", str(report), "--base-branch", "dev"])

    assert rc == 1


def test_exits_zero_when_nothing_drifted(drift_to_pr, tmp_path):
    """Stub and probe_failed are signals, not failures, and open no pull request."""
    report = _write_report(
        tmp_path,
        [
            {"dataset_name": "gencc:submissions", "status": "clean"},
            {
                "dataset_name": "gevir:metrics",
                "status": "stub",
                "observed": {"reason": "supplementary data; no programmatic URL"},
            },
            {
                "dataset_name": "cptac:phospho",
                "status": "probe_failed",
                "probe_error": "cptac package not installed",
            },
        ],
    )

    assert drift_to_pr.main(["--report", str(report), "--base-branch", "dev"]) == 0


def test_dry_run_makes_no_side_effects_and_succeeds(drift_to_pr, tmp_path):
    """Drift itself is a signal, not a failure: a clean dry run still exits 0."""
    report = _write_report(tmp_path, [DRIFTED])

    rc = drift_to_pr.main(
        ["--report", str(report), "--base-branch", "dev", "--dry-run"]
    )

    assert rc == 0


def test_branch_name_rejects_a_bare_provider(drift_to_pr):
    """A name without a dataset would collide across a provider's datasets."""
    assert drift_to_pr.branch_name_for("clinvar:variants") == "drift/clinvar-variants"
    with pytest.raises(ValueError):
        drift_to_pr.branch_name_for("clinvar")


# --- the daily-churn guard ---------------------------------------------------------
#
# Every open drift PR was force-pushed and its body re-edited once per DAY for a week --
# nine PRs, ~63 notification events, none carrying new information. `--regenerate`
# rewrites `fetched_at` on every run, and the pre-existing emptiness check compares
# against BASE_BRANCH, so for a still-drifted dataset it could never fire. Nothing
# compared the new fingerprint against what the branch already proposed.


def test_fingerprints_match_ignores_only_volatile_keys(drift_to_pr):
    """Same data, later probe -> match. Different source_version -> no match."""
    base = {
        "source_version": "2026-05-15",
        "checksums": {"f.tsv": "abc"},
        "fetched_at": "2026-08-01T00:00:00+00:00",
        "probe_version": 1,
    }
    later = dict(base, fetched_at="2026-08-04T06:00:00+00:00", probe_version=2)
    moved = dict(base, source_version="2026-07-31")

    assert drift_to_pr.fingerprints_match(json.dumps(base), json.dumps(later))
    assert not drift_to_pr.fingerprints_match(json.dumps(base), json.dumps(moved))


def test_fingerprints_match_treats_a_content_change_as_material(drift_to_pr):
    """ClinGen's real case: checksum identical, only `extras.content_length` moved.

    That must still count as a change -- this guard is about suppressing *timestamp*
    churn, not about suppressing upstream content revisions.
    """
    base = {"checksums": {"f.csv": "abc"}, "extras": {"content_length": "1113685"},
            "fetched_at": "2026-08-01T00:00:00+00:00"}
    grown = {"checksums": {"f.csv": "abc"}, "extras": {"content_length": "1114672"},
             "fetched_at": "2026-08-04T00:00:00+00:00"}

    assert not drift_to_pr.fingerprints_match(json.dumps(base), json.dumps(grown))


def test_unparseable_fingerprint_is_never_treated_as_matching(drift_to_pr):
    """"I cannot tell" must mean "push", never "skip".

    Returning True here would silently suppress a real drift PR whenever a fingerprint
    was malformed -- the failure mode this whole script exists to avoid.
    """
    assert not drift_to_pr.fingerprints_match("{not json", "{}")
    assert not drift_to_pr.fingerprints_match("{}", "")


def test_ignored_keys_stay_in_sync_with_the_package(drift_to_pr):
    """The script duplicates PROBE_FINGERPRINT_IGNORED_KEYS because it runs from a
    checkout where hvantk may not be importable. Duplication needs a guard, or the two
    drift apart and the bot starts disagreeing with the detector."""
    from hvantk.core.plugin.api import PROBE_FINGERPRINT_IGNORED_KEYS

    assert drift_to_pr.FINGERPRINT_IGNORED_KEYS == PROBE_FINGERPRINT_IGNORED_KEYS


def test_branch_needs_update_is_true_when_the_branch_is_new(drift_to_pr, monkeypatch):
    """No branch yet -> always push. Cheap, but it is the common first-drift path."""
    monkeypatch.setattr(drift_to_pr, "remote_branch_exists", lambda b, dry_run: False)
    assert drift_to_pr.branch_needs_update("drift/x-y", dry_run=False) is True
