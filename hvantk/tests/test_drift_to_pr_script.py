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
    spec = importlib.util.spec_from_file_location("drift_to_pr", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def drift_to_pr():
    return _load_module()


def _write_report(tmp_path: Path, entries: list[dict]) -> Path:
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
    report = _write_report(tmp_path, [DRIFTED])

    def _explode(cmd, **kwargs):
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
