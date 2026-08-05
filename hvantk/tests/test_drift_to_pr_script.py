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


# --- one PR per drift SIGNAL, not per dataset --------------------------------------
#
# ucsc-cellbrowser declares three datasets (default / adult-ctx / dev-ctx) that are
# distinct SCHEMA variants -- their obs cell-type column is `celltype`, `Class` and
# `Type_v2` -- but all three share one drift_fingerprint and a zero-argument probe that
# fingerprints the provider-wide catalog. One upstream event therefore produced three
# identical PRs, and merging any one made the other two conflict (#241/#242/#243).


def _drifted(name, path):
    return {"dataset_name": name, "status": "drifted", "fingerprint_path": path,
            "observed": {}, "expected": {}, "diff": {}, "probe_error": None}


UCSC_FP = "hvantk/skills/ucsc_cellbrowser/tests/drift_fingerprint.json"


def test_datasets_sharing_a_baseline_collapse_to_one_entry(drift_to_pr):
    groups = drift_to_pr.group_by_drift_signal([
        _drifted("ucsc-cellbrowser:default", UCSC_FP),
        _drifted("ucsc-cellbrowser:adult-ctx", UCSC_FP),
        _drifted("ucsc-cellbrowser:dev-ctx", UCSC_FP),
    ])

    assert len(groups) == 1
    assert groups[0]["datasets"] == [
        "ucsc-cellbrowser:default",
        "ucsc-cellbrowser:adult-ctx",
        "ucsc-cellbrowser:dev-ctx",
    ]


def test_distinct_baselines_are_never_merged(drift_to_pr):
    """The guard that keeps this from over-collapsing: different file, different PR."""
    groups = drift_to_pr.group_by_drift_signal([
        _drifted("clinvar:variants", "hvantk/skills/clinvar/tests/drift_fingerprint.json"),
        _drifted("hgnc:lookup", "hvantk/skills/hgnc/tests/drift_fingerprint.json"),
    ])
    assert len(groups) == 2


def test_entries_without_a_fingerprint_path_are_never_grouped(drift_to_pr):
    """An older report shape must not silently collapse unrelated datasets into one PR."""
    groups = drift_to_pr.group_by_drift_signal([
        {"dataset_name": "a:x", "status": "drifted"},
        {"dataset_name": "b:y", "status": "drifted"},
    ])
    assert len(groups) == 2
    assert [g["datasets"] for g in groups] == [["a:x"], ["b:y"]]


def test_group_branch_is_provider_level_but_single_datasets_keep_their_name(drift_to_pr):
    """A lone dataset must keep its historical branch, or open PRs stop being matched."""
    assert drift_to_pr.branch_name_for_signal(["clinvar:variants"]) == "drift/clinvar-variants"
    assert drift_to_pr.branch_name_for_signal([
        "ucsc-cellbrowser:default", "ucsc-cellbrowser:adult-ctx",
    ]) == "drift/ucsc-cellbrowser"


def test_group_branch_is_stable_regardless_of_member_order(drift_to_pr):
    """Branch names must not move between runs, or every run orphans yesterday's PR."""
    a = drift_to_pr.branch_name_for_signal(["ucsc:default", "ucsc:adult", "ucsc:dev"])
    b = drift_to_pr.branch_name_for_signal(["ucsc:dev", "ucsc:default", "ucsc:adult"])
    assert a == b


def test_grouped_pr_names_every_dataset_it_covers(drift_to_pr):
    body = drift_to_pr.build_pr_body(
        "ucsc-cellbrowser:default", {}, None, [],
        covers=["ucsc-cellbrowser:default", "ucsc-cellbrowser:adult-ctx"],
    )
    assert "ucsc-cellbrowser:adult-ctx" in body

    title = drift_to_pr.pr_title_for(
        "ucsc-cellbrowser:default",
        covers=["ucsc-cellbrowser:default", "ucsc-cellbrowser:adult-ctx"],
    )
    assert "2 datasets" in title
    # A single dataset keeps the original title verbatim.
    assert drift_to_pr.pr_title_for("clinvar:variants") == (
        "chore(drift): clinvar:variants snapshot regeneration"
    )
# --- the skip path itself, end to end -----------------------------------------------
#
# Adversarial review caught that NO test exercised the behaviour this change adds:
# replacing the whole body of `branch_needs_update` with `return True` -- a complete
# neutering of the fix -- left every test green. These drive handle_drifted and assert
# on the git/gh commands actually issued.


def _capture_handle_drifted(
    drift_to_pr, monkeypatch, tmp_path, *, needs_update, pr_exists, staged=True
):
    """Run handle_drifted with git/gh stubbed.

    Returns (commands issued, cleanup-call count, step-summary text). The cleanup count
    and summary are captured deliberately: an earlier version of this helper recorded
    only `_run` and passed step_summary=None, so deleting the
    `_discard_staged_fingerprints(...)` CALL SITES -- a full revert of the fix they
    belong to -- left the whole suite green. `_discard_staged_fingerprints` calls
    `subprocess.run` directly, not `_run`, so it was invisible here.
    """
    cmds: list[list[str]] = []
    cleanups: list[dict] = []

    def _record(cmd, **kwargs):
        cmds.append(list(cmd))
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(drift_to_pr, "_run", _record)
    monkeypatch.setattr(drift_to_pr, "branch_needs_update", lambda b, dry_run: needs_update)
    monkeypatch.setattr(
        drift_to_pr, "pr_exists_for_branch", lambda b, dry_run: "7" if pr_exists else None
    )
    # Record the KWARGS, not just the fact of a call. `_discard_staged_fingerprints`
    # is a no-op under dry_run=True (it prints and returns), so a spy that discards its
    # arguments cannot tell a real cleanup from a neutered one -- a call site changed to
    # `dry_run=True` would restore the contamination bug with the suite still green.
    monkeypatch.setattr(
        drift_to_pr, "_discard_staged_fingerprints", lambda **kw: cleanups.append(kw)
    )
    # `git diff --cached --quiet` -> returncode 1 means "there is something staged",
    # which is the state after a real regeneration. 0 means nothing to commit.
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *a, **k: subprocess.CompletedProcess(a[0] if a else [], 1 if staged else 0),
    )

    summary = tmp_path / "summary.md"
    drift_to_pr.handle_drifted(
        dict(DRIFTED), base_branch="dev", dry_run=False, step_summary=summary
    )
    return cmds, cleanups, (summary.read_text() if summary.exists() else "")


def test_skip_issues_no_commit_push_or_pr(drift_to_pr, monkeypatch, tmp_path):
    """The behaviour the whole PR exists for. Fails if branch_needs_update is neutered
    to `return True`, which is exactly the hole review found."""
    cmds, cleanups, summary = _capture_handle_drifted(
        drift_to_pr, monkeypatch, tmp_path, needs_update=False, pr_exists=True
    )
    joined = [" ".join(c) for c in cmds]
    assert not any(c.startswith("git commit") for c in joined), joined
    assert not any(c.startswith("git push") for c in joined), joined
    assert not any(c.startswith("gh pr") for c in joined), joined
    # The staged fingerprint MUST be discarded before returning, or it is committed
    # onto the next dataset's branch. Deleting this call site was a silent revert.
    assert cleanups == [{"dry_run": False}], (
        "skip path must discard the staged fingerprint, and must do it for real -- "
        "dry_run=True would print and return, leaving the index contaminated"
    )
    # ...and a still-drifting dataset must not vanish from the rendered report.
    assert "DRIFT (unchanged)" in summary, summary


def test_update_still_commits_pushes_and_edits(drift_to_pr, monkeypatch, tmp_path):
    """The complement: when the branch DOES need updating, nothing is suppressed."""
    cmds, _, _ = _capture_handle_drifted(
        drift_to_pr, monkeypatch, tmp_path, needs_update=True, pr_exists=True
    )
    joined = [" ".join(c) for c in cmds]
    assert any(c.startswith("git commit") for c in joined), joined
    assert any(c.startswith("git push") for c in joined), joined
    assert any(c.startswith("gh pr edit") for c in joined), joined


def test_matching_branch_with_no_open_pr_is_never_skipped(drift_to_pr, monkeypatch, tmp_path):
    """A branch can outlive its PR -- closing a PR leaves the head branch, and a run
    whose push succeeded while `gh pr create` failed leaves a branch with no PR at all.
    Skipping on branch content alone would suppress that dataset's drift forever."""
    cmds, _, _ = _capture_handle_drifted(
        drift_to_pr, monkeypatch, tmp_path, needs_update=False, pr_exists=False
    )
    joined = [" ".join(c) for c in cmds]
    assert any(c.startswith("gh pr create") for c in joined), joined


# --- branch_needs_update against a REAL git repo -------------------------------------
#
# The stubbed tests above drive handle_drifted's branching, but they monkeypatch
# branch_needs_update itself, so they cannot detect that function being wrong. Neutering
# it to `return True` left them all green -- the same vacuity adversarial review found in
# the first version of this test file. This exercises the real thing over real git.


def _git(repo, *args):
    return subprocess.run(
        ["git", *args], cwd=repo, check=True, text=True, capture_output=True
    )


@pytest.fixture
def repo_with_drift_branch(tmp_path, drift_to_pr, monkeypatch):
    """A repo whose `drift/p-d` branch already carries a fingerprint, with `origin`
    pointing back at itself so ls-remote/fetch behave as they do in CI."""
    repo = tmp_path / "repo"
    (repo / "hvantk" / "skills").mkdir(parents=True)
    fp = repo / "hvantk" / "skills" / "drift_fingerprint.json"

    _git(repo.parent, "init", "-q", str(repo))
    _git(repo, "config", "user.email", "t@t")
    _git(repo, "config", "user.name", "t")
    fp.write_text(json.dumps({"source_version": "v1", "fetched_at": "2026-08-01"}))
    _git(repo, "add", "-A")
    _git(repo, "commit", "-qm", "base")
    _git(repo, "branch", "-M", "dev")
    _git(repo, "checkout", "-qb", "drift/p-d")
    _git(repo, "commit", "-q", "--allow-empty", "-m", "branch head")
    _git(repo, "checkout", "-q", "dev")
    _git(repo, "remote", "add", "origin", str(repo))

    monkeypatch.setattr(drift_to_pr, "REPO_ROOT", repo)
    monkeypatch.chdir(repo)
    return repo, fp


def test_branch_needs_update_false_when_only_the_timestamp_moved(
    repo_with_drift_branch, drift_to_pr
):
    """The regression this PR exists for, against real git."""
    repo, fp = repo_with_drift_branch
    fp.write_text(json.dumps({"source_version": "v1", "fetched_at": "2026-08-04"}))
    _git(repo, "add", "hvantk/skills")

    assert drift_to_pr.branch_needs_update("drift/p-d", dry_run=False) is False


def test_branch_needs_update_true_when_the_source_actually_moved(
    repo_with_drift_branch, drift_to_pr
):
    """The signal must survive: a real upstream change still pushes."""
    repo, fp = repo_with_drift_branch
    fp.write_text(json.dumps({"source_version": "v2", "fetched_at": "2026-08-04"}))
    _git(repo, "add", "hvantk/skills")

    assert drift_to_pr.branch_needs_update("drift/p-d", dry_run=False) is True


def test_discard_staged_fingerprints_clears_index_and_worktree(
    repo_with_drift_branch, drift_to_pr
):
    """The contamination fix: a skipped dataset must leave nothing behind for the next
    one. `git checkout -B` does not clear the index, so a leftover staged fingerprint
    would be committed onto an unrelated dataset's branch."""
    repo, fp = repo_with_drift_branch
    original = fp.read_text()
    fp.write_text(json.dumps({"source_version": "LEAKED", "fetched_at": "x"}))
    _git(repo, "add", "hvantk/skills")
    assert _git(repo, "diff", "--cached", "--name-only").stdout.strip()

    drift_to_pr._discard_staged_fingerprints(dry_run=False)

    assert _git(repo, "diff", "--cached", "--name-only").stdout.strip() == ""
    assert fp.read_text() == original, "working tree must be restored too, or the next "\
        "`git add hvantk/skills` re-stages the leak"


def test_two_signals_from_one_provider_do_not_collide_on_one_branch(drift_to_pr):
    """Review finding on #263: keying a group's branch on the provider alone discards the
    baseline path that DEFINES the signal. Multi-dataset providers suffix their baselines
    per dataset, so a provider can own two independent signals; collapsing them onto one
    branch means the second overwrites the first's commit and rewrites its PR."""
    a = drift_to_pr.branch_name_for_signal(
        ["p:one", "p:two"], "hvantk/skills/p/tests/drift_fingerprint.json"
    )
    b = drift_to_pr.branch_name_for_signal(
        ["p:three", "p:four"], "hvantk/skills/p/tests/drift_fingerprint_samples.json"
    )
    assert a != b, f"distinct signals collided on {a}"
    assert a == "drift/p"                 # the unsuffixed baseline keeps the plain name
    assert b == "drift/p-samples"


def test_group_branch_still_stable_across_member_order(drift_to_pr):
    fp = "hvantk/skills/ucsc_cellbrowser/tests/drift_fingerprint.json"
    assert drift_to_pr.branch_name_for_signal(["u:a", "u:b", "u:c"], fp) == \
           drift_to_pr.branch_name_for_signal(["u:c", "u:a", "u:b"], fp)

def test_nothing_staged_path_also_discards(drift_to_pr, monkeypatch, tmp_path):
    """The OTHER early return after `git add`. Both must clean up, or whichever is left
    uncovered reintroduces contamination on its own path."""
    _, cleanups, _ = _capture_handle_drifted(
        drift_to_pr, monkeypatch, tmp_path, needs_update=True, pr_exists=True, staged=False
    )
    assert cleanups == [{"dry_run": False}], (
        "the 'nothing to commit' return must discard too, and for real"
    )
