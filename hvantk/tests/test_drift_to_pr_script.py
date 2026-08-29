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
from datetime import datetime, timedelta, timezone
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


@pytest.fixture(autouse=True)
def isolated_ledger(drift_to_pr, monkeypatch, tmp_path):
    """Give every test its own throwaway ledger file, never the real
    ``hvantk/resources/drift_ledger.json``.

    `record_in_ledger` writes straight to `LEDGER_PATH` with `Path.write_text` -- a
    plain file write, not a subprocess call -- so it is invisible to the `_run` /
    `subprocess.run` monkeypatching the rest of this suite already relies on to stay
    side-effect-free (see `_capture_handle_drifted` and `test_exits_nonzero_when_a_pr_
    cannot_be_opened`, which mock exactly those two for the same reason). Without
    this, any test that drives `handle_drifted` / `handle_routine_batch` past both
    anti-churn guards with `dry_run=False` -- several pre-existing tests do -- writes
    real dataset entries into the tracked ledger file on every test run.
    """
    monkeypatch.setattr(drift_to_pr, "LEDGER_PATH", tmp_path / "drift_ledger.json")


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

#: A grouped entry, as `group_by_drift_signal` emits it: several datasets sharing one
#: baseline, with the first as the anchor. `DRIFTED` carries neither `datasets` nor
#: `fingerprint_path`, so it only ever exercises the single-dataset path.
GROUPED = {
    **DRIFTED,
    "dataset_name": "ucsc-cellbrowser:default",
    "datasets": ["ucsc-cellbrowser:default", "ucsc-cellbrowser:adult-ctx"],
    "fingerprint_path": "hvantk/skills/ucsc_cellbrowser/tests/drift_fingerprint.json",
}


def test_exits_nonzero_when_a_pr_cannot_be_opened(drift_to_pr, tmp_path, monkeypatch):
    """The regression: drift detected, branch pushed, `gh pr create` refused, exit 0.

    Reproduces the exact stderr the repository setting produced in run 30255526454.
    """
    report = _write_report(tmp_path, [DRIFTED])

    def _explode(cmd, **kwargs):
        """Succeed for every git call; fail only on `gh pr create`.

        `DRIFTED`'s diff has no schema-key changes, so `classify_risk` reads it as
        routine and this report is routed through `handle_routine_batch`, whose own
        emptiness check (`git diff --cached --quiet`) goes through `_run` -- unlike
        `handle_drifted`'s, which shells out to `subprocess.run` directly. Answer that
        one command with "there IS something staged" (nonzero), matching this test's
        premise -- drift was detected and a branch was pushed -- so the guard does not
        short-circuit before `gh pr create` gets a chance to fail.
        """
        if cmd[:3] == ["gh", "pr", "create"]:
            raise subprocess.CalledProcessError(
                1,
                cmd,
                output="",
                stderr="GitHub Actions is not permitted to create or approve "
                "pull requests (createPullRequest)",
            )
        if cmd == ["git", "diff", "--cached", "--quiet"]:
            return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(drift_to_pr, "_run", _explode)
    monkeypatch.setattr(drift_to_pr, "pr_exists_for_branch", lambda b, dry_run: None)
    # `handle_drifted`'s own emptiness check shells out to `subprocess.run` directly
    # (rather than `_run`); the branch under test is what happens after a commit was
    # made, so it also must see "something is staged" (nonzero).
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
    drift_to_pr, monkeypatch, tmp_path, *, needs_update, pr_exists, staged=True, entry=None
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
        dict(entry or DRIFTED), base_branch="dev", dry_run=False, step_summary=summary
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


def test_grouped_entry_checks_out_the_signal_branch_and_titles_the_group(
    drift_to_pr, monkeypatch, tmp_path
):
    """`handle_drifted`'s grouping wiring, not just the helpers it calls.

    `branch_name_for_signal` and `pr_title_for` are unit-tested on their own, but every
    other test here passes `DRIFTED`, which has no `datasets` and no `fingerprint_path`
    -- so the single-dataset path is the only one they run. Passing `[dataset]` instead
    of `covers`, or dropping `fingerprint_path`, would leave those unit tests green
    while the bot silently opened a per-dataset PR for a grouped signal: the exact
    regression #263 exists to prevent.
    """
    cmds, _, _ = _capture_handle_drifted(
        drift_to_pr, monkeypatch, tmp_path, needs_update=True, pr_exists=False, entry=GROUPED
    )
    joined = [" ".join(c) for c in cmds]

    checkout = next(c for c in joined if c.startswith("git checkout -B"))
    assert "drift/ucsc-cellbrowser" in checkout, checkout
    # The anchor's per-dataset name is what a dropped `covers` would produce.
    assert "drift/ucsc-cellbrowser-default" not in checkout, checkout

    create = next(c for c in joined if c.startswith("gh pr create"))
    assert "(2 datasets)" in create, create


def test_matching_branch_with_no_open_pr_is_never_skipped(drift_to_pr, monkeypatch, tmp_path):
    """A branch can outlive its PR -- closing a PR leaves the head branch, and a run
    whose push succeeded while `gh pr create` failed leaves a branch with no PR at all.
    Skipping on branch content alone would suppress that dataset's drift forever."""
    cmds, _, _ = _capture_handle_drifted(
        drift_to_pr, monkeypatch, tmp_path, needs_update=False, pr_exists=False
    )
    joined = [" ".join(c) for c in cmds]
    assert any(c.startswith("gh pr create") for c in joined), joined


# --- partial-batch contamination: a failure must not leave a dirty working tree -----
#
# `handle_routine_batch` regenerates every routine dataset in a loop, staging only AFTER
# the loop finishes. If dataset N of M fails, datasets before it have already written
# modified fingerprint files to the working tree. `main()` catches the resulting
# CalledProcessError, logs it, and continues into the schema loop, whose `git checkout -B`
# does not touch an already-dirty working tree, and whose `git add hvantk/skills` would
# stage -- and then commit -- those leftovers onto an unrelated PR: a schema-change PR,
# invisible in its diff and its ledger entry. The same risk exists one level up in
# `handle_drifted`'s own single regenerate call, contaminating the NEXT schema-loop entry.


def test_routine_batch_discards_working_tree_on_mid_loop_regenerate_failure(
    drift_to_pr, monkeypatch
):
    """A regenerate failure partway through the batch must discard whatever earlier
    datasets already wrote to the working tree before the exception is allowed to
    propagate -- or `main()`'s except-and-continue leaves it for the schema loop to
    silently absorb."""
    calls: list[list[str]] = []

    def _regen_fails_on_third(cmd, **kwargs):
        calls.append(list(cmd))
        if cmd[:4] == ["python", "-m", "hvantk.hvantk", "drift"] and cmd[-1] == "c:three":
            raise subprocess.CalledProcessError(1, cmd, output="", stderr="boom")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    def _record_subprocess_run(cmd, **kwargs):
        calls.append(list(cmd))
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(drift_to_pr, "_run", _regen_fails_on_third)
    monkeypatch.setattr(subprocess, "run", _record_subprocess_run)

    entries = [
        {"dataset_name": "a:one", "diff": {}},
        {"dataset_name": "b:two", "diff": {}},
        {"dataset_name": "c:three", "diff": {}},
    ]

    with pytest.raises(subprocess.CalledProcessError):
        drift_to_pr.handle_routine_batch(
            entries, base_branch="dev", dry_run=False, step_summary=None
        )

    joined = [" ".join(c) for c in calls]
    # The cleanup must actually run -- and run BEFORE the exception escapes, i.e. as the
    # very next command after the failing regenerate call.
    assert joined[-1] == "git checkout HEAD -- hvantk/skills", joined
    # And the leftover must never reach the index in the first place on this path.
    assert not any(c.startswith("git add") for c in joined), joined


def test_handle_drifted_discards_working_tree_when_regenerate_fails(
    drift_to_pr, monkeypatch
):
    """Same contamination mechanism one level up: `handle_drifted` regenerates only ONE
    dataset, but a failure here still leaves that dataset's fingerprint dirty in the
    working tree. `main()`'s per-entry except logs it and moves on to the NEXT schema
    entry, whose `git checkout -B` inherits the dirty file and whose `git add
    hvantk/skills` would stage it into a commit that has nothing to do with it."""
    calls: list[list[str]] = []

    def _regen_fails(cmd, **kwargs):
        calls.append(list(cmd))
        if cmd[:4] == ["python", "-m", "hvantk.hvantk", "drift"]:
            raise subprocess.CalledProcessError(1, cmd, output="", stderr="boom")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    def _record_subprocess_run(cmd, **kwargs):
        calls.append(list(cmd))
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(drift_to_pr, "_run", _regen_fails)
    monkeypatch.setattr(subprocess, "run", _record_subprocess_run)

    with pytest.raises(subprocess.CalledProcessError):
        drift_to_pr.handle_drifted(
            dict(DRIFTED), base_branch="dev", dry_run=False, step_summary=None
        )

    joined = [" ".join(c) for c in calls]
    assert joined[-1] == "git checkout HEAD -- hvantk/skills", joined
    assert not any(c.startswith("git add") for c in joined), joined


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


def test_same_baseline_different_probes_are_not_grouped(drift_to_pr):
    """The bot must apply the runner's rule, not a weaker one.

    `run_drift_checks` keys on (baseline, probe callable) so that two datasets sharing a
    baseline while declaring DIFFERENT probes -- a manifest error -- are never merged.
    Grouping on the baseline alone here would re-merge them, opening one PR whose
    regeneration covers only the first dataset and leaving the second's drift silently
    unaddressed.
    """
    groups = drift_to_pr.group_by_drift_signal([
        {**_drifted("p:one", UCSC_FP), "probe_ref": "mod:probe_a"},
        {**_drifted("p:two", UCSC_FP), "probe_ref": "mod:probe_b"},
    ])
    assert len(groups) == 2, "different probes must not share a PR"


def test_same_baseline_same_probe_still_groups(drift_to_pr):
    """The complement: the ucsc case must keep collapsing."""
    groups = drift_to_pr.group_by_drift_signal([
        {**_drifted("u:a", UCSC_FP), "probe_ref": "mod:fetch_fingerprint"},
        {**_drifted("u:b", UCSC_FP), "probe_ref": "mod:fetch_fingerprint"},
    ])
    assert len(groups) == 1
    assert groups[0]["datasets"] == ["u:a", "u:b"]


def test_an_undecodable_staged_file_does_not_abort_the_run(
    repo_with_drift_branch, drift_to_pr
):
    """`read_text()` raises UnicodeDecodeError on non-UTF-8, which is not a
    CalledProcessError, so `main` would not catch it and every remaining dataset would be
    skipped. `paths` is every staged path, not only JSON this script wrote."""
    repo, fp = repo_with_drift_branch
    fp.write_bytes(b'\xff\xfe{"source_version": "v1"}')
    _git(repo, "add", "hvantk/skills")

    # Must return a verdict rather than raising. Undecodable content cannot match, so
    # the safe answer is "push".
    assert drift_to_pr.branch_needs_update("drift/p-d", dry_run=False) is True


# --- risk classification ------------------------------------------------------------
#
# The routine/schema split is the whole basis for batching: routine datasets share one
# PR, schema changes get their own. Misclassifying a schema change as routine would bury
# a builder-breaking change inside a batch nobody reads closely.

def test_classify_risk_routine_when_headers_unchanged(drift_to_pr):
    diff = {
        "changed": {
            "extras": {"expected": {"content_length": "100"},
                       "observed": {"content_length": "205"}},
        },
        "added": {}, "removed": {},
    }
    assert drift_to_pr.classify_risk(diff) == "routine"


def test_classify_risk_schema_when_headers_changed(drift_to_pr):
    diff = {
        "changed": {
            "headers": {"expected": {"f.txt": ["a", "b"]},
                        "observed": {"f.txt": ["a", "b", "c"]}},
        },
        "added": {}, "removed": {},
    }
    assert drift_to_pr.classify_risk(diff) == "schema"


def test_classify_risk_schema_when_checksums_changed(drift_to_pr):
    """For header-hashing probes (clingen, gencc, hgnc) the checksum IS the schema
    signal, so a moved checksum is a schema change, not routine content drift."""
    diff = {
        "changed": {
            "checksums": {"expected": {"f.txt": "aaa"}, "observed": {"f.txt": "bbb"}},
        },
        "added": {}, "removed": {},
    }
    assert drift_to_pr.classify_risk(diff) == "schema"


def test_classify_risk_schema_when_keys_added_or_removed(drift_to_pr):
    """A probe that gained or lost a top-level key changed shape; treat as schema so a
    human looks. Cheap to be wrong in this direction."""
    assert drift_to_pr.classify_risk({"added": {"extras": {}}, "removed": {}, "changed": {}}) == "schema"
    assert drift_to_pr.classify_risk({"added": {}, "removed": {"checksums": {}}, "changed": {}}) == "schema"


def test_classify_risk_unparseable_diff_is_schema(drift_to_pr):
    """"I cannot tell" must mean "show a human", never "batch it silently"."""
    assert drift_to_pr.classify_risk(None) == "schema"
    assert drift_to_pr.classify_risk({}) == "schema"


def test_classify_risk_source_version_alone_is_routine(drift_to_pr):
    """A version string moving with no schema signal is routine -- e.g. gtex-eqtl's
    portal version cl361->cl362, which is a website redeploy, not a data schema change."""
    diff = {
        "changed": {"source_version": {"expected": "cl361", "observed": "cl362"}},
        "added": {}, "removed": {},
    }
    assert drift_to_pr.classify_risk(diff) == "routine"


# --- routine batching ---------------------------------------------------------------

def test_routine_datasets_share_one_branch(drift_to_pr):
    """All routine drift lands on a single branch so it becomes one reviewable PR.
    34 PRs in the first month came from one-PR-per-dataset; batching is what takes
    that to ~1 per run."""
    assert drift_to_pr.ROUTINE_BRANCH == "drift/routine-batch"


def test_partition_by_risk_splits_the_report(drift_to_pr):
    drifted = [
        {"dataset_name": "hgnc:lookup",
         "diff": {"changed": {"extras": {}}, "added": {}, "removed": {}}},
        {"dataset_name": "gtex-eqtl:eqtls",
         "diff": {"changed": {"headers": {}}, "added": {}, "removed": {}}},
        {"dataset_name": "clinvar:variants",
         "diff": {"changed": {"extras": {}}, "added": {}, "removed": {}}},
    ]
    routine, schema = drift_to_pr.partition_by_risk(drifted)

    assert [d["dataset_name"] for d in routine] == ["hgnc:lookup", "clinvar:variants"]
    assert [d["dataset_name"] for d in schema] == ["gtex-eqtl:eqtls"]


def test_partition_by_risk_handles_an_empty_report(drift_to_pr):
    assert drift_to_pr.partition_by_risk([]) == ([], [])


def test_handle_routine_batch_is_a_noop_on_empty_input(drift_to_pr, monkeypatch):
    """No routine drift must mean no branch, no commit, no PR -- not an empty PR."""
    calls = []
    monkeypatch.setattr(drift_to_pr, "_run", lambda cmd, **kw: calls.append(cmd))
    drift_to_pr.handle_routine_batch([], base_branch="dev", dry_run=True, step_summary=None)
    assert calls == []


def test_handle_routine_batch_does_not_commit_when_nothing_staged(
    drift_to_pr, monkeypatch, tmp_path
):
    """The no-op-commit guard, mirrored from handle_drifted: if regenerating every
    routine dataset in the batch produces no staged diff at all -- a re-run after a
    manual merge, or a race where the branch already carries these fingerprints --
    `git commit` must not be invoked. Nothing staged means `git commit` fails outright
    (nothing to commit), which would take the whole batch down with an unhandled
    CalledProcessError instead of a clean skip.
    """
    cmds: list[list[str]] = []
    cleanups: list[dict] = []

    def _record(cmd, **kwargs):
        cmds.append(list(cmd))
        # Every command "succeeds", including `git diff --cached --quiet`: a returncode
        # of 0 from THAT specific command is exactly the "nothing staged" signal this
        # test means to drive.
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(drift_to_pr, "_run", _record)
    monkeypatch.setattr(
        drift_to_pr, "_discard_staged_fingerprints", lambda **kw: cleanups.append(kw)
    )

    summary = tmp_path / "summary.md"
    entries = [
        {"dataset_name": "hgnc:lookup", "diff": {}},
        {"dataset_name": "clinvar:variants", "diff": {}},
    ]
    drift_to_pr.handle_routine_batch(
        entries, base_branch="dev", dry_run=False, step_summary=summary
    )

    joined = [" ".join(c) for c in cmds]
    assert any(c == "git diff --cached --quiet" for c in joined), joined
    assert not any(c.startswith("git commit") for c in joined), joined
    assert not any(c.startswith("git push") for c in joined), joined
    assert not any(c.startswith("gh pr") for c in joined), joined
    # Discard must still run: `git diff --cached --quiet` only proves the INDEX
    # matches HEAD, not that the working tree has nothing left over outside it, and
    # `git checkout -B` (the schema loop's very next move) carries the working tree
    # forward regardless.
    assert cleanups == [{"dry_run": False}]
    assert "DRIFT (unchanged)" in summary.read_text()


# --- visibility requirements (V1 ready-for-review, V3 assigned, V4 labelled, V5 table)

def test_pr_is_created_ready_for_review_not_draft(drift_to_pr):
    """V1. A draft PR cannot be merged and is filtered out of most review queues and
    notification defaults, so the old shape hid the work it was asking for -- which is
    why five PRs sat unreviewed for 3 days on 2026-08-28."""
    argv = drift_to_pr.pr_create_argv(
        base_branch="dev", branch="drift/routine-batch",
        title="t", body="b", risk="routine", assignees=["enriquea"],
    )
    assert "--draft" not in argv


def test_pr_create_applies_the_risk_label(drift_to_pr):
    """V4."""
    argv = drift_to_pr.pr_create_argv(
        base_branch="dev", branch="b", title="t", body="b",
        risk="schema", assignees=[],
    )
    assert "--label" in argv
    assert "drift:schema" in argv


def test_pr_create_assigns_maintainers(drift_to_pr):
    """V3. Nobody was assigned, so nothing appeared on anyone's list."""
    argv = drift_to_pr.pr_create_argv(
        base_branch="dev", branch="b", title="t", body="b",
        risk="routine", assignees=["enriquea", "ypriverol"],
    )
    assert "--assignee" in argv
    assert "enriquea,ypriverol" in argv


def test_pr_create_omits_assignee_when_no_maintainers(drift_to_pr):
    """`gh pr create --assignee ''` errors, so the flag must be absent, not empty."""
    argv = drift_to_pr.pr_create_argv(
        base_branch="dev", branch="b", title="t", body="b",
        risk="routine", assignees=[],
    )
    assert "--assignee" not in argv


def test_pr_create_never_auto_merges(drift_to_pr):
    """Explicit project constraint: nothing in this pipeline may auto-merge."""
    argv = drift_to_pr.pr_create_argv(
        base_branch="dev", branch="b", title="t", body="b",
        risk="routine", assignees=[],
    )
    assert "--auto" not in argv
    assert "merge" not in argv


def test_classification_table_marks_schema_changes(drift_to_pr):
    """A reviewer must be able to spot a schema change without reading JSON."""
    table = drift_to_pr.classification_table([
        {"dataset_name": "hgnc:lookup",
         "diff": {"changed": {"extras": {}}, "added": {}, "removed": {}}},
        {"dataset_name": "gtex-eqtl:eqtls",
         "diff": {"changed": {"headers": {}}, "added": {}, "removed": {}}},
    ])
    assert "`hgnc:lookup`" in table
    assert "routine — schema unchanged" in table
    assert "**SCHEMA CHANGE**" in table


def test_batch_body_leads_with_the_table_not_the_json(drift_to_pr):
    """V5. The old body opened with per-dataset JSON, which is why five PRs were
    indistinguishable at a glance."""
    body = drift_to_pr.build_batch_pr_body([
        {"dataset_name": "hgnc:lookup",
         "diff": {"changed": {"extras": {}}, "added": {}, "removed": {}}},
    ])
    assert body.index("| Dataset |") < body.index("```json")


def test_schema_change_pr_is_also_not_a_draft(drift_to_pr):
    """V1 applies to the schema path too -- arguably more so, since that is the PR that
    most needs a human to look at it."""
    import re
    src = (drift_to_pr.__file__ or "")
    assert src, "could not locate the script source"
    text = open(src).read()
    assert "--draft" not in text


# --- default assignee fallback ------------------------------------------------------
#
# `read_maintainers` reads plugin.yaml's `maintainers:` field, which no manifest
# currently declares -- so without a fallback the assignment requirement is dead code
# and drift PRs go on nobody's list, which is half of why five sat unreviewed for
# three days.

def test_resolve_assignees_prefers_declared_maintainers(drift_to_pr, monkeypatch):
    monkeypatch.setattr(drift_to_pr, "DEFAULT_ASSIGNEE", "fallback-user")
    assert drift_to_pr.resolve_assignees(["alice", "bob"]) == ["alice", "bob"]


def test_resolve_assignees_falls_back_when_none_declared(drift_to_pr, monkeypatch):
    monkeypatch.setattr(drift_to_pr, "DEFAULT_ASSIGNEE", "enriquea")
    assert drift_to_pr.resolve_assignees([]) == ["enriquea"]


def test_resolve_assignees_empty_when_no_fallback_configured(drift_to_pr, monkeypatch):
    """`gh pr create --assignee ''` errors, so with nothing configured the result must
    be an empty list (the flag is then omitted entirely), never [''] ."""
    monkeypatch.setattr(drift_to_pr, "DEFAULT_ASSIGNEE", "")
    assert drift_to_pr.resolve_assignees([]) == []


def test_resolve_assignees_rejects_a_malformed_fallback(drift_to_pr, monkeypatch):
    """A junk env value must not become a --assignee argument."""
    monkeypatch.setattr(drift_to_pr, "DEFAULT_ASSIGNEE", "not a valid handle!")
    assert drift_to_pr.resolve_assignees([]) == []


# --- rebuild ledger -----------------------------------------------------------------
#
# Merging a fingerprint accepts a new baseline; without a ledger, the fact that a built
# artifact is now stale survives only in git history. ClinVar gained ~408 KB of variants
# across 2026-08 with nothing recording that a rebuild was due.

def test_ledger_entry_records_the_accepted_change(drift_to_pr):
    entry = drift_to_pr.ledger_entry(
        dataset="clinvar:variants",
        diff={"changed": {"headers": {}}, "added": {}, "removed": {}},
        pr_ref="PR #288",
        now="2026-08-23T16:56:29+00:00",
    )
    assert entry == {
        "last_upstream_change": "2026-08-23T16:56:29+00:00",
        "accepted_in": "PR #288",
        "signal": "schema",
        "rebuilt_at": None,
    }


def test_ledger_entry_signal_matches_classify_risk(drift_to_pr):
    """The ledger's `signal` must agree with the PR's own classification, or the two
    tell different stories about the same change."""
    diff = {"changed": {"extras": {}}, "added": {}, "removed": {}}
    entry = drift_to_pr.ledger_entry(
        dataset="hgnc:lookup", diff=diff, pr_ref="PR #1", now="2026-08-01T00:00:00+00:00"
    )
    assert entry["signal"] == drift_to_pr.classify_risk(diff) == "routine"


def test_ledger_update_preserves_rebuilt_at_of_other_datasets(drift_to_pr):
    """Updating one dataset must not clear another's rebuild record."""
    ledger = {"hgnc:lookup": {"last_upstream_change": "x", "accepted_in": "PR #1",
                              "signal": "routine", "rebuilt_at": "2026-08-01T00:00:00+00:00"}}
    out = drift_to_pr.ledger_update(
        ledger, dataset="clinvar:variants",
        diff={"changed": {"extras": {}}, "added": {}, "removed": {}},
        pr_ref="PR #2", now="2026-08-23T00:00:00+00:00",
    )
    assert out["hgnc:lookup"]["rebuilt_at"] == "2026-08-01T00:00:00+00:00"
    assert out["clinvar:variants"]["signal"] == "routine"


def test_ledger_update_does_not_mutate_its_input(drift_to_pr):
    ledger = {}
    drift_to_pr.ledger_update(
        ledger, dataset="a:b", diff={"changed": {}, "added": {}, "removed": {}},
        pr_ref="PR #1", now="2026-08-01T00:00:00+00:00",
    )
    assert ledger == {}


def test_load_ledger_returns_empty_dict_on_corrupt_file(drift_to_pr, monkeypatch, tmp_path):
    """A broken ledger must not block a drift PR -- it just starts recording afresh."""
    bad = tmp_path / "drift_ledger.json"
    bad.write_text("{ not json")
    monkeypatch.setattr(drift_to_pr, "LEDGER_PATH", bad)
    assert drift_to_pr.load_ledger() == {}


def test_load_ledger_returns_empty_dict_when_absent(drift_to_pr, monkeypatch, tmp_path):
    monkeypatch.setattr(drift_to_pr, "LEDGER_PATH", tmp_path / "nope.json")
    assert drift_to_pr.load_ledger() == {}


def test_load_ledger_returns_empty_dict_on_truthy_non_dict_json(drift_to_pr, monkeypatch, tmp_path):
    """`json.loads(text) or {}` only substitutes `{}` for FALSY JSON (`[]`, `0`, `""`,
    `null`) -- truthy non-dict JSON (a populated list, a bare string) passes straight
    through unchanged, breaking the docstring's "a missing or corrupt file yields {}"
    promise for exactly the shapes a half-written or hand-edited ledger is likely to
    produce. `record_in_ledger` then calls `dict(ledger)` on the result, which raises
    `TypeError` for a list."""
    non_dict = tmp_path / "drift_ledger.json"
    monkeypatch.setattr(drift_to_pr, "LEDGER_PATH", non_dict)

    non_dict.write_text("[1, 2, 3]")
    assert drift_to_pr.load_ledger() == {}

    non_dict.write_text('"x"')
    assert drift_to_pr.load_ledger() == {}


# --- ledger write ORDERING: staged only after both anti-churn guards ---------------
#
# `record_in_ledger`'s docstring explains at length why it must be called and staged
# only AFTER both anti-churn guards (the `git diff --cached --quiet` emptiness check
# and `branch_needs_update`) have already passed: its `last_upstream_change` is
# `datetime.now(...)`, a different value on literally every invocation, so staging it
# before either guard runs would make both report "needs update" unconditionally --
# reviving the #268 incident (nine PRs force-pushed daily, ~63 notifications/week).
#
# Every other test in this module mocks git/gh and asserts on the SET of commands
# issued, or on the final outcome -- none of them observe WHERE `git add
# hvantk/resources/drift_ledger.json` lands relative to the emptiness guard. Moving
# the ledger write to the top of `handle_drifted` -- reintroducing the #268 bug
# outright -- left a previous 71-test version of this suite 71/71 green.


def test_handle_drifted_stages_the_ledger_only_after_the_emptiness_guard(
    drift_to_pr, monkeypatch
):
    """`handle_drifted`'s own emptiness check shells out to `subprocess.run` directly
    (see `_capture_handle_drifted`'s docstring above) rather than through `_run`, so
    BOTH are patched here into one shared, order-preserving list -- patching `_run`
    alone would leave the `git diff --cached --quiet` call invisible and this
    ordering assertion vacuous.
    """
    calls: list[list[str]] = []

    def _record_run(cmd, **kwargs):
        calls.append(list(cmd))
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    def _record_subprocess_run(cmd, **kwargs):
        calls.append(list(cmd))
        # nonzero == "something is staged", the state after a real regenerate. 0
        # would take the early "nothing to commit" return before the ledger is ever
        # touched -- a different path, already covered by
        # test_skip_issues_no_commit_push_or_pr / test_nothing_staged_path_also_discards.
        if cmd == ["git", "diff", "--cached", "--quiet"]:
            return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(drift_to_pr, "_run", _record_run)
    monkeypatch.setattr(subprocess, "run", _record_subprocess_run)
    monkeypatch.setattr(drift_to_pr, "branch_needs_update", lambda b, dry_run: True)
    monkeypatch.setattr(drift_to_pr, "pr_exists_for_branch", lambda b, dry_run: None)

    drift_to_pr.handle_drifted(
        dict(DRIFTED), base_branch="dev", dry_run=False, step_summary=None
    )

    diff_idx = calls.index(["git", "diff", "--cached", "--quiet"])
    add_ledger_idx = calls.index(["git", "add", "hvantk/resources/drift_ledger.json"])
    assert add_ledger_idx > diff_idx, calls


def test_handle_routine_batch_stages_the_ledger_only_after_the_emptiness_guard(
    drift_to_pr, monkeypatch
):
    """The routine-batch mirror. Its emptiness check DOES go through `_run` (unlike
    `handle_drifted`'s), so patching `_run` alone is enough to see both commands."""
    calls: list[list[str]] = []

    def _record_run(cmd, **kwargs):
        calls.append(list(cmd))
        if cmd == ["git", "diff", "--cached", "--quiet"]:
            return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(drift_to_pr, "_run", _record_run)
    monkeypatch.setattr(drift_to_pr, "branch_needs_update", lambda b, dry_run: True)
    monkeypatch.setattr(drift_to_pr, "pr_exists_for_branch", lambda b, dry_run: None)

    entries = [{"dataset_name": "a:one", "diff": {}}, {"dataset_name": "b:two", "diff": {}}]
    drift_to_pr.handle_routine_batch(
        entries, base_branch="dev", dry_run=False, step_summary=None
    )

    diff_idx = calls.index(["git", "diff", "--cached", "--quiet"])
    add_ledger_idx = calls.index(["git", "add", "hvantk/resources/drift_ledger.json"])
    assert add_ledger_idx > diff_idx, calls


# --- stale-PR escalation ------------------------------------------------------------
#
# A drift PR still open after a full regeneration cycle was not acted on. The bot must
# say so on the EXISTING PR, never by opening another: #268 fixed the inverse failure,
# where nine PRs were force-pushed every morning producing ~63 notification events a
# week carrying no new information.


def test_pr_older_than_one_cycle_is_escalated(drift_to_pr):
    """A PR still open after a full regeneration cycle was not acted on."""
    assert drift_to_pr.should_escalate(
        pr_created_at="2026-08-01T06:00:00Z", now="2026-08-16T06:00:00Z"
    ) is True


def test_pr_within_one_cycle_is_not_escalated(drift_to_pr):
    assert drift_to_pr.should_escalate(
        pr_created_at="2026-08-01T06:00:00Z", now="2026-08-10T06:00:00Z"
    ) is False


def test_pr_exactly_at_the_cycle_boundary_is_escalated(drift_to_pr):
    assert drift_to_pr.should_escalate(
        pr_created_at="2026-08-01T06:00:00Z", now="2026-08-15T06:00:00Z"
    ) is True


def test_unparseable_timestamp_does_not_escalate(drift_to_pr):
    """Escalation is a notification; a parse failure must not spam the PR."""
    assert drift_to_pr.should_escalate(pr_created_at="not-a-date", now="2026-08-16T06:00:00Z") is False
    assert drift_to_pr.should_escalate(pr_created_at="2026-08-01T06:00:00Z", now="") is False
    assert drift_to_pr.should_escalate(pr_created_at=None, now="2026-08-16T06:00:00Z") is False


# `maybe_escalate` itself: the plumbing from a PR number to (at most) one `gh pr
# comment`. `_pr_created_at` is stubbed directly rather than faking `subprocess.run`,
# because -- like `pr_exists_for_branch` / `remote_branch_exists` -- it is a
# value-returning read with its own dry_run gate, not a fire-and-forget action routed
# through `_run`.


def test_maybe_escalate_comments_once_when_pr_outlived_a_cycle(drift_to_pr, monkeypatch):
    """The plumbing: a stale createdAt drives exactly one `gh pr comment`."""
    stale = (datetime.now(timezone.utc) - timedelta(days=20)).isoformat()
    monkeypatch.setattr(drift_to_pr, "_pr_created_at", lambda pr, *, dry_run: stale)
    calls: list[list[str]] = []
    monkeypatch.setattr(
        drift_to_pr, "_run",
        lambda cmd, **kw: calls.append(list(cmd)) or subprocess.CompletedProcess(cmd, 0, stdout="", stderr=""),
    )

    drift_to_pr.maybe_escalate("55", dry_run=False)

    assert len(calls) == 1, calls
    assert calls[0][:3] == ["gh", "pr", "comment"], calls
    assert "55" in calls[0], calls


def test_maybe_escalate_does_not_comment_when_pr_is_recent(drift_to_pr, monkeypatch):
    recent = (datetime.now(timezone.utc) - timedelta(days=1)).isoformat()
    monkeypatch.setattr(drift_to_pr, "_pr_created_at", lambda pr, *, dry_run: recent)
    calls: list[list[str]] = []
    monkeypatch.setattr(
        drift_to_pr, "_run",
        lambda cmd, **kw: calls.append(list(cmd)) or subprocess.CompletedProcess(cmd, 0, stdout="", stderr=""),
    )

    drift_to_pr.maybe_escalate("55", dry_run=False)

    assert calls == []


def test_maybe_escalate_does_nothing_when_pr_created_at_is_unknown(drift_to_pr, monkeypatch):
    """`_pr_created_at` returns "" whenever gh could not answer. No timestamp means no
    verdict, so no comment -- never a crash."""
    monkeypatch.setattr(drift_to_pr, "_pr_created_at", lambda pr, *, dry_run: "")
    calls: list[list[str]] = []
    monkeypatch.setattr(
        drift_to_pr, "_run",
        lambda cmd, **kw: calls.append(list(cmd)) or subprocess.CompletedProcess(cmd, 0, stdout="", stderr=""),
    )

    drift_to_pr.maybe_escalate("55", dry_run=False)

    assert calls == []


def test_maybe_escalate_is_a_noop_under_dry_run(drift_to_pr, monkeypatch):
    """A previous agent found `_run` fabricates `returncode=0` with EMPTY (not None)
    stdout under --dry-run. `maybe_escalate` must degrade safely regardless: no crash,
    no false escalation, and -- since the read is a value-returning query like
    `pr_exists_for_branch` / `remote_branch_exists` -- no subprocess call of any kind,
    mirroring how those two also go silent under dry_run without touching `_run`.
    """
    run_calls: list[list[str]] = []
    monkeypatch.setattr(
        drift_to_pr, "_run",
        lambda cmd, **kw: run_calls.append(list(cmd)) or subprocess.CompletedProcess(cmd, 0, stdout="", stderr=""),
    )
    direct_calls: list = []
    monkeypatch.setattr(
        subprocess, "run",
        lambda *a, **k: direct_calls.append(a) or subprocess.CompletedProcess(a[0] if a else [], 0),
    )

    drift_to_pr.maybe_escalate("999", dry_run=True)

    assert run_calls == []
    assert direct_calls == []


def test_pr_created_at_is_empty_under_dry_run(drift_to_pr):
    """Direct unit test of the dry_run gate `maybe_escalate` relies on."""
    assert drift_to_pr._pr_created_at("1", dry_run=True) == ""


# --- escalation must not repeat forever -----------------------------------------------
#
# `should_escalate` is a pure function of the PR's creation time, which never changes:
# once a PR passes 14 days it is True on EVERY subsequent run, and nothing recorded that
# a comment was already posted. A PR left open six months would collect roughly one
# near-identical comment per run -- the notification-fatigue failure #268 fixed for
# force-pushes, recurring here on a fortnightly cadence instead of daily.


def test_pr_has_escalation_comment_true_when_marker_present(drift_to_pr, monkeypatch):
    monkeypatch.setattr(
        subprocess, "run",
        lambda *a, **k: subprocess.CompletedProcess(
            a[0], 0, stdout=f"unrelated\n{drift_to_pr.ESCALATION_MARKER}\nmore text", stderr=""
        ),
    )
    assert drift_to_pr._pr_has_escalation_comment("1", dry_run=False) is True


def test_pr_has_escalation_comment_false_when_absent(drift_to_pr, monkeypatch):
    monkeypatch.setattr(
        subprocess, "run",
        lambda *a, **k: subprocess.CompletedProcess(a[0], 0, stdout="just a normal comment", stderr=""),
    )
    assert drift_to_pr._pr_has_escalation_comment("1", dry_run=False) is False


def test_pr_has_escalation_comment_false_under_dry_run(drift_to_pr, monkeypatch):
    """Mirrors `_pr_created_at`: a value-returning read with its own dry_run gate, no
    subprocess call of any kind under --dry-run."""
    calls: list = []
    monkeypatch.setattr(
        subprocess, "run",
        lambda *a, **k: calls.append(a) or subprocess.CompletedProcess(a[0] if a else [], 0),
    )
    assert drift_to_pr._pr_has_escalation_comment("1", dry_run=True) is False
    assert calls == []


def test_maybe_escalate_includes_the_marker_in_the_comment_body(drift_to_pr, monkeypatch):
    """The marker must actually be IN the posted comment, or the next run's
    `_pr_has_escalation_comment` check can never find it."""
    stale = (datetime.now(timezone.utc) - timedelta(days=20)).isoformat()
    monkeypatch.setattr(drift_to_pr, "_pr_created_at", lambda pr, *, dry_run: stale)
    monkeypatch.setattr(drift_to_pr, "_pr_has_escalation_comment", lambda pr, *, dry_run: False)
    calls: list[list[str]] = []
    monkeypatch.setattr(
        drift_to_pr, "_run",
        lambda cmd, **kw: calls.append(list(cmd)) or subprocess.CompletedProcess(cmd, 0, stdout="", stderr=""),
    )

    drift_to_pr.maybe_escalate("55", dry_run=False)

    assert len(calls) == 1, calls
    body = calls[0][calls[0].index("--body") + 1]
    assert drift_to_pr.ESCALATION_MARKER in body, body


def test_maybe_escalate_is_idempotent_across_repeated_calls(drift_to_pr, monkeypatch):
    """The regression test: simulate two runs against the same stale PR. The first
    finds no existing comments; the second finds the marker that this test's own `gh pr
    comment` stub "posted" on the first call -- exactly as a real escalation comment
    would leave a marker for the next run's `gh pr view` to find. Only ONE `gh pr
    comment` call must ever be issued across both runs.
    """
    stale = (datetime.now(timezone.utc) - timedelta(days=20)).isoformat()
    monkeypatch.setattr(drift_to_pr, "_pr_created_at", lambda pr, *, dry_run: stale)

    posted_comments: list[str] = []

    def _fake_subprocess_run(cmd, **kwargs):
        if cmd[:3] == ["gh", "pr", "view"]:
            return subprocess.CompletedProcess(cmd, 0, stdout="\n".join(posted_comments), stderr="")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", _fake_subprocess_run)

    calls: list[list[str]] = []

    def _record(cmd, **kw):
        calls.append(list(cmd))
        if cmd[:3] == ["gh", "pr", "comment"]:
            posted_comments.append(cmd[cmd.index("--body") + 1])
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(drift_to_pr, "_run", _record)

    drift_to_pr.maybe_escalate("55", dry_run=False)
    drift_to_pr.maybe_escalate("55", dry_run=False)

    comment_calls = [c for c in calls if c[:3] == ["gh", "pr", "comment"]]
    assert len(comment_calls) == 1, calls


# --- wiring: escalation fires only on the "leaving it untouched" skip path ----------
#
# Placed after `_discard_staged_fingerprints` and before `_summary_line`/`return` in
# BOTH `handle_drifted` and `handle_routine_batch` -- that early return is precisely
# "this PR is still sitting there unmerged". The sibling early return in each function
# (nothing staged at all) is a different situation -- there may be no PR yet -- and
# must not escalate.


def test_handle_drifted_skip_path_calls_maybe_escalate(drift_to_pr, monkeypatch, tmp_path):
    escalated: list[tuple] = []
    monkeypatch.setattr(
        drift_to_pr, "maybe_escalate",
        lambda pr, *, dry_run: escalated.append((pr, dry_run)),
    )
    _capture_handle_drifted(
        drift_to_pr, monkeypatch, tmp_path, needs_update=False, pr_exists=True
    )
    assert escalated == [("7", False)]


def test_handle_drifted_update_path_never_escalates(drift_to_pr, monkeypatch, tmp_path):
    """The complement: a branch that DOES get updated is not "left untouched", so it
    must not also escalate."""
    escalated: list[tuple] = []
    monkeypatch.setattr(
        drift_to_pr, "maybe_escalate",
        lambda pr, *, dry_run: escalated.append((pr, dry_run)),
    )
    _capture_handle_drifted(
        drift_to_pr, monkeypatch, tmp_path, needs_update=True, pr_exists=True
    )
    assert escalated == []


def test_handle_drifted_nothing_staged_path_never_escalates(drift_to_pr, monkeypatch, tmp_path):
    """The OTHER early return (nothing to commit) is not "an existing PR left
    untouched" -- there may be no PR at all yet -- so it must not escalate."""
    escalated: list[tuple] = []
    monkeypatch.setattr(
        drift_to_pr, "maybe_escalate",
        lambda pr, *, dry_run: escalated.append((pr, dry_run)),
    )
    _capture_handle_drifted(
        drift_to_pr, monkeypatch, tmp_path, needs_update=True, pr_exists=True, staged=False
    )
    assert escalated == []


def test_handle_routine_batch_skip_path_calls_maybe_escalate(drift_to_pr, monkeypatch):
    """The routine-batch mirror of the handle_drifted wiring above."""
    escalated: list[tuple] = []

    def _record(cmd, **kwargs):
        if cmd == ["git", "diff", "--cached", "--quiet"]:
            return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(drift_to_pr, "_run", _record)
    monkeypatch.setattr(drift_to_pr, "branch_needs_update", lambda b, dry_run: False)
    monkeypatch.setattr(drift_to_pr, "pr_exists_for_branch", lambda b, dry_run: "42")
    monkeypatch.setattr(drift_to_pr, "_discard_staged_fingerprints", lambda **kw: None)
    monkeypatch.setattr(
        drift_to_pr, "maybe_escalate",
        lambda pr, *, dry_run: escalated.append((pr, dry_run)),
    )

    drift_to_pr.handle_routine_batch(
        [{"dataset_name": "hgnc:lookup", "diff": {}}],
        base_branch="dev", dry_run=False, step_summary=None,
    )

    assert escalated == [("42", False)]


def test_nothing_in_the_drift_pipeline_auto_merges():
    """Hard project constraint, enforced rather than conventional.

    Auto-merge was considered and rejected: batching and cadence alone take drift
    volume from ~34 PRs/month to ~2, so auto-merge would only be the step from 2 to 0
    -- and that step would require carving an exception into CLAUDE.md's "PR must be
    approved before merging to `dev`". The 2026-08-28 backlog was a visibility failure,
    not a review-burden one, so auto-merging would route around the problem instead of
    fixing it.

    This asserts on the real files rather than on a constant, because the failure mode
    is someone adding `gh pr merge --auto` to a workflow in six months without reading
    that reasoning.
    """
    import pathlib
    import re

    root = pathlib.Path(__file__).resolve().parents[2] / ".github"
    banned = re.compile(r"gh\s+pr\s+merge|--auto\b|\"merge\"")
    offenders = []
    for path in sorted(root.rglob("*")):
        if path.suffix not in {".yml", ".yaml", ".py"} or not path.is_file():
            continue
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            stripped = line.strip()
            # Prose is fine -- the constraint is documented in several places.
            if stripped.startswith("#") or stripped.startswith("Never emits"):
                continue
            if banned.search(line):
                offenders.append(f"{path.relative_to(root.parent)}:{lineno}: {stripped}")

    assert not offenders, "auto-merge found in the drift pipeline:\n" + "\n".join(offenders)
