#!/usr/bin/env python3
"""Open or update a draft PR per drifted dataset reported by
``hvantk drift --all --json``.

This script is invoked by ``.github/workflows/drift.yml``. It is also
runnable locally with ``--dry-run`` for offline testing against a synthetic
report; the dry-run path makes zero git/gh side effects.

Drift report shape (one element per dataset)::

    [
      {
        "dataset_name": "clinvar:variants",
        "status": "clean" | "drifted" | "probe_failed" | "stub",
        "observed": {...} | null,
        "expected": {...} | null,
        "diff": {...} | null,
        "probe_error": str | null
      },
      ...
    ]

For each ``status == "drifted"`` entry the script:

1. Computes a branch name ``drift/<provider>-<dataset>`` (the colon in
   ``provider:dataset`` is replaced with ``-`` so it is a valid ref).
2. Branches off ``BASE_BRANCH``.
3. Runs ``hvantk drift --regenerate <provider:dataset>`` to overwrite the
   committed ``drift_fingerprint.json``.
4. Commits with a plain Conventional Commits message.
5. Force-pushes with lease (the branch is bot-owned) -- but only if the branch does
   not already propose the same fingerprint. ``--regenerate`` rewrites ``fetched_at``
   every run, so without that check each open PR was re-pushed and its body re-edited
   once a day forever: nine PRs churned daily for a week, ~63 notifications, none of
   them new information.
6. Opens a draft PR (or updates the body of an existing one).

Probe-failed entries are recorded in the GitHub step summary but never
produce a PR (they are infrastructure failures, not data drift). ``stub``
entries (documentation-only sources with no programmatic probe) are likewise
never turned into PRs; they are surfaced in the summary count so they are not
silently dropped.

Drift and probe-failure are signals, not workflow failures, so neither sets a
nonzero exit. Failing to *act* on a signal is a different matter: if a drifted
dataset gets as far as needing a PR and the ``gh`` call errors, this script
exits 1. It previously caught that error per dataset, logged it, and still
returned 0 -- so when the repository's "Allow GitHub Actions to create and
approve pull requests" setting was off, every run pushed its branches, failed
all 11 ``gh pr create`` calls, and reported success. The daily job was green for
months while opening nothing.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

try:  # PyYAML ships in the conda env; soft-import so --dry-run still works
    import yaml  # type: ignore
except ImportError:  # pragma: no cover - exercised only in stripped envs
    yaml = None


REPO_ROOT = Path(__file__).resolve().parents[2]
SKILLS_ROOT = REPO_ROOT / "hvantk" / "skills"

# Mirrors hvantk.core.plugin.api.PROBE_FINGERPRINT_IGNORED_KEYS. Duplicated rather than
# imported because this script runs from a checkout where the package may not be
# importable, and `_fingerprints_match` must not become a no-op if the import fails --
# a silently-empty ignore set would make every comparison "different" and restore the
# exact churn this guards against. Kept in sync by test_drift_to_pr_script.py.
FINGERPRINT_IGNORED_KEYS = frozenset({"fetched_at", "probe_version"})


# --------------------------------------------------------------------------- #
# Pure helpers (no side effects, easy to reason about / test).
# --------------------------------------------------------------------------- #


def branch_name_for(dataset: str) -> str:
    """``clinvar:variants`` -> ``drift/clinvar-variants``.

    Refuses any name that isn't ``<provider>:<dataset>`` so we don't end up
    with surprises like ``drift/foo`` (no dataset) silently colliding.
    """
    if ":" not in dataset:
        raise ValueError(f"expected '<provider>:<dataset>', got: {dataset!r}")
    return "drift/" + dataset.replace(":", "-")


def strip_ignored(fingerprint: dict) -> dict:
    """Fingerprint minus the keys that change on every probe regardless of upstream."""
    return {k: v for k, v in fingerprint.items() if k not in FINGERPRINT_IGNORED_KEYS}


def fingerprints_match(a: str, b: str) -> bool:
    """True if two fingerprint JSON blobs agree once volatile keys are dropped.

    Unparseable input returns False -- "I cannot tell" must mean "push it", never
    "skip it", or a malformed fingerprint would silently suppress a real drift PR.
    """
    try:
        return strip_ignored(json.loads(a)) == strip_ignored(json.loads(b))
    except (ValueError, AttributeError, TypeError):
        return False


def branch_name_for_signal(datasets: list[str]) -> str:
    """Branch name for a drift signal covering one or more datasets.

    A single dataset keeps its historical name exactly (``drift/<provider>-<dataset>``),
    so existing branches and their open PRs are still matched. A group that is entirely
    one provider's collapses to ``drift/<provider>`` -- stable across runs, and it does
    not privilege whichever member happened to be listed first in the manifest.
    """
    if len(datasets) == 1:
        return branch_name_for(datasets[0])
    providers = {split_dataset(d)[0] for d in datasets}
    if len(providers) == 1:
        return "drift/" + providers.pop()
    # A baseline shared ACROSS providers should be impossible -- the path lives inside
    # one plugin directory -- but fall back to something deterministic rather than
    # picking arbitrarily, so the branch does not move between runs.
    return branch_name_for(sorted(datasets)[0])


def group_by_drift_signal(drifted: list[dict]) -> list[dict]:
    """Collapse drifted entries that report the SAME drift signal into one.

    Datasets sharing a ``fingerprint_path`` share a baseline file, so they cannot have
    independent drift: one upstream event produces N identical reports, N branches
    writing the same file, and N mutually conflicting PRs. `ucsc-cellbrowser` is the
    worked case -- `default`, `adult-ctx` and `dev-ctx` are distinct *schema* variants
    (obs cell-type column `celltype` / `Class` / `Type_v2`), but the probe fingerprints
    the provider-wide catalog, so all three always agree. Merging one made the other two
    conflict.

    Returns one entry per signal, each carrying a ``datasets`` list naming every member,
    so the PR can say what it covers. Entries WITHOUT a fingerprint_path are never
    grouped -- an older report shape, or a runner that did not populate it, must not
    silently collapse unrelated datasets into one PR.

    Order is preserved so branch names stay stable across runs.
    """
    grouped: dict[str, dict] = {}
    out: list[dict] = []
    for entry in drifted:
        path = entry.get("fingerprint_path")
        if not path:
            out.append({**entry, "datasets": [entry.get("dataset_name")]})
            continue
        existing = grouped.get(path)
        if existing is None:
            merged = {**entry, "datasets": [entry.get("dataset_name")]}
            grouped[path] = merged
            out.append(merged)
        else:
            existing["datasets"].append(entry.get("dataset_name"))
    return out


def split_dataset(dataset: str) -> tuple[str, str]:
    provider, _, name = dataset.partition(":")
    return provider, name


def find_plugin_dir(provider: str) -> Path | None:
    candidate = SKILLS_ROOT / provider
    return candidate if candidate.is_dir() else None


def find_skill_md(provider: str, dataset: str) -> Path | None:
    """Best-effort lookup. Multi-dataset providers nest SKILL.md under
    ``<provider>/<dataset>/SKILL.md``; single-dataset providers keep it at
    the provider root. Returns the first existing path or None.
    """
    plugin = find_plugin_dir(provider)
    if plugin is None:
        return None
    nested = plugin / dataset / "SKILL.md"
    if nested.is_file():
        return nested
    flat = plugin / "SKILL.md"
    return flat if flat.is_file() else None


def read_maintainers(provider: str) -> list[str]:
    """Pull ``maintainers:`` out of the provider's plugin.yaml. Returns an
    empty list if the file or field is missing, or PyYAML isn't available.
    """
    if yaml is None:
        return []
    plugin = find_plugin_dir(provider)
    if plugin is None:
        return []
    manifest = plugin / "plugin.yaml"
    if not manifest.is_file():
        return []
    try:
        data = yaml.safe_load(manifest.read_text()) or {}
    except yaml.YAMLError:
        return []
    maint = data.get("maintainers") or []
    return [str(m).strip() for m in maint if str(m).strip()]


_GITHUB_HANDLE = re.compile(r"^[A-Za-z0-9](?:[A-Za-z0-9-]{0,38})$")


def format_cc_line(maintainers: list[str]) -> str | None:
    """Convert ``["alice@example.com", "@bob", "carol"]`` into a
    ``cc @bob @carol`` line. Email-only entries are skipped (we can't mention
    them on GitHub). Returns None if nothing to mention.
    """
    handles: list[str] = []
    for entry in maintainers:
        stripped = entry.lstrip("@")
        if "@" in stripped:  # looks like email
            continue
        if _GITHUB_HANDLE.match(stripped):
            handles.append("@" + stripped)
    if not handles:
        return None
    return "cc " + " ".join(handles)


def build_pr_body(
    dataset: str,
    diff: dict | None,
    skill_md: Path | None,
    maintainers: list[str],
    covers: list[str] | None = None,
) -> str:
    lines: list[str] = []
    others = [d for d in (covers or []) if d and d != dataset]
    lines.append(
        f"Automated drift detection found upstream changes for `{dataset}`."
    )
    if others:
        listed = ", ".join(f"`{d}`" for d in others)
        lines.append("")
        lines.append(
            f"This also covers {listed}. Those datasets declare the same "
            "`drift_fingerprint` baseline and the same probe, so they report one drift "
            "signal between them, not one each -- they are distinct *schema* variants "
            "of the same upstream resource. Previously each opened its own PR proposing "
            "byte-identical content, and merging any one made the rest conflict."
        )
    lines.append("")
    lines.append(
        "This PR regenerates `drift_fingerprint.json` so the test suite "
        "matches the live upstream snapshot. Review the structured diff "
        "below to decide whether this is:"
    )
    lines.append("")
    lines.append(
        "- a compatible upstream update (merge the snapshot bump as-is),"
    )
    lines.append(
        "- a breaking schema change (also update the plugin's `builder.py`),"
    )
    lines.append(
        "- or a spurious probe difference (fix the probe instead of merging)."
    )
    lines.append("")
    if skill_md is not None:
        rel = skill_md.relative_to(REPO_ROOT).as_posix()
        lines.append(f"Plugin skill: [`{rel}`]({rel})")
        lines.append("")
    lines.append("## Diff")
    lines.append("")
    lines.append("```json")
    lines.append(json.dumps(diff or {}, indent=2, sort_keys=True, default=str))
    lines.append("```")
    cc = format_cc_line(maintainers)
    if cc:
        lines.append("")
        lines.append(cc)
    return "\n".join(lines) + "\n"


def pr_title_for(dataset: str, covers: list[str] | None = None) -> str:
    extra = len([d for d in (covers or []) if d and d != dataset])
    if extra:
        provider = split_dataset(dataset)[0]
        return f"chore(drift): {provider} snapshot regeneration ({extra + 1} datasets)"
    return f"chore(drift): {dataset} snapshot regeneration"


def commit_message_for(dataset: str) -> str:
    return f"chore(drift): regenerate fingerprint for {dataset}"


# --------------------------------------------------------------------------- #
# Side-effecting helpers. Each accepts ``dry_run`` and short-circuits cleanly.
# --------------------------------------------------------------------------- #


def _run(cmd: list[str], *, dry_run: bool, check: bool = True) -> subprocess.CompletedProcess:
    """Wrapper around subprocess that prints (and optionally skips) the
    command. ``check=True`` is the default; pass ``check=False`` for probes
    like ``git ls-remote`` whose nonzero exit is meaningful.
    """
    printable = " ".join(cmd)
    if dry_run:
        print(f"[dry-run] {printable}")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")
    print(f"$ {printable}")
    return subprocess.run(
        cmd,
        check=check,
        text=True,
        capture_output=True,
    )


def remote_branch_exists(branch: str, *, dry_run: bool) -> bool:
    if dry_run:
        # In dry-run we have no remote; report False so the code path
        # exercises the "create" branch.
        return False
    result = subprocess.run(
        ["git", "ls-remote", "--heads", "origin", branch],
        check=False,
        text=True,
        capture_output=True,
    )
    # `or ""` rather than a bare .strip(): stdout is None whenever the call was made
    # without capture_output, and a crash here would abort a drift run over a branch
    # existence check.
    return bool((result.stdout or "").strip())


def branch_needs_update(branch: str, *, dry_run: bool) -> bool:
    """Should the staged fingerprints be pushed to ``branch``?

    False only when the branch already exists AND every staged fingerprint is
    materially identical to the one it already carries -- i.e. the sole difference is
    `fetched_at`. Everything else returns True, deliberately: a branch that does not
    exist, a file the branch lacks, an unreadable blob or a git failure all mean "I
    cannot prove this is redundant", and the safe answer is to push. Suppressing a real
    drift PR is far worse than one redundant force-push.
    """
    if dry_run:
        return True
    if not remote_branch_exists(branch, dry_run=dry_run):
        return True

    staged = subprocess.run(
        ["git", "diff", "--cached", "--name-only"],
        check=False, text=True, capture_output=True,
    )
    paths = [p for p in (staged.stdout or "").split("\n") if p.strip()]
    if staged.returncode != 0 or not paths:
        return True

    # The branch ref may not exist locally -- the run checked out BASE_BRANCH, not this
    # one -- so fetch it and read blobs out of FETCH_HEAD rather than assuming
    # origin/<branch> is present.
    fetched = subprocess.run(
        ["git", "fetch", "origin", branch],
        check=False, text=True, capture_output=True,
    )
    if fetched.returncode != 0:
        return True

    for path in paths:
        current = (REPO_ROOT / path).read_text() if (REPO_ROOT / path).exists() else None
        previous = subprocess.run(
            ["git", "show", f"FETCH_HEAD:{path}"],
            check=False, text=True, capture_output=True,
        )
        if current is None or previous.returncode != 0:
            return True
        if not fingerprints_match(current, previous.stdout):
            return True

    return False


def pr_exists_for_branch(branch: str, *, dry_run: bool) -> str | None:
    """Return the existing PR number (as a string) for ``branch`` if any,
    else None. Uses ``gh pr list`` so it tolerates rate limits gracefully.
    """
    if dry_run:
        return None
    result = subprocess.run(
        [
            "gh", "pr", "list",
            "--head", branch,
            "--state", "open",
            "--json", "number",
            "--jq", ".[0].number // empty",
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    out = result.stdout.strip()
    return out or None


# --------------------------------------------------------------------------- #
# Per-dataset workflow.
# --------------------------------------------------------------------------- #


def handle_drifted(
    entry: dict,
    *,
    base_branch: str,
    dry_run: bool,
    step_summary: Path | None,
) -> None:
    dataset = entry["dataset_name"]
    # Datasets sharing this entry's drift signal (see group_by_drift_signal). Defaults
    # to just this one, so a report without the field behaves exactly as before.
    covers = entry.get("datasets") or [dataset]
    provider, dataset_short = split_dataset(dataset)
    branch = branch_name_for_signal(covers)
    skill_md = find_skill_md(provider, dataset_short)
    maintainers = read_maintainers(provider)
    diff = entry.get("diff") or {}

    print(f"\n=== Drifted: {dataset} -> branch {branch} ===")

    # 1) branch off base
    _run(["git", "fetch", "origin", base_branch], dry_run=dry_run, check=False)
    _run(["git", "checkout", "-B", branch, f"origin/{base_branch}"], dry_run=dry_run)

    # 2) regenerate fingerprint via the CLI we already ship
    _run(
        ["python", "-m", "hvantk.hvantk", "drift", "--regenerate", dataset],
        dry_run=dry_run,
    )

    # 3) stage only fingerprint JSONs (avoid accidentally pulling in unrelated
    #    working-tree noise from the runner).
    _run(
        ["git", "add", "hvantk/skills"],
        dry_run=dry_run,
    )
    # No-op commits should not fail the job — check first.
    if dry_run:
        print("[dry-run] (skipping diff/commit emptiness check)")
    else:
        status = subprocess.run(
            ["git", "diff", "--cached", "--quiet"],
            check=False,
        )
        if status.returncode == 0:
            print(
                f"  no fingerprint changes to commit for {dataset}; "
                "drift may have already been addressed. Skipping PR."
            )
            return

        # ...and neither should a re-push that changes nothing but a timestamp.
        #
        # The check above compares against BASE_BRANCH, so for a dataset that is still
        # drifted it can never fire: `--regenerate` rewrites `fetched_at` on every run,
        # which alone guarantees a non-empty diff. The result was that each of the nine
        # open drift PRs got a fresh force-push and a body edit EVERY morning for a
        # week -- ~63 notification events, none of them carrying new information --
        # because nothing compared the new fingerprint against what the branch already
        # proposed. `fetched_at` is excluded from drift comparison
        # (PROBE_FINGERPRINT_IGNORED_KEYS) but still written into the committed file,
        # so it is invisible to the detector and load-bearing for the diff.
        #
        # So: if the branch already exists and already proposes materially the same
        # fingerprint, leave it alone. The PR stays open with its original body; only a
        # genuine upstream change re-pushes.
        if not branch_needs_update(branch, dry_run=dry_run):
            print(
                f"  branch {branch} already proposes this fingerprint "
                f"(only volatile keys differ); leaving it untouched."
            )
            return

    _run(
        ["git", "commit", "-m", commit_message_for(dataset)],
        dry_run=dry_run,
    )

    # 4) push (force-with-lease; branch is bot-managed)
    _run(
        ["git", "push", "--force-with-lease", "--set-upstream", "origin", branch],
        dry_run=dry_run,
    )

    # 5) open or update PR
    body = build_pr_body(dataset, diff, skill_md, maintainers, covers=covers)
    title = pr_title_for(dataset, covers=covers)

    existing = pr_exists_for_branch(branch, dry_run=dry_run)
    if existing:
        print(f"  updating existing PR #{existing}")
        _run(
            [
                "gh", "pr", "edit", existing,
                "--title", title,
                "--body", body,
            ],
            dry_run=dry_run,
        )
    else:
        print("  creating new draft PR")
        _run(
            [
                "gh", "pr", "create",
                "--draft",
                "--base", base_branch,
                "--head", branch,
                "--title", title,
                "--body", body,
            ],
            dry_run=dry_run,
        )

    if dry_run:
        print("\n--- PR body preview ---")
        print(body)

    _summary_line(
        step_summary,
        f"- DRIFT: `{dataset}` -> branch `{branch}` (draft PR opened/updated)",
    )


def handle_probe_failed(
    entry: dict,
    *,
    step_summary: Path | None,
) -> None:
    dataset = entry["dataset_name"]
    err = entry.get("probe_error") or "(no error message)"
    print(f"\n=== Probe failed: {dataset} ===")
    print(f"  {err}")
    _summary_line(
        step_summary,
        f"- PROBE_FAILED: `{dataset}` -- {err}",
    )


def handle_stub(
    entry: dict,
    *,
    step_summary: Path | None,
) -> None:
    """Surface a documentation-only stub in the rendered step summary so it is
    not silently dropped. Stubs never open a PR (they are not data drift) and
    never fail the workflow (they exit clean)."""
    dataset = entry["dataset_name"]
    reason = (entry.get("observed") or {}).get(
        "reason", "documentation-only source; no programmatic probe"
    )
    print(f"\n=== Stub (no real drift detection): {dataset} ===")
    print(f"  {reason}")
    _summary_line(
        step_summary,
        f"- STUB: `{dataset}` -- {reason}",
    )


def _summary_line(step_summary: Path | None, line: str) -> None:
    if step_summary is None:
        return
    try:
        with step_summary.open("a", encoding="utf-8") as f:
            f.write(line + "\n")
    except OSError as exc:
        print(f"warning: could not write to step summary: {exc}", file=sys.stderr)


# --------------------------------------------------------------------------- #
# Entry point.
# --------------------------------------------------------------------------- #


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("/tmp/drift_report.json"),
        help="Path to the JSON report emitted by `hvantk drift --all --json`.",
    )
    parser.add_argument(
        "--base-branch",
        default=os.environ.get("BASE_BRANCH", "dev"),
        help="Base branch the bot branches off and targets PRs against.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print actions but never invoke git/gh. Safe to run locally.",
    )
    args = parser.parse_args(argv)

    if not args.report.is_file():
        print(f"error: report not found at {args.report}", file=sys.stderr)
        return 1

    try:
        report = json.loads(args.report.read_text())
    except json.JSONDecodeError as exc:
        print(f"error: report is not valid JSON: {exc}", file=sys.stderr)
        return 1

    if not isinstance(report, list):
        print(
            f"error: expected a list of entries, got {type(report).__name__}",
            file=sys.stderr,
        )
        return 1

    summary_path_str = os.environ.get("GITHUB_STEP_SUMMARY")
    step_summary = Path(summary_path_str) if summary_path_str else None
    if step_summary:
        _summary_line(step_summary, "# Drift report")
        _summary_line(step_summary, "")

    drifted = [e for e in report if e.get("status") == "drifted"]
    probe_failed = [e for e in report if e.get("status") == "probe_failed"]
    clean = [e for e in report if e.get("status") == "clean"]
    stub = [e for e in report if e.get("status") == "stub"]

    print(
        f"summary: {len(drifted)} drifted, {len(probe_failed)} probe_failed, "
        f"{len(clean)} clean, {len(stub)} stub"
    )

    groups = group_by_drift_signal(drifted)
    if len(groups) < len(drifted):
        print(
            f"  {len(drifted)} drifted dataset(s) share {len(groups)} distinct drift "
            f"signal(s); opening one PR per signal."
        )

    failed: list[str] = []
    for entry in groups:
        try:
            handle_drifted(
                entry,
                base_branch=args.base_branch,
                dry_run=args.dry_run,
                step_summary=step_summary,
            )
        except subprocess.CalledProcessError as exc:
            # One drifted dataset failing to PR shouldn't stop the others, but it
            # must not vanish either: collected here and re-raised as a nonzero
            # exit once every dataset has had its turn.
            failed.append(str(entry.get("dataset_name")))
            print(
                f"error handling {entry.get('dataset_name')}: "
                f"{exc.cmd} exited {exc.returncode}\n"
                f"stdout: {exc.stdout}\nstderr: {exc.stderr}",
                file=sys.stderr,
            )
            _summary_line(
                step_summary,
                f"- ERROR: `{entry.get('dataset_name')}` -- "
                f"{exc.cmd[0]} exited {exc.returncode}",
            )

    for entry in probe_failed:
        handle_probe_failed(entry, step_summary=step_summary)

    for entry in stub:
        handle_stub(entry, step_summary=step_summary)

    if failed:
        joined = ", ".join(sorted(failed))
        print(
            f"\nERROR: {len(failed)} drifted dataset(s) were detected but no pull "
            f"request could be opened for them: {joined}.\n"
            "If the failure above is 'GitHub Actions is not permitted to create or "
            "approve pull requests', enable Settings -> Actions -> General -> "
            "'Allow GitHub Actions to create and approve pull requests'.",
            file=sys.stderr,
        )
        _summary_line(
            step_summary,
            f"\n**{len(failed)} drifted dataset(s) could not be turned into a PR.** "
            "The branches were pushed; the `gh pr create` calls failed.",
        )
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
