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
5. Force-pushes with lease (the branch is bot-owned).
6. Opens a draft PR (or updates the body of an existing one).

Probe-failed entries are recorded in the GitHub step summary but never
produce a PR (they are infrastructure failures, not data drift). ``stub``
entries (documentation-only sources with no programmatic probe) are likewise
never turned into PRs; they are surfaced in the summary count so they are not
silently dropped.

Exit code is always 0 unless a wholly unexpected error escapes; drift /
probe-failed are not workflow failures.
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
) -> str:
    lines: list[str] = []
    lines.append(
        f"Automated drift detection found upstream changes for `{dataset}`."
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


def pr_title_for(dataset: str) -> str:
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
    return bool(result.stdout.strip())


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
    provider, dataset_short = split_dataset(dataset)
    branch = branch_name_for(dataset)
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
    body = build_pr_body(dataset, diff, skill_md, maintainers)
    title = pr_title_for(dataset)

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

    for entry in drifted:
        try:
            handle_drifted(
                entry,
                base_branch=args.base_branch,
                dry_run=args.dry_run,
                step_summary=step_summary,
            )
        except subprocess.CalledProcessError as exc:
            # One drifted dataset failing to PR shouldn't stop the others.
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

    return 0


if __name__ == "__main__":
    sys.exit(main())
