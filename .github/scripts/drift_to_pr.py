#!/usr/bin/env python3
"""Open or update a PR per drifted dataset reported by
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
6. Opens a PR ready for review -- labelled with its risk classification and
   assigned to the plugin's declared maintainers, if any -- or updates the body
   of an existing one.

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
from datetime import datetime, timedelta, timezone
from pathlib import Path

try:  # PyYAML ships in the conda env; soft-import so --dry-run still works
    import yaml  # type: ignore
except ImportError:  # pragma: no cover - exercised only in stripped envs
    yaml = None


REPO_ROOT = Path(__file__).resolve().parents[2]
SKILLS_ROOT = REPO_ROOT / "hvantk" / "skills"
LEDGER_PATH = REPO_ROOT / "hvantk" / "resources" / "drift_ledger.json"

# Mirrors hvantk.core.plugin.api.PROBE_FINGERPRINT_IGNORED_KEYS. Duplicated rather than
# imported because this script runs from a checkout where the package may not be
# importable, and `_fingerprints_match` must not become a no-op if the import fails --
# a silently-empty ignore set would make every comparison "different" and restore the
# exact churn this guards against. Kept in sync by test_drift_to_pr_script.py.
FINGERPRINT_IGNORED_KEYS = frozenset({"fetched_at", "probe_version", "informational"})


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


# Fingerprint keys that carry the SCHEMA signal. `headers` is the column list;
# `checksums` is a hash of the column-header row for the header-hashing probes
# (clingen, gencc, hgnc), so a moved checksum means the columns moved.
SCHEMA_KEYS = frozenset({"headers", "checksums"})


def classify_risk(diff: dict | None) -> str:
    """Classify a drift diff as ``"routine"`` or ``"schema"``.

    Routine means the schema signal is unchanged and only content/version moved --
    safe to batch with other routine datasets into one PR. Schema means the column
    list or header hash moved, or the fingerprint gained/lost a top-level key, and
    the plugin's builder.py may need a change.

    Defaults to ``"schema"`` for anything it cannot read. Misclassifying a schema
    change as routine would bury it in a batch; the reverse just opens one extra PR.
    """
    if not isinstance(diff, dict) or not diff:
        return "schema"
    if diff.get("added") or diff.get("removed"):
        return "schema"
    changed = diff.get("changed") or {}
    if not isinstance(changed, dict):
        return "schema"
    if SCHEMA_KEYS & set(changed):
        return "schema"
    return "routine"


# All routine drift shares one branch, so it becomes one reviewable PR per run rather
# than one per dataset. Schema changes keep their own per-dataset branches.
ROUTINE_BRANCH = "drift/routine-batch"


def partition_by_risk(drifted: list[dict]) -> tuple[list[dict], list[dict]]:
    """Split drifted entries into (routine, schema), preserving report order."""
    routine, schema = [], []
    for entry in drifted:
        target = routine if classify_risk(entry.get("diff")) == "routine" else schema
        target.append(entry)
    return routine, schema


# --------------------------------------------------------------------------- #
# Rebuild ledger. A fingerprint bump is also the signal that a built artifact may now
# be stale; without this, that fact survives only in git history. Recorded in the same
# commit as the fingerprints it describes -- see record_in_ledger below for why it must
# be written and staged only once handle_drifted / handle_routine_batch have already
# decided a commit is happening, never earlier.
# --------------------------------------------------------------------------- #


def ledger_entry(*, dataset: str, diff: dict | None, pr_ref: str, now: str) -> dict:
    """One ledger row. ``rebuilt_at`` starts None and is set by whoever rebuilds;
    a value older than ``last_upstream_change`` is the stale-artifact condition.

    ``dataset`` is accepted for a call signature symmetric with ``ledger_update``
    (and to self-document each call site); it is not part of the row itself, since
    the dataset name is already the ledger's key at the ``ledger_update`` level.
    """
    return {
        "last_upstream_change": now,
        "accepted_in": pr_ref,
        "signal": classify_risk(diff),
        "rebuilt_at": None,
    }


def ledger_update(ledger: dict, *, dataset: str, diff: dict | None, pr_ref: str, now: str) -> dict:
    """Return a copy of ``ledger`` with ``dataset`` updated. Never touches other rows."""
    out = dict(ledger)
    out[dataset] = ledger_entry(dataset=dataset, diff=diff, pr_ref=pr_ref, now=now)
    return out


def fingerprints_match(a: str, b: str) -> bool:
    """True if two fingerprint JSON blobs agree once volatile keys are dropped.

    Unparseable input returns False -- "I cannot tell" must mean "push it", never
    "skip it", or a malformed fingerprint would silently suppress a real drift PR.
    """
    try:
        return strip_ignored(json.loads(a)) == strip_ignored(json.loads(b))
    except (ValueError, AttributeError, TypeError):
        return False


# One full regeneration cycle. A PR still open after this was not acted on.
ESCALATE_AFTER = timedelta(days=14)


def should_escalate(*, pr_created_at: str, now: str) -> bool:
    """True if an open drift PR has outlived one full cycle.

    Returns False on any unparseable input: escalation is a notification, and a parse
    failure must not turn into repeated comments on the PR.
    """
    def _parse(s: str) -> datetime | None:
        try:
            return datetime.fromisoformat(s.replace("Z", "+00:00")).astimezone(timezone.utc)
        except (ValueError, AttributeError, TypeError):
            return None

    created, current = _parse(pr_created_at), _parse(now)
    if created is None or current is None:
        return False
    return (current - created) >= ESCALATE_AFTER


def branch_name_for_signal(datasets: list[str], fingerprint_path: str = "") -> str:
    """Branch name for a drift signal covering one or more datasets.

    A single dataset keeps its historical name exactly (``drift/<provider>-<dataset>``),
    so existing branches and their open PRs are still matched.

    A group takes ``drift/<provider>``, plus the baseline's distinguishing suffix when it
    has one. The suffix matters: the branch must identify the SIGNAL, and a provider can
    own more than one. Multi-dataset providers suffix their baselines per dataset
    (``drift_fingerprint_samples.json`` alongside ``drift_fingerprint.json``), so keying
    on the provider alone would collapse two independent signals onto one branch, where
    the second would overwrite the first's commit and rewrite its PR body. Deriving from
    the path rather than from the member list also keeps the name stable across runs, and
    does not privilege whichever member the manifest happens to list first.
    """
    if len(datasets) == 1:
        return branch_name_for(datasets[0])

    providers = {split_dataset(d)[0] for d in datasets}
    if len(providers) != 1:
        # A baseline shared ACROSS providers should be impossible -- the path lives
        # inside one plugin directory -- but stay deterministic rather than arbitrary.
        return branch_name_for(sorted(datasets)[0])

    base = "drift/" + providers.pop()
    stem = Path(fingerprint_path).stem if fingerprint_path else ""
    suffix = stem[len("drift_fingerprint"):].strip("_-") if stem.startswith("drift_fingerprint") else stem
    return f"{base}-{suffix}" if suffix else base


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
    grouped: dict[tuple, dict] = {}
    out: list[dict] = []
    for entry in drifted:
        path = entry.get("fingerprint_path")
        if not path:
            out.append({**entry, "datasets": [entry.get("dataset_name")]})
            continue
        # Key on (baseline, probe), exactly as `run_drift_checks` does. Grouping on the
        # baseline alone would re-merge what the runner deliberately kept apart: two
        # datasets sharing a baseline while declaring DIFFERENT probes is a manifest
        # error, and merging them would open a single PR whose regeneration covers only
        # `entry["dataset_name"]`, leaving the other dataset's drift silently
        # unaddressed. A report without probe_ref (older shape) falls back to the path,
        # which is the pre-existing behaviour rather than a new risk.
        key = (path, entry.get("probe_ref"))
        existing = grouped.get(key)
        if existing is None:
            merged = {**entry, "datasets": [entry.get("dataset_name")]}
            grouped[key] = merged
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


# Fallback reviewer when a plugin declares no `maintainers:` in its plugin.yaml.
# No manifest currently declares that field, so without this the --assignee flag is
# never emitted and drift PRs land on nobody's list -- half of why five of them sat
# unreviewed for three days. Set via DRIFT_DEFAULT_ASSIGNEE in the workflow env; a
# per-plugin `maintainers:` still wins wherever one is declared.
DEFAULT_ASSIGNEE = os.environ.get("DRIFT_DEFAULT_ASSIGNEE", "").strip()


def resolve_assignees(maintainers: list[str]) -> list[str]:
    """Declared maintainers if any, else the configured fallback, else nothing.

    Every returned handle is validated against ``_GITHUB_HANDLE``: `gh pr create`
    errors on a malformed or empty --assignee, and a junk env value must degrade to
    "no assignee" rather than fail the whole run.
    """
    valid = [m for m in maintainers if _GITHUB_HANDLE.match(m)]
    if valid:
        return valid
    if DEFAULT_ASSIGNEE and _GITHUB_HANDLE.match(DEFAULT_ASSIGNEE):
        return [DEFAULT_ASSIGNEE]
    return []


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
    # Lead with the verdict table (V5) -- built from what this function already has
    # (`dataset`, `diff`), not from the full report entry, which build_pr_body's
    # signature does not carry.
    lines.append(classification_table([{"dataset_name": dataset, "diff": diff}]))
    lines.append("")
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


def _discard_staged_fingerprints(*, dry_run: bool) -> None:
    """Return index AND working tree for ``hvantk/skills`` to HEAD.

    Every early return out of `handle_drifted` after `git add` must call this. The
    regeneration is per-dataset but the checkout is not isolated: the next iteration
    does `git checkout -B <next-branch> origin/<base>`, which leaves both the index and
    the working tree untouched. A leftover staged fingerprint would therefore be picked
    up by the next dataset's `git add hvantk/skills` and committed onto ITS branch --
    contaminating an unrelated PR, and masking that dataset's own skip check because
    the diff is no longer empty.

    `git checkout HEAD -- <path>` rather than `git reset`: reset alone unstages but
    leaves the modified file in the working tree, where the next `git add` re-stages it.
    """
    if dry_run:
        print("[dry-run] (would discard staged fingerprint changes)")
        return
    result = subprocess.run(
        ["git", "checkout", "HEAD", "--", "hvantk/skills"],
        check=False, text=True, capture_output=True,
    )
    if result.returncode != 0:
        # `check=False` keeps one failed cleanup from aborting the whole run, but
        # swallowing the output would hide the exact state this function exists to
        # prevent: the fingerprint stays staged and the next dataset commits it onto
        # ITS branch. Surface it so the job log names the contaminating run.
        print(
            "  WARNING: could not discard staged fingerprints; the next dataset may "
            f"commit them onto its branch: {(result.stderr or '').strip()}",
            file=sys.stderr,
        )


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
        # `errors="replace"` rather than a bare read_text(): a staged file that is not
        # valid UTF-8 raises UnicodeDecodeError, which is not a CalledProcessError, so
        # `main` would not catch it and every remaining dataset would be skipped. And
        # `paths` is every staged path, not only JSON under hvantk/skills, so that
        # depends on runner state rather than on this script's own staging. A mangled
        # decode can only make the comparison unequal, which pushes -- the safe answer.
        target = REPO_ROOT / path
        current = target.read_text(errors="replace") if target.exists() else None
        # errors="replace" here too: text=True decodes strictly, so a blob that is not
        # valid UTF-8 would raise UnicodeDecodeError out of subprocess itself.
        previous = subprocess.run(
            ["git", "show", f"FETCH_HEAD:{path}"],
            check=False, text=True, errors="replace", capture_output=True,
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


def _pr_created_at(pr_number: str, *, dry_run: bool) -> str:
    """The PR's ``createdAt`` timestamp, or ``""`` if unknown.

    A value-returning read, so -- like ``pr_exists_for_branch`` / ``remote_branch_exists``
    above -- it gates its own ``dry_run`` short-circuit and calls ``subprocess.run``
    directly rather than going through ``_run``. ``_run`` is for the fire-and-forget
    action commands (commit, push, ``gh pr create``/``edit``/``comment``), each already
    gated uniformly by its own caller; routing a query through it here would make this
    read visible to every test that captures ``_run`` calls to assert "no PR command
    was issued" on a path that only ever meant "no PR was created or edited" -- e.g.
    the skip path's own no-churn guarantee.
    """
    if dry_run:
        return ""
    result = subprocess.run(
        ["gh", "pr", "view", pr_number, "--json", "createdAt", "--jq", ".createdAt"],
        check=False, text=True, capture_output=True,
    )
    return (result.stdout or "").strip()


def maybe_escalate(pr_number: str, *, dry_run: bool) -> None:
    """Comment once on a drift PR that has outlived a full cycle.

    Deliberately a comment on the EXISTING PR, never a new PR: #268 fixed the inverse
    failure where nine PRs were force-pushed every morning, ~63 notification events a
    week carrying no new information.
    """
    created_at = _pr_created_at(pr_number, dry_run=dry_run)
    if not created_at:
        return
    if not should_escalate(
        pr_created_at=created_at, now=datetime.now(timezone.utc).isoformat()
    ):
        return
    _run(
        ["gh", "pr", "comment", pr_number, "--body",
         "This drift PR has been open for a full regeneration cycle (14 days). "
         "Upstream is still drifted and the baseline here is still unmerged."],
        dry_run=dry_run, check=False,
    )


def load_ledger() -> dict:
    """Read the ledger. A missing or corrupt file yields {} rather than raising --
    a broken ledger must not block a drift PR, it just starts recording afresh.
    """
    if not LEDGER_PATH.is_file():
        return {}
    try:
        return json.loads(LEDGER_PATH.read_text()) or {}
    except ValueError:
        return {}


def record_in_ledger(entries: list[dict], *, pr_ref: str, dry_run: bool) -> None:
    """Update the ledger for every dataset in ``entries`` and write it back.

    Callers must invoke this only once they have already decided a commit is
    happening -- i.e. AFTER both anti-churn guards in handle_drifted /
    handle_routine_batch (the `git diff --cached --quiet` emptiness check and
    `branch_needs_update`) have passed, and the caller must stage the ledger file
    itself in a separate `git add` right after calling this, rather than folding it
    into the earlier `git add hvantk/skills`.
    That earlier add is exactly what both guards inspect: `last_upstream_change` is
    `datetime.now(...)` at call time, so it is a different value on literally every
    invocation. Staging the ledger before either guard runs would make
    `git diff --cached --quiet` non-empty even when the fingerprint content did not
    change, and would make `branch_needs_update`'s per-path `fingerprints_match` check
    see the ledger's own timestamp move on every run -- both report "needs update"
    unconditionally, silently restoring the exact daily re-push/notification churn
    those guards exist to prevent (see the module docstring and the comment beside
    `branch_needs_update`). Calling it here, right before the commit, still gets the
    ledger change into the SAME commit as the fingerprints it describes -- just via
    its own `git add` rather than the earlier one.
    """
    if dry_run:
        print(f"  [dry-run] would record {len(entries)} ledger entries")
        return
    now = datetime.now(timezone.utc).isoformat()
    ledger = load_ledger()
    for entry in entries:
        ledger = ledger_update(
            ledger,
            dataset=entry["dataset_name"],
            diff=entry.get("diff"),
            pr_ref=pr_ref,
            now=now,
        )
    LEDGER_PATH.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")


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
    branch = branch_name_for_signal(covers, entry.get("fingerprint_path") or "")
    skill_md = find_skill_md(provider, dataset_short)
    maintainers = read_maintainers(provider)
    # Validation (bare GitHub handles only) and the default-assignee fallback both
    # live in resolve_assignees now; see its docstring.
    assignees = resolve_assignees(maintainers)
    diff = entry.get("diff") or {}

    print(f"\n=== Drifted: {dataset} -> branch {branch} ===")

    # 1) branch off base
    _run(["git", "fetch", "origin", base_branch], dry_run=dry_run, check=False)
    _run(["git", "checkout", "-B", branch, f"origin/{base_branch}"], dry_run=dry_run)

    # 2) regenerate fingerprint via the CLI we already ship
    try:
        _run(
            ["python", "-m", "hvantk.hvantk", "drift", "--regenerate", dataset],
            dry_run=dry_run,
        )
    except Exception:
        # A failed regenerate can still leave THIS dataset's fingerprint dirty in the
        # working tree (partially written, or left over from a prior attempt). `main()`
        # catches the CalledProcessError this raises, logs it, and moves on to the NEXT
        # schema-loop entry -- whose `git checkout -B` does not touch an already-dirty
        # working tree, and whose `git add hvantk/skills` would stage (and then commit)
        # this leftover onto a PR that has nothing to do with it. Discard before the
        # failure is allowed to propagate. Mirrors the same guard around
        # handle_routine_batch's regenerate loop, one dataset at a time instead of N.
        _discard_staged_fingerprints(dry_run=dry_run)
        raise

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
            _discard_staged_fingerprints(dry_run=dry_run)
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
        # Only skip when an open PR actually exists to be left alone. A branch can
        # outlive its PR -- closing a PR does not delete the head branch, and a run
        # whose push succeeded while `gh pr create` failed leaves a branch with no PR
        # at all (that exact failure is why this script exits nonzero on gh errors; see
        # the module docstring). In either case the branch content matches, so a
        # content-only check would skip forever and the dataset's drift would never be
        # surfaced again -- the precise outcome `branch_needs_update` promises cannot
        # happen. Re-opening a PR for an existing branch is cheap; silence is not.
        # `existing_pr` is captured via walrus so the escalation call below can reuse
        # it, while preserving the original short-circuit: `pr_exists_for_branch`
        # (a `gh pr list` call) still runs only when `branch_needs_update` is False,
        # exactly as before this change.
        if not branch_needs_update(branch, dry_run=dry_run) and (
            existing_pr := pr_exists_for_branch(branch, dry_run=dry_run)
        ):
            print(
                f"  branch {branch} already proposes this fingerprint "
                f"(only volatile keys differ); leaving it untouched."
            )
            # The regenerated fingerprint is still staged at this point. Leaving it
            # there would carry THIS dataset's baseline into the NEXT dataset's branch:
            # `git checkout -B` does not clear the index, so the next iteration's
            # `git add hvantk/skills` would stage both, and the next PR would commit a
            # bump it has nothing to do with. It would also defeat this very skip for
            # every dataset processed after a skipped one.
            _discard_staged_fingerprints(dry_run=dry_run)
            # This IS "this PR is still sitting there unmerged": say so on the PR
            # itself once it has outlived a full cycle, rather than opening another.
            maybe_escalate(existing_pr, dry_run=dry_run)
            _summary_line(
                step_summary,
                f"- DRIFT (unchanged): `{dataset}` -> branch `{branch}` still open; "
                f"nothing new to push",
            )
            return

    # Record the accepted change in the rebuild ledger, and stage it on its own --
    # only now, after both guards above have passed, so its ever-moving
    # `last_upstream_change` timestamp cannot defeat either of them. See
    # record_in_ledger's docstring.
    record_in_ledger([entry], pr_ref=branch, dry_run=dry_run)
    _run(["git", "add", "hvantk/resources/drift_ledger.json"], dry_run=dry_run)

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
        print("  creating new PR")
        _run(
            pr_create_argv(
                base_branch=base_branch, branch=branch, title=title, body=body,
                risk="schema", assignees=assignees,
            ),
            dry_run=dry_run,
        )

    if dry_run:
        print("\n--- PR body preview ---")
        print(body)

    _summary_line(
        step_summary,
        f"- DRIFT: `{dataset}` -> branch `{branch}` (PR opened/updated)",
    )


def handle_routine_batch(
    entries: list[dict],
    *,
    base_branch: str,
    dry_run: bool,
    step_summary: Path | None,
) -> None:
    """Regenerate every routine fingerprint onto ONE branch and open ONE PR.

    Mirrors handle_drifted's sequence, but the regenerate step loops over datasets
    before a single commit, so N routine datasets cost one PR instead of N.
    """
    if not entries:
        return

    datasets = [e["dataset_name"] for e in entries]
    branch = ROUTINE_BRANCH
    print(f"\n=== Routine batch ({len(datasets)}) -> branch {branch} ===")
    for name in datasets:
        print(f"      {name}")

    _run(["git", "fetch", "origin", base_branch], dry_run=dry_run, check=False)
    _run(["git", "checkout", "-B", branch, f"origin/{base_branch}"], dry_run=dry_run)

    try:
        for name in datasets:
            _run(
                ["python", "-m", "hvantk.hvantk", "drift", "--regenerate", name],
                dry_run=dry_run,
            )
    except Exception:
        # A mid-loop failure leaves every dataset regenerated BEFORE the failing one
        # sitting dirty in the working tree (staging happens only once, after the whole
        # loop finishes). `main()` catches the CalledProcessError this raises, logs it,
        # and continues into the schema loop -- whose `git checkout -B` does not touch
        # an already-dirty working tree, and whose `git add hvantk/skills` would stage
        # (and then commit) these leftovers onto an unrelated PR: a schema-change PR,
        # invisible in its diff and its ledger entry. Discard before the failure is
        # allowed to propagate.
        _discard_staged_fingerprints(dry_run=dry_run)
        raise

    _run(["git", "add", "hvantk/skills"], dry_run=dry_run)

    # No-op commits should not fail the job -- check first, exactly like
    # handle_drifted. Reachable in practice: a re-run after someone already merged the
    # fix manually, or a race where this branch already carries every regenerated
    # fingerprint byte-for-byte. Routed through `_run` (rather than a direct
    # subprocess.run call, unlike handle_drifted) so dry_run short-circuits it the same
    # way as every other command here -- calling it unconditionally would make `_run`
    # fabricate a returncode of 0 under --dry-run, which would misreport "nothing to
    # commit" on every dry run and silently swallow the demonstration path.
    if dry_run:
        print("[dry-run] (skipping diff/commit emptiness check)")
    else:
        status = _run(
            ["git", "diff", "--cached", "--quiet"], dry_run=dry_run, check=False
        )
        if status.returncode == 0:
            print(
                f"  no fingerprint changes to commit for the routine batch of "
                f"{len(datasets)}; drift may have already been addressed. "
                "Skipping PR."
            )
            # `git diff --cached --quiet` only proves the INDEX matches HEAD -- it
            # says nothing about whether the working tree has a stray modification
            # sitting outside the index (a .gitignore quirk, a partial `add`, or any
            # other way a regenerated file could escape being staged). Discarding is
            # cheap and resets both the index AND the working tree, so call it
            # unconditionally rather than assume "nothing staged" implies "nothing to
            # clean up" -- `git checkout -B` carries the working tree forward across
            # branches regardless, and this function returns right into the schema
            # loop's first `git checkout -B`, which would inherit anything left behind.
            _discard_staged_fingerprints(dry_run=dry_run)
            _summary_line(
                step_summary,
                f"- DRIFT (unchanged): routine batch of {len(datasets)} had nothing "
                "to commit; drift may have already been addressed",
            )
            return

    # See the matching comment in handle_drifted: walrus preserves the original
    # short-circuit so `pr_exists_for_branch` still runs only when needed, while making
    # the PR number available to the escalation call below.
    if not branch_needs_update(branch, dry_run=dry_run) and (
        existing_pr := pr_exists_for_branch(branch, dry_run=dry_run)
    ):
        print(f"  branch {branch} already proposes these fingerprints; leaving it.")
        _discard_staged_fingerprints(dry_run=dry_run)
        # This IS "this PR is still sitting there unmerged": say so on the PR itself
        # once it has outlived a full cycle, rather than opening another.
        maybe_escalate(existing_pr, dry_run=dry_run)
        _summary_line(
            step_summary,
            f"- DRIFT (unchanged): routine batch of {len(datasets)} still open",
        )
        return

    # See the matching comment in handle_drifted: the ledger is written and staged in
    # its own `git add` only after both anti-churn guards above have passed, never
    # before, or its ever-changing `last_upstream_change` timestamp would defeat them.
    record_in_ledger(entries, pr_ref=branch, dry_run=dry_run)
    _run(["git", "add", "hvantk/resources/drift_ledger.json"], dry_run=dry_run)

    _run(
        ["git", "commit", "-m",
         f"chore(drift): refresh {len(datasets)} snapshots\n\n"
         + "\n".join(f"- {n}" for n in datasets)],
        dry_run=dry_run,
    )
    _run(
        ["git", "push", "--force-with-lease", "--set-upstream", "origin", branch],
        dry_run=dry_run,
    )

    title = f"chore(drift): refresh {len(datasets)} snapshots"
    body = build_batch_pr_body(entries)
    # Collect maintainers across every dataset in the batch, then resolve ONCE on the
    # combined list -- not per dataset, or a fallback would land on the batch as many
    # times as it has entries with no declared maintainer of their own.
    maintainers = sorted({
        m
        for e in entries
        for m in read_maintainers(split_dataset(e["dataset_name"])[0])
    })
    assignees = resolve_assignees(maintainers)

    existing = pr_exists_for_branch(branch, dry_run=dry_run)
    if existing:
        _run(["gh", "pr", "edit", existing, "--title", title, "--body", body],
             dry_run=dry_run)
    else:
        _run(
            pr_create_argv(
                base_branch=base_branch, branch=branch, title=title, body=body,
                risk="routine", assignees=assignees,
            ),
            dry_run=dry_run,
        )

    _summary_line(
        step_summary,
        f"- DRIFT: routine batch of {len(datasets)} -> branch `{branch}`",
    )


def pr_create_argv(
    *,
    base_branch: str,
    branch: str,
    title: str,
    body: str,
    risk: str,
    assignees: list[str],
) -> list[str]:
    """Build the ``gh pr create`` argv.

    Deliberately NOT a draft (V1): a draft cannot be merged and is filtered out of
    review queues and notification defaults, which is how five drift PRs sat
    unreviewed for three days. Extracted as a pure function so the flags are
    testable without invoking gh.

    Never emits ``--auto``: nothing in this pipeline may auto-merge.
    """
    argv = [
        "gh", "pr", "create",
        "--base", base_branch,
        "--head", branch,
        "--title", title,
        "--body", body,
        "--label", f"drift:{risk}",
    ]
    if assignees:
        argv += ["--assignee", ",".join(assignees)]
    return argv


def classification_table(entries: list[dict]) -> str:
    """Markdown table summarising what moved per dataset, and the verdict.

    Leads the PR body (V5) so a reviewer sees the judgement before the raw JSON
    diffs. The old body opened with per-dataset JSON, which is why five PRs were
    indistinguishable at a glance.
    """
    rows = [
        "| Dataset | What moved | Verdict |",
        "| --- | --- | --- |",
    ]
    for entry in entries:
        diff = entry.get("diff") or {}
        changed = diff.get("changed") or {}
        signal = ", ".join(f"`{k}`" for k in sorted(changed)) or "—"
        verdict = (
            "routine — schema unchanged"
            if classify_risk(diff) == "routine"
            else "**SCHEMA CHANGE** — check `builder.py`"
        )
        rows.append(f"| `{entry['dataset_name']}` | {signal} | {verdict} |")
    return "\n".join(rows)


def build_batch_pr_body(entries: list[dict]) -> str:
    """PR body for the routine batch: verdict table first, raw diffs collapsed below."""
    parts = [
        f"Automated drift detection found upstream changes for "
        f"{len(entries)} dataset(s). The schema signal is unchanged for every one "
        f"below — only content and/or version moved.",
        "",
        classification_table(entries),
        "",
        "<details><summary>Raw fingerprint diffs</summary>",
        "",
    ]
    for entry in entries:
        parts += [
            f"### `{entry['dataset_name']}`",
            "",
            "```json",
            json.dumps(entry.get("diff") or {}, indent=2, sort_keys=True),
            "```",
            "",
        ]
    parts.append("</details>")
    return "\n".join(parts)


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

    # Routine drift (schema signal unchanged) batches onto ONE branch/PR; schema-risk
    # drift keeps the existing one-branch-per-signal path so a human looks at it on
    # its own. See classify_risk / partition_by_risk.
    routine, schema = partition_by_risk(groups)
    if routine:
        print(
            f"  {len(routine)} of {len(groups)} signal(s) are routine (schema "
            f"unchanged); batching them onto one branch: {ROUTINE_BRANCH}."
        )

    failed: list[str] = []
    try:
        handle_routine_batch(
            routine,
            base_branch=args.base_branch,
            dry_run=args.dry_run,
            step_summary=step_summary,
        )
    except subprocess.CalledProcessError as exc:
        # Same contract as the per-dataset loop below: a git/gh failure here must not
        # be swallowed, or the batch silently stops opening PRs while the job still
        # reports success -- the exact regression this script exists to prevent (see
        # the module docstring).
        covered = [str(e.get("dataset_name")) for e in routine]
        failed.extend(covered)
        label = ", ".join(covered)
        print(
            f"error handling routine batch ({label}): "
            f"{exc.cmd} exited {exc.returncode}\n"
            f"stdout: {exc.stdout}\nstderr: {exc.stderr}",
            file=sys.stderr,
        )
        _summary_line(
            step_summary,
            f"- ERROR: routine batch (`{label}`) -- {exc.cmd[0]} exited {exc.returncode}",
        )

    for entry in schema:
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
            #
            # Name every dataset the signal covers, not just the anchor. A grouped
            # entry regenerates one baseline on behalf of several datasets, so
            # reporting `dataset_name` alone would show the operator one dataset when
            # several are left unaddressed.
            covered = [str(d) for d in (entry.get("datasets") or [entry.get("dataset_name")])]
            failed.extend(covered)
            label = ", ".join(covered)
            print(
                f"error handling {label}: "
                f"{exc.cmd} exited {exc.returncode}\n"
                f"stdout: {exc.stdout}\nstderr: {exc.stderr}",
                file=sys.stderr,
            )
            _summary_line(
                step_summary,
                f"- ERROR: `{label}` -- {exc.cmd[0]} exited {exc.returncode}",
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
