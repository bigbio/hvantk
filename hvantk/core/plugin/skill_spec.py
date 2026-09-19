"""The per-resource ``SKILL.md`` contract, as code rather than as prose.

``hvantk/skills/_conventions/SKILL.md`` s 2 says every per-resource ``SKILL.md``
**MUST** carry nine sections, in order, with exact headings. Nothing checked that,
so the tree and the contract drifted apart silently: 9 of 23 provider specs
contained zero of the nine headings and no YAML frontmatter (issue #334).

This module is the single implementation of that check. Both consumers import it
rather than re-deriving the rules:

* ``hvantk/tests/test_plugin_skill_conformance.py`` -- enforces it on every run.
* ``hvantk plugins validate <manifest>`` -- reports it for one plugin.

Imports nothing heavier than ``yaml`` (already a hard dependency -- the plugin
loader parses every ``plugin.yaml`` with it) and never touches Hail, so the test
that uses it runs in the default ``pytest`` selection rather than only in the
``hail``-marked job.

Why frontmatter is required as well as the headings: the portable Agent Skills
convention that external harnesses use to discover a skill directory needs at
minimum ``name`` and ``description``. A spec without them is invisible to those
harnesses even when its prose is perfect.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

import yaml

#: The nine mandatory headings, verbatim and in order (``_conventions`` s 2).
REQUIRED_SECTIONS: tuple[str, ...] = (
    "## 1. Status & scope",
    "## 2. Source identity",
    "## 3. Backend choice + reasoning",
    "## 4. Raw format & gotchas",
    "## 5. Output contract",
    "## 6. hvantk integration points",
    "## 7. Workflow steps",
    "## 8. Update playbook",
    "## 9. Validation contract",
)

#: Frontmatter keys every per-resource spec must declare. ``name`` and
#: ``description`` are what an external skill harness reads to register the
#: directory at all; the rest are hvantk's own routing metadata.
REQUIRED_FRONTMATTER: tuple[str, ...] = ("name", "description")

#: ``_conventions/SKILL.md`` IS the contract and has its own structure, so it is
#: not measured against itself.
EXEMPT_DIRS: frozenset[str] = frozenset({"_conventions"})


def iter_skill_specs(skills_dir: Path) -> Iterator[Path]:
    """Yield every per-resource ``SKILL.md`` under ``skills_dir``, sorted.

    Covers both layouts from ``_conventions`` s 1: ``<provider>/SKILL.md`` for a
    single-dataset provider and ``<provider>/<dataset>/SKILL.md`` for a
    multi-dataset one. Exempt directories are skipped.
    """
    seen: set[Path] = set()
    for pattern in ("*/SKILL.md", "*/*/SKILL.md"):
        for path in sorted(skills_dir.glob(pattern)):
            if EXEMPT_DIRS.intersection(path.relative_to(skills_dir).parts):
                continue
            if path not in seen:
                seen.add(path)
                yield path


def parse_frontmatter(text: str) -> dict[str, object] | None:
    """Return the leading YAML frontmatter as a mapping, or ``None`` if absent/invalid.

    Parsed with ``yaml.safe_load`` rather than by hand. An earlier hand-rolled
    version split on the first ``:`` of each line, which accepted things an Agent
    Skills harness could not load: ``name: [`` came back as the string ``"["``, and
    any file merely *starting with* the characters ``---`` (e.g. ``---not-frontmatter``)
    was treated as opening a block. The opening line must now equal ``---`` exactly.

    PyYAML resolves duplicate keys silently, last-one-wins, so a spec declaring
    ``name`` twice would load clean and quietly drop the first value. That is checked
    separately here rather than delegated, since it cannot be delegated.

    Returns ``None`` when there is no well-formed frontmatter block at all;
    ``check_skill_spec`` turns that into the user-facing message.
    """
    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        return None
    try:
        end = next(
            i for i, line in enumerate(lines[1:], start=1) if line.strip() == "---"
        )
    except StopIteration:
        return None
    block = "\n".join(lines[1:end])
    try:
        loaded = yaml.safe_load(block)
    except yaml.YAMLError:
        return None
    if not isinstance(loaded, dict):
        return None
    return loaded


def duplicate_frontmatter_keys(text: str) -> list[str]:
    """Top-level keys declared more than once in the frontmatter block.

    Separate from :func:`parse_frontmatter` because ``yaml.safe_load`` cannot report
    this -- it silently keeps the last value, so the duplicate has to be found before
    parsing collapses it.
    """
    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        return []
    try:
        end = next(
            i for i, line in enumerate(lines[1:], start=1) if line.strip() == "---"
        )
    except StopIteration:
        return []
    seen: set[str] = set()
    dupes: list[str] = []
    for line in lines[1:end]:
        # Top-level keys only: an indented line belongs to a nested structure.
        if not line.strip() or line.startswith((" ", "\t", "#")):
            continue
        key, sep, _ = line.partition(":")
        if not sep:
            continue
        key = key.strip()
        if key in seen and key not in dupes:
            dupes.append(key)
        seen.add(key)
    return dupes


def _body_heading_lines(text: str) -> list[str]:
    """Heading lines of the document body, excluding fenced code blocks.

    Both exclusions matter. Matching a heading as a line *prefix* accepted
    ``## 1. Status & scope (draft)`` as though it were the mandated heading, and
    matching anywhere in the file accepted a heading that only ever appeared quoted
    inside a ``\u0060\u0060\u0060`` fence -- including, awkwardly, the fence in this repo's own
    conventions document.
    """
    out: list[str] = []
    in_fence = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("```") or stripped.startswith("~~~"):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        if stripped.startswith("#"):
            out.append(stripped)
    return out


def check_skill_spec(path: Path) -> list[str]:
    """Return a list of contract violations for one ``SKILL.md``; empty means conforming.

    Checks, in the order a reader would hit them:

    1. YAML frontmatter is present, parses as YAML, declares no duplicate top-level
       key, and carries every key in ``REQUIRED_FRONTMATTER``.
    2. Every heading in ``REQUIRED_SECTIONS`` appears as a complete body heading
       line -- not as a prefix of a longer one, and not inside a fenced code block.
    3. Those headings appear in the declared order.

    Ordering is checked separately from presence so the message says which of the
    two failed -- a spec with all nine headings shuffled is a different (and much
    smaller) fix than one missing six of them.
    """
    problems: list[str] = []
    try:
        text = path.read_text()
    except OSError as exc:  # pragma: no cover - unreadable file is its own failure
        return [f"could not read: {exc}"]

    frontmatter = parse_frontmatter(text)
    if frontmatter is None:
        problems.append(
            "missing or unparseable YAML frontmatter (the file must open with a line "
            "that is exactly '---', and the block must be valid flat YAML)"
        )
    else:
        absent = [k for k in REQUIRED_FRONTMATTER if not frontmatter.get(k)]
        if absent:
            problems.append(f"frontmatter missing or empty: {', '.join(absent)}")
    dupes = duplicate_frontmatter_keys(text)
    if dupes:
        # Worth its own message: YAML keeps the LAST value, so the author's first
        # (probably intended) one is discarded without any error anywhere.
        problems.append(
            f"frontmatter declares duplicate key(s): {', '.join(dupes)} "
            "(YAML silently keeps the last value)"
        )

    headings = _body_heading_lines(text)
    positions: list[tuple[str, int]] = []
    for section in REQUIRED_SECTIONS:
        # Whole-line equality. A prefix match accepted '## 1. Status & scope (draft)'.
        try:
            idx = headings.index(section)
        except ValueError:
            problems.append(f"missing required section: {section!r}")
        else:
            positions.append((section, idx))

    if len(positions) == len(REQUIRED_SECTIONS):
        ordered = [s for s, _ in sorted(positions, key=lambda pair: pair[1])]
        if ordered != list(REQUIRED_SECTIONS):
            problems.append(
                "required sections are out of order; expected "
                + " -> ".join(s.split(". ", 1)[-1] for s in REQUIRED_SECTIONS)
                + "; found "
                + " -> ".join(s.split(". ", 1)[-1] for s in ordered)
            )

    return problems
