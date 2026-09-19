"""The per-resource ``SKILL.md`` contract, as code rather than as prose.

``hvantk/skills/_conventions/SKILL.md`` s 2 says every per-resource ``SKILL.md``
**MUST** carry nine sections, in order, with exact headings. Nothing checked that,
so the tree and the contract drifted apart silently: 9 of 23 provider specs
contained zero of the nine headings and no YAML frontmatter (issue #334).

This module is the single implementation of that check. Both consumers import it
rather than re-deriving the rules:

* ``hvantk/tests/test_plugin_skill_conformance.py`` -- enforces it on every run.
* ``hvantk plugins validate <manifest>`` -- reports it for one plugin.

Deliberately dependency-free (stdlib only, no yaml, no Hail), so the test that
uses it runs in the default ``pytest`` selection rather than only in the
``hail``-marked job.

Why frontmatter is required as well as the headings: the portable Agent Skills
convention that external harnesses use to discover a skill directory needs at
minimum ``name`` and ``description``. A spec without them is invisible to those
harnesses even when its prose is perfect.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

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


def parse_frontmatter(text: str) -> dict[str, str] | None:
    """Return the leading YAML frontmatter as a flat mapping, or ``None``.

    Hand-rolled rather than ``yaml.safe_load`` to keep this module stdlib-only.
    The frontmatter block in these specs is flat ``key: value`` lines, so a
    parser that handles exactly that is sufficient -- and a spec whose
    frontmatter is NOT flat is itself a problem worth surfacing.
    """
    if not text.startswith("---"):
        return None
    lines = text.splitlines()
    try:
        end = next(
            i for i, line in enumerate(lines[1:], start=1) if line.strip() == "---"
        )
    except StopIteration:
        return None
    fields: dict[str, str] = {}
    for line in lines[1:end]:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        key, sep, value = line.partition(":")
        if sep:
            fields[key.strip()] = value.strip()
    return fields


def check_skill_spec(path: Path) -> list[str]:
    """Return a list of contract violations for one ``SKILL.md``; empty means conforming.

    Checks, in the order a reader would hit them:

    1. YAML frontmatter is present and carries every key in ``REQUIRED_FRONTMATTER``.
    2. Every heading in ``REQUIRED_SECTIONS`` is present, verbatim.
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
            "missing YAML frontmatter (a '---' delimited block must open the file)"
        )
    else:
        absent = [k for k in REQUIRED_FRONTMATTER if not frontmatter.get(k)]
        if absent:
            problems.append(f"frontmatter missing or empty: {', '.join(absent)}")

    positions: list[tuple[str, int]] = []
    for section in REQUIRED_SECTIONS:
        # Anchored to line starts so a heading quoted mid-paragraph does not count.
        idx = text.find(f"\n{section}")
        if idx == -1 and text.startswith(section):
            idx = 0
        if idx == -1:
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
