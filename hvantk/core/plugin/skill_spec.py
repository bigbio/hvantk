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

import re
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

#: The line every spec opens with, directing a reader to the shared contract before the
#: provider-specific detail. 21 of 23 already carried it in two wordings; #356 settled on
#: the majority one and made it checkable, because a spec read in isolation otherwise
#: looks self-contained when it is not. Checked as the first line ``_first_body_line``
#: returns, not as a substring anywhere in the file (#364): a whole-file ``in text`` test
#: passed a spec with the line moved to the end, or quoted only inside a code fence.
CONVENTIONS_PREAMBLE: str = (
    "Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its "
    "repository map, helpers, keying conventions, builder pattern, and validation "
    "contract."
)

#: What s 6 must let a reader reach without opening plugin.yaml. Each entry is
#: (label, substrings that satisfy it). Before #356, seven specs named no drift probe
#: and five named no tests -- so the section existed while the wiring it documents did
#: not, which is the same failure as a heading with nothing under it.
SECTION_6_REFERENCES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("the builder", (r"\bbuilder\b", r"build_")),
    (
        "the drift probe",
        (r"drift_probe", r"drift probe", r"fetch_fingerprint", r"hvantk drift"),
    ),
    # Word-bounded: a bare "test" substring is satisfied by "latest", which these
    # specs say constantly because that is what drift probes fingerprint.
    ("its tests", (r"\btests?\b", r"\bpytest\b")),
)

#: The five test-artifact keys s 9 must name, spelled as the manifest's `tests:` keys.
#: `command`, NOT `test_command`: the schema sets additionalProperties false, so a
#: manifest authored from the wrong spelling fails validation at load. 22 of 23 specs
#: used the wrong one until #350 (24 occurrences).
SECTION_9_ARTIFACT_KEYS: tuple[str, ...] = (
    "fixture",
    "schema_snapshot",
    "row_snapshot",
    "drift_fingerprint",
    "command",
)

#: Tokens that open a stand-in body rather than real content. Checked as a LEAD --
#: the body need only start this way, not equal it exactly -- because the rot #364
#: item 8 found was not only a bare ``TODO``/``TBD``/``...``/``n/a`` but also a longer
#: sentence that opens by admitting it is one (``TODO. This section will be fleshed
#: out once ...``). Before this, only s 6 and s 9 had any content rule at all, so
#: ss 1-5, 7 and 8 could be emptied outright and the checker kept passing; and s 6 / s 9's
#: own rules checked for keywords rather than substance, so a body that only *named* the
#: required words satisfied them too -- ``TODO. builder, drift probe, tests.`` passes
#: ``SECTION_6_REFERENCES``, and a sentence merely listing the five ``SECTION_9_ARTIFACT_KEYS``
#: passes that check, without either being an actual section.
_PLACEHOLDER_LEADS: tuple[str, ...] = ("todo", "tbd", "...")

#: Tokens that are a stand-in only when they are the WHOLE body. ``n/a`` is also how a
#: real answer opens -- "N/A for this plugin: the catalog is static" explains why the
#: section has nothing to add -- so only a bare ``n/a`` (punctuation aside) is rejected.
_PLACEHOLDER_WHOLE: tuple[str, ...] = ("n/a",)


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
            parts = path.relative_to(skills_dir).parts
            if EXEMPT_DIRS.intersection(parts):
                continue
            # Track the loader, which skips every `_`-prefixed directory
            # (`load_from_skills_root`). Without this, adding `_hooks/SKILL.md` would
            # fail the conformance test for a directory that is not a plugin and that
            # the registry never reads.
            if any(part.startswith("_") for part in parts[:-1]):
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


def _first_body_line(text: str) -> str | None:
    """The first line of prose in the body: after the frontmatter, ignoring blank lines,
    headings and fenced code blocks. ``None`` when there is no such line (an empty body,
    or a frontmatter block that never closes).

    The preamble is checked HERE rather than by substring over the whole file: ``in text``
    accepted a spec with the line deleted from the top and pasted at the bottom, or present
    only inside a fence, while ``CONVENTIONS_PREAMBLE``'s own docstring says it is the line
    every spec *opens* with (#364).
    """
    lines = text.splitlines()
    start = 0
    if lines and lines[0].strip() == "---":
        try:
            start = (
                next(
                    i
                    for i, line in enumerate(lines[1:], start=1)
                    if line.strip() == "---"
                )
                + 1
            )
        except StopIteration:
            return None
    in_fence = False
    for line in lines[start:]:
        stripped = line.strip()
        if stripped.startswith("```") or stripped.startswith("~~~"):
            in_fence = not in_fence
            continue
        if in_fence or not stripped or stripped.startswith("#"):
            continue
        return stripped
    return None


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

    if _first_body_line(text) != CONVENTIONS_PREAMBLE:
        where = (
            " (it appears later in the file, or only inside a code fence)"
            if CONVENTIONS_PREAMBLE in text
            else ""
        )
        problems.append(
            "the body must open with the shared preamble line directing the reader to "
            f"_conventions/SKILL.md (see CONVENTIONS_PREAMBLE){where}"
        )

    problems.extend(_section_content_problems(text))

    return problems


def _section_body(text: str, section: str) -> str:
    """The text under ``section``, up to the next ``## `` heading or end of file.

    Fence-aware, to agree with ``_body_heading_lines``. A raw ``str.find`` would treat a
    heading quoted inside a ``` block as the real one -- and ``_conventions`` now ships a
    template containing exactly those headings in a fence. A spec that quoted the
    template above its real section 6 would then have its content checked against the
    template, which by construction mentions the builder, the drift probe and tests, so
    an empty real section would report conforming.
    """
    lines = text.splitlines()
    in_fence = False
    start = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("```") or stripped.startswith("~~~"):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        if start is None:
            if stripped == section:
                start = i + 1
        elif stripped.startswith("## "):
            return "\n".join(lines[start:i])
    return "\n".join(lines[start:]) if start is not None else ""


def _is_placeholder_body(body: str) -> bool:
    """True if the text under one heading is empty, or a stand-in rather than content.

    A leading list/quote marker (``- ``, ``* ``, ``> ``) and surrounding Markdown
    emphasis (backticks, asterisks, underscores) are stripped before matching, so
    ``- **TODO**`` and `` `TBD` `` are recognised the same as bare ``TODO`` / ``TBD``,
    and a lone ``-`` counts as empty. A typographic ellipsis (``…``) counts as ``...``.
    The character immediately following a lead token must not be alphanumeric, so the
    check fires on ``TODO``, ``TODO.`` and ``TODO: builder, drift probe, tests.``
    alike, but not on a word that merely starts with the same letters (``TBDone``).
    It does fire on a sentence whose first word is the token itself (``Todo lists
    are ...``): the price of checking a lead, paid because no real section opens that
    way.
    """
    normalized = re.sub(r"^[-*>\s]+", "", body.strip()).strip("`*_ \t")
    lowered = normalized.lower().replace("…", "...")
    if not lowered or lowered.rstrip(".:;!- \t") in _PLACEHOLDER_WHOLE:
        return True
    for token in _PLACEHOLDER_LEADS:
        if lowered.startswith(token):
            tail = lowered[len(token) : len(token) + 1]
            if not tail.isalnum():
                return True
    return False


def _section_content_problems(text: str) -> list[str]:
    """Check that every section carries a real body, and that s 6 / s 9 carry their
    specific required content, not just their headings.

    The per-section check runs against all nine ``REQUIRED_SECTIONS``, but only the
    ones actually present as a heading: a heading that is missing or malformed is
    already reported once, by the heading check in ``check_skill_spec``, and re-testing
    an empty string for a heading that was never found would only duplicate that
    message under a different name. For every heading that IS present, an empty or
    placeholder-led body (``_is_placeholder_body``) is rejected regardless of section --
    which is what makes this the fix for #364 item 8: before it, ss 1-5, 7 and 8 carried
    no content rule at all, and s 6 / s 9's own rules (below) checked for keywords
    rather than substance, so a body that only named the required words satisfied them
    without being a real section.

    The s 6 / s 9 keyword and artifact-key checks are deliberately a word-bounded regex
    search rather than a structural parse. These sections are prose with per-provider
    shape -- some use ``- **Label:**``, some ``- Label:`` -- and a parser strict enough
    to read the shape would fail on formatting instead of substance, which is the
    opposite of useful. The question asked here is only "can a reader get from this
    section to that thing", which a search answers honestly. The word boundaries
    matter: they are what keeps the s 9 ``command`` needle from matching inside
    ``test_command`` (``_`` is a word character, so there is no boundary between them).
    """
    problems: list[str] = []
    headings = _body_heading_lines(text)

    for section in REQUIRED_SECTIONS:
        if section not in headings:
            continue  # already reported as missing/malformed, above
        if _is_placeholder_body(_section_body(text, section)):
            problems.append(
                f"{section!r} has no real body -- it is empty, or only a "
                "TODO/TBD/.../n/a-style placeholder"
            )

    six = _section_body(text, "## 6. hvantk integration points").lower()
    for label, needles in SECTION_6_REFERENCES:
        if not any(re.search(n, six) for n in needles):
            problems.append(f"section 6 never references {label}")

    nine = _section_body(text, "## 9. Validation contract")
    absent = [k for k in SECTION_9_ARTIFACT_KEYS if not re.search(rf"\b{k}\b", nine)]
    if absent:
        problems.append(
            "section 9 does not name these test artifacts: " + ", ".join(absent)
        )
    if "test_command" in nine:
        problems.append(
            "section 9 says 'test_command'; the manifest key is 'command' "
            "(additionalProperties is false, so the wrong spelling fails at load)"
        )

    return problems
