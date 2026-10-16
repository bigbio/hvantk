"""Enforce the workflow-skill contract (``skills/workflow-conventions/SKILL.md``) on every
``skills/<name>/SKILL.md``.

Sibling of ``test_plugin_skill_conformance.py``, which does the same for the resource
specs under ``hvantk/skills/``. The parsing helpers are imported from
``hvantk.core.plugin.skill_spec`` rather than re-derived, so the two contracts cannot
disagree on what a heading or a placeholder body is. Stdlib + yaml only: this runs in
the default ``pytest`` selection.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from hvantk.core.plugin.skill_spec import (
    _body_heading_lines,
    _is_placeholder_body,
    _section_body,
    duplicate_frontmatter_keys,
    parse_frontmatter,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
SKILLS_DIR = REPO_ROOT / "skills"

#: The conventions skill itself: not a workflow skill, so it is not a workflow-skill
#: catalog entry and is not checked against the six-section contract it defines.
CONVENTIONS_SKILL_NAME = "workflow-conventions"

WORKFLOW_SECTIONS: tuple[str, ...] = (
    "## 1. Scope and outcome",
    "## 2. Preflight",
    "## 3. Procedure",
    "## 4. Scale-up",
    "## 5. Traps",
    "## 6. Verification",
)
#: Top-level frontmatter keys, exactly: the two the Agent Skills layout requires, the
#: two Claude Code reads, and ``metadata`` for hvantk's own routing keys.
REQUIRED_FRONTMATTER: tuple[str, ...] = (
    "name",
    "description",
    "user-invocable",
    "argument-hint",
    "metadata",
)
METADATA_KEYS: tuple[str, ...] = ("stage", "status", "scale")
STAGES = {"substrate", "cohort", "annotation", "analysis"}
STATUSES = {"provisional", "stable"}
SCALES = {"laptop", "node", "cluster"}
EVAL_FILES: tuple[str, ...] = (
    "task.md",
    "stage.sh",
    "control.sh",
    "check.py",
    "reference.json",
)
_LOCAL_REF = re.compile(r"`((?:references|templates)/[^`\s]+)`")
#: A line that is nothing but an angle-bracket placeholder, optionally as a list item.
_PLACEHOLDER_LINE = re.compile(r"^\s*(?:[-*]\s+|\d+\.\s+)?<[^<>\n]+>\.?\s*$", re.M)


def iter_workflow_skills(root: Path = SKILLS_DIR):
    for path in sorted(root.glob("*/SKILL.md")):
        if path.parent.name.startswith("_") or path.parent.name == CONVENTIONS_SKILL_NAME:
            continue
        yield path


def check_workflow_skill(path: Path, *, require_files: bool = True) -> list[str]:
    """Contract violations for one workflow ``SKILL.md``; empty means conforming.

    ``require_files=False`` is for the copy-me template in ``workflow-conventions``: it still
    carries its ``<...>`` placeholders and has no folder of its own.
    """
    problems: list[str] = []
    text = path.read_text()
    fm = parse_frontmatter(text)
    if fm is None:
        return ["missing or unparseable YAML frontmatter (file must open with '---')"]
    dupes = duplicate_frontmatter_keys(text)
    if dupes:
        problems.append(f"frontmatter declares duplicate key(s): {', '.join(dupes)}")
    absent = [k for k in REQUIRED_FRONTMATTER if k not in fm or fm[k] in ("", None)]
    if absent:
        problems.append(f"frontmatter missing or empty: {', '.join(absent)}")
    extra = sorted(set(fm) - set(REQUIRED_FRONTMATTER))
    if extra:
        problems.append(
            f"frontmatter has keys outside the contract: {', '.join(extra)} "
            "(hvantk keys go under metadata)"
        )
    if fm.get("name") != path.parent.name:
        problems.append(
            f"frontmatter name {fm.get('name')!r} must equal the folder name "
            f"{path.parent.name!r} (Agent Skills rule)"
        )
    if not str(fm.get("description", "")).startswith("Use when"):
        problems.append("description must start with 'Use when' (it is the trigger)")
    if fm.get("user-invocable") is not True:
        problems.append("user-invocable must be true")
    meta = fm.get("metadata")
    if not isinstance(meta, dict):
        problems.append("metadata must be a mapping with stage, status and scale")
        meta = {}
    extra_meta = sorted(set(meta) - set(METADATA_KEYS))
    if extra_meta:
        problems.append(
            f"metadata has keys outside the contract: {', '.join(extra_meta)}"
        )
    if meta.get("stage") not in STAGES:
        problems.append(f"metadata.stage must be one of {sorted(STAGES)}")
    if meta.get("status") not in STATUSES:
        problems.append(f"metadata.status must be one of {sorted(STATUSES)}")
    if meta.get("scale") not in SCALES:
        problems.append(f"metadata.scale must be one of {sorted(SCALES)}")

    headings = _body_heading_lines(text)
    positions = []
    for section in WORKFLOW_SECTIONS:
        try:
            positions.append((section, headings.index(section)))
        except ValueError:
            problems.append(f"missing required section: {section!r}")
    if len(positions) == len(WORKFLOW_SECTIONS):
        ordered = [s for s, _ in sorted(positions, key=lambda p: p[1])]
        if ordered != list(WORKFLOW_SECTIONS):
            problems.append("required sections are out of order")
    for section in WORKFLOW_SECTIONS:
        if section not in headings:
            continue
        body = _section_body(text, section)
        if _is_placeholder_body(body):
            problems.append(f"{section!r} has no real body")
        elif require_files and _PLACEHOLDER_LINE.search(body):
            problems.append(f"{section!r} still contains a template placeholder line")

    if "reference.json" in text or re.search(r"\beval/", text):
        problems.append(
            "SKILL.md names eval material (eval/ or reference.json); it stays out of "
            "the skill body"
        )

    if require_files:
        for rel in sorted(set(_LOCAL_REF.findall(text))):
            if not (path.parent / rel).exists():
                problems.append(f"body names {rel} but the file does not exist")
        eval_dir = path.parent / "eval"
        for name in EVAL_FILES:
            f = eval_dir / name
            if not f.is_file() or f.stat().st_size == 0:
                problems.append(f"eval/{name} is missing or empty")
    return problems


SKILLS = list(iter_workflow_skills())


def test_at_least_one_workflow_skill_is_discovered():
    assert SKILLS, f"no skills/<name>/SKILL.md found under {SKILLS_DIR}"


@pytest.mark.parametrize("skill", SKILLS, ids=lambda p: p.parent.name)
def test_workflow_skill_conforms(skill: Path):
    problems = check_workflow_skill(skill)
    assert not problems, (
        f"skills/{skill.parent.name}/SKILL.md does not meet the workflow-skill "
        "contract:\n  "
        + "\n  ".join(problems)
        + "\n\nSee skills/workflow-conventions/SKILL.md."
    )


def test_contract_matches_the_conventions_document():
    conventions = (SKILLS_DIR / CONVENTIONS_SKILL_NAME / "SKILL.md").read_text()
    listed = re.findall(r"^\s*\d+\.\s+`(##\s[^`]+)`\s*$", conventions, re.MULTILINE)
    assert listed == list(WORKFLOW_SECTIONS), (
        "skills/workflow-conventions/SKILL.md and WORKFLOW_SECTIONS have diverged:\n"
        f"  prose lists: {listed}\n  code enforces: {list(WORKFLOW_SECTIONS)}"
    )


def test_the_shipped_template_satisfies_its_own_checker(tmp_path):
    conventions = (SKILLS_DIR / CONVENTIONS_SKILL_NAME / "SKILL.md").read_text()
    blocks = re.findall(r"```markdown\n(.*?)```", conventions, re.S)
    assert len(blocks) == 1, (
        f"expected one ```markdown template block, found {len(blocks)}"
    )
    fm = parse_frontmatter(blocks[0])
    assert fm and fm.get("name"), "template frontmatter must declare a name"
    folder = tmp_path / str(fm["name"])
    folder.mkdir()
    (folder / "SKILL.md").write_text(blocks[0])
    problems = check_workflow_skill(folder / "SKILL.md", require_files=False)
    assert problems == [], "\n  ".join(problems)


def test_an_unfilled_copy_of_the_template_is_rejected(tmp_path):
    """Copying the template and changing only the name must not pass."""
    conventions = (SKILLS_DIR / CONVENTIONS_SKILL_NAME / "SKILL.md").read_text()
    block = re.findall(r"```markdown\n(.*?)```", conventions, re.S)[0]
    folder = tmp_path / "copied-skill"
    (folder / "eval").mkdir(parents=True)
    for name in EVAL_FILES:
        (folder / "eval" / name).write_text("x\n")
    (folder / "SKILL.md").write_text(block.replace("<folder-name>", "copied-skill", 1))
    problems = check_workflow_skill(folder / "SKILL.md")
    assert any("template placeholder line" in p for p in problems), problems


def test_readme_catalog_matches_the_folders():
    """The README's catalog table and the folders must list the same skills."""
    readme = (SKILLS_DIR / "README.md").read_text()
    listed = set(re.findall(r"^\|\s*`hvantk:([a-z0-9-]+)`\s*\|", readme, re.MULTILINE))
    folders = {p.parent.name for p in SKILLS}
    assert listed == folders, (
        f"README lists {sorted(listed)}, folders are {sorted(folders)}"
    )
