"""Enforce the ``SKILL.md`` nine-section contract on every per-resource plugin spec.

``hvantk/skills/_conventions/SKILL.md`` s 2 has said "MUST" since the plugin system
landed, but nothing checked it, so the two drifted apart in silence: 9 of 23 provider
specs carried zero of the nine required headings and no YAML frontmatter (issue #334).
Documentation that is not executed is a wish.

This test closes that loop. It is intentionally a hard failure rather than a shrinking
ledger (contrast ``test_plugin_contract_artifacts.py``, which ratchets): a missing
snapshot needs upstream data that may not exist or may not be redistributable, whereas
a conforming spec needs only that someone writes it. There is no legitimate reason for
a new plugin to land non-conforming, so there is nothing to grandfather.

The rules live in ``hvantk.core.plugin.skill_spec`` -- one implementation, shared with
``hvantk plugins validate``. Stdlib-only, so this runs in the default ``pytest``
selection rather than only in the ``hail``-marked job.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml

from hvantk.core.plugin.skill_spec import (
    REQUIRED_SECTIONS,
    check_skill_spec,
    iter_skill_specs,
)

SKILLS_DIR = Path(__file__).resolve().parents[1] / "skills"

SPECS = list(iter_skill_specs(SKILLS_DIR))


def _spec_id(path: Path) -> str:
    return path.relative_to(SKILLS_DIR).as_posix()


def _manifest_declared_specs() -> set[Path]:
    """Every SKILL.md path declared by a `skill:` entry in some plugin.yaml.

    This is the independent inventory the glob is measured against. Read with
    ``yaml`` directly rather than through the registry so it needs no Hail.
    """
    declared: set[Path] = set()
    for manifest_path in sorted(SKILLS_DIR.glob("*/plugin.yaml")):
        manifest = yaml.safe_load(manifest_path.read_text())
        for dataset in manifest.get("datasets", []):
            rel = dataset.get("skill")
            if rel:
                declared.add((manifest_path.parent / rel).resolve())
    return declared


def test_discovery_matches_what_the_manifests_declare():
    """Guard the guard: a broken glob would make every check below vacuously pass.

    A bare lower bound (the first version asserted ``len(SPECS) >= 20``) does not do
    this: with 23 specs in the tree, up to three could move out of the two searched
    layouts and the count would still clear the bar, silently un-parametrizing them.
    Comparing against the paths the manifests themselves declare is exact.

    Only a subset relation is asserted in the declared->discovered direction, because
    several providers ship one `SKILL.md` shared by datasets that declare it
    individually, and a provider-level spec (``ensembl_gene/SKILL.md``) is discovered
    without any manifest naming it.
    """
    discovered = {p.resolve() for p in SPECS}
    declared = _manifest_declared_specs()
    assert (
        declared
    ), "no plugin.yaml declared a skill: -- has the manifest schema moved?"
    missing = sorted(str(p.relative_to(SKILLS_DIR)) for p in declared - discovered)
    assert not missing, (
        "these SKILL.md files are declared by a plugin.yaml but were NOT discovered by "
        f"iter_skill_specs, so nothing checks them: {', '.join(missing)}"
    )


@pytest.mark.parametrize("spec", SPECS, ids=_spec_id)
def test_skill_spec_conforms(spec: Path):
    """Every per-resource SKILL.md carries frontmatter and the nine sections, in order."""
    problems = check_skill_spec(spec)
    assert not problems, (
        f"{_spec_id(spec)} does not meet the SKILL.md contract:\n"
        + "\n".join(f"  - {p}" for p in problems)
        + "\n\nSee hvantk/skills/_conventions/SKILL.md s 2 for the contract and "
        "hvantk/skills/hgnc/SKILL.md for a conforming example."
    )


def test_contract_matches_the_conventions_document():
    """The checker's heading list must stay identical to the prose it enforces.

    Without this, someone could edit ``_conventions/SKILL.md`` to add a tenth
    mandatory section and the checker would keep passing on nine -- reintroducing
    exactly the documentation-vs-reality gap this module exists to close.
    """
    conventions = (SKILLS_DIR / "_conventions" / "SKILL.md").read_text()

    # _conventions s 2 lists them as a numbered list of backtick-quoted headings:
    #     1. `## 1. Status & scope`
    # Extract that list in document order and compare BOTH directions. The first
    # version of this test only asserted every enforced heading appeared somewhere in
    # the prose, which is the one direction that cannot catch the failure the
    # docstring describes: adding a tenth requirement to the prose left it green.
    listed = re.findall(r"^\s*\d+\.\s+`(##\s[^`]+)`\s*$", conventions, re.MULTILINE)

    assert listed == list(REQUIRED_SECTIONS), (
        "hvantk/skills/_conventions/SKILL.md s 2 and skill_spec.REQUIRED_SECTIONS have "
        "diverged -- the contract and its checker must stay identical.\n"
        f"  prose lists: {listed}\n"
        f"  code enforces: {list(REQUIRED_SECTIONS)}"
    )


# --- Negative cases. Each of these passed the first version of the checker, which is
# --- why they are pinned here: a conformance check that cannot fail is decoration.
@pytest.mark.parametrize(
    "label, mutate",
    [
        # A prefix match accepted a heading with anything appended to it.
        (
            "suffixed heading",
            lambda t: t.replace("## 1. Status & scope", "## 1. Status & scope (draft)"),
        ),
        # Searching the whole file accepted a heading that only appears quoted.
        (
            "heading only inside a fence",
            lambda t: t.replace(
                "## 2. Source identity", "```\n## 2. Source identity\n```"
            ),
        ),
        # startswith("---") accepted any line merely beginning with three dashes.
        ("malformed opening delimiter", lambda t: "---not-frontmatter\n" + t[4:]),
        # Splitting on the first ":" turned "name: [" into the string "[".
        (
            "unparseable frontmatter value",
            lambda t: t.replace("name: hvantk:resource-hgnc", "name: ["),
        ),
        # YAML keeps the last value, silently discarding the author's first one.
        (
            "duplicate frontmatter key",
            lambda t: t.replace(
                "status: provisional", "status: provisional\nname: sneaky"
            ),
        ),
    ],
)
def test_checker_rejects_non_conforming_variants(tmp_path, label, mutate):
    """A mutated copy of a known-good spec must be rejected, not silently accepted."""
    conforming = (SKILLS_DIR / "hgnc" / "SKILL.md").read_text()
    assert not check_skill_spec(
        SKILLS_DIR / "hgnc" / "SKILL.md"
    ), "baseline must conform"

    broken = tmp_path / "SKILL.md"
    broken.write_text(mutate(conforming))
    assert check_skill_spec(broken), f"checker accepted a spec with a {label}"
