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

from pathlib import Path

import pytest

from hvantk.core.plugin.skill_spec import (
    REQUIRED_SECTIONS,
    check_skill_spec,
    iter_skill_specs,
)

SKILLS_DIR = Path(__file__).resolve().parents[1] / "skills"

SPECS = list(iter_skill_specs(SKILLS_DIR))


def _spec_id(path: Path) -> str:
    return path.relative_to(SKILLS_DIR).as_posix()


def test_specs_were_discovered():
    """Guard the guard: a broken glob would make every check below vacuously pass."""
    assert len(SPECS) >= 20, (
        f"expected the skills tree to hold at least 20 per-resource SKILL.md files, "
        f"found {len(SPECS)} -- has the layout moved?"
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
    for section in REQUIRED_SECTIONS:
        # _conventions lists them as numbered backtick-quoted items, not as headings.
        assert f"`{section}`" in conventions, (
            f"{section!r} is enforced by skill_spec.py but is not listed in "
            "_conventions/SKILL.md s 2; the contract and its checker have diverged."
        )
