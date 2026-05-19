"""Assert that hvantk's four-package layout (core / algorithms / skills / tools)
honors the one-way dependency rule documented in the design spec:

    skills/      ----+
                     +-->  algorithms/  -->  core/
    tools/       ----+

Implementation note: each layer-pair is checked independently and gated with
xfail until the corresponding migration phase fixes the violations. Phases
remove xfail markers as they land. Phase 7 deletes every xfail marker --
after that, this test enforces the contract for the future.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def _imports_in(layer: str) -> list[tuple[Path, str]]:
    """Return (file, dotted-import) pairs for every import statement in <layer>."""
    out: list[tuple[Path, str]] = []
    root = PACKAGE_ROOT / layer
    if not root.is_dir():
        return out
    for py in root.rglob("*.py"):
        if "__pycache__" in py.parts:
            continue
        try:
            tree = ast.parse(py.read_text())
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    out.append((py, alias.name))
            elif isinstance(node, ast.ImportFrom) and node.module:
                out.append((py, node.module))
    return out


def _forbidden_matches(layer: str, forbidden_prefixes: list[str]) -> list[tuple[Path, str]]:
    bad: list[tuple[Path, str]] = []
    for file, dotted in _imports_in(layer):
        for prefix in forbidden_prefixes:
            if dotted == prefix or dotted.startswith(prefix + "."):
                bad.append((file, dotted))
    return bad


def test_core_does_not_import_upward():
    bad = _forbidden_matches(
        "core", ["hvantk.algorithms", "hvantk.skills", "hvantk.tools"]
    )
    assert not bad, (
        "core/ must not import from algorithms/skills/tools. Offenders:\n"
        + "\n".join(f"  {p.relative_to(PACKAGE_ROOT)} -> {d}" for p, d in bad)
    )


def test_algorithms_does_not_import_skills_or_tools():
    bad = _forbidden_matches("algorithms", ["hvantk.skills", "hvantk.tools"])
    assert not bad, (
        "algorithms/ must not import from skills/ or tools/. Offenders:\n"
        + "\n".join(f"  {p.relative_to(PACKAGE_ROOT)} -> {d}" for p, d in bad)
    )


def test_skills_does_not_import_algorithms_or_tools():
    """Skills are siblings -- they meet only through core/. Already true today
    (Phase 1 of the original plugin migration enforced this for the 13
    migrated plugins). Should stay green; no xfail."""
    bad = _forbidden_matches("skills", ["hvantk.algorithms", "hvantk.tools"])
    sibling_bad = [
        (file, dotted)
        for file, dotted in _imports_in("skills")
        if dotted.startswith("hvantk.skills.")
        and dotted.split(".")[2] != file.relative_to(PACKAGE_ROOT / "skills").parts[0]
    ]
    assert not bad and not sibling_bad, (
        "skills/ must not import from algorithms/, tools/, or sibling skills/. "
        f"Offenders: {bad + sibling_bad}"
    )
