"""Assert intra-core directional rules:

    core/io      may import   core/models
    core/models  may NOT import core/io
    core/utils   may NOT import core/models or core/io  (utils stays generic)

Same shape as test_dependency_directions.py: AST walk every import.
"""
from __future__ import annotations

import ast
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def _imports_in(subpath: str) -> list[tuple[Path, str]]:
    out: list[tuple[Path, str]] = []
    root = PACKAGE_ROOT / "core" / subpath
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


def _has(forbidden: str, imports: list[tuple[Path, str]]) -> list[tuple[Path, str]]:
    return [(p, m) for p, m in imports if m == forbidden or m.startswith(forbidden + ".")]


def test_core_models_does_not_import_core_io():
    bad = _has("hvantk.core.io", _imports_in("models"))
    assert not bad, (
        "core/models must not import from core/io. Offenders:\n"
        + "\n".join(f"  {p.relative_to(PACKAGE_ROOT)} -> {m}" for p, m in bad)
    )


def test_core_utils_does_not_import_core_models_or_io():
    bad = _has("hvantk.core.models", _imports_in("utils")) + \
          _has("hvantk.core.io", _imports_in("utils"))
    assert not bad, (
        "core/utils must stay generic (no core/models or core/io imports). "
        f"Offenders: {bad}"
    )


def test_provenance_unknown_is_quarantined():
    """Provenance.unknown is for the legacy shim and tests only.

    If a real plugin or algorithm calls it, drift detection silently fails.
    Allowlist: hvantk/core/io/_legacy.py and anything under hvantk/tests/.
    """
    pkg = PACKAGE_ROOT
    allowed = {
        pkg / "core" / "io" / "_legacy.py",
        pkg / "core" / "models" / "provenance.py",  # the definition itself
    }
    offenders: list[Path] = []
    for py in pkg.rglob("*.py"):
        if "__pycache__" in py.parts:
            continue
        if py in allowed:
            continue
        if "tests" in py.parts:
            continue  # tests can use it freely
        source = py.read_text()
        if "Provenance.unknown" in source:
            offenders.append(py)
    assert not offenders, (
        "Provenance.unknown is restricted to the legacy shim and tests. "
        "Offenders:\n" + "\n".join(f"  {p.relative_to(pkg)}" for p in offenders)
    )
