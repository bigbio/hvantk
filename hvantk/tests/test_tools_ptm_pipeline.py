"""Smoke test: the ptm workflow wrapper exists in tools/ptm/ after Phase C."""
from __future__ import annotations


def test_tools_ptm_exports_download_and_pipeline():
    from hvantk.tools.ptm.pipeline import (
        download_uniprot_ptm,
        ptm_build_pipeline,
    )
    assert callable(download_uniprot_ptm)
    assert callable(ptm_build_pipeline)


def test_algorithms_ptm_has_no_skill_imports():
    """AST-walk hvantk/algorithms/ptm/ asserting no hvantk.skills import.

    This is a tighter check than test_dependency_directions; it specifically
    catches any regression in the ptm subpackage.
    """
    import ast
    from pathlib import Path

    root = Path(__file__).resolve().parents[1] / "algorithms" / "ptm"
    bad = []
    for py in root.rglob("*.py"):
        if "__pycache__" in py.parts:
            continue
        tree = ast.parse(py.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name.startswith("hvantk.skills"):
                        bad.append((py.name, alias.name))
            elif isinstance(node, ast.ImportFrom) and node.module:
                if node.module.startswith("hvantk.skills"):
                    bad.append((py.name, node.module))
    assert not bad, f"algorithms/ptm/ has skill imports: {bad}"
