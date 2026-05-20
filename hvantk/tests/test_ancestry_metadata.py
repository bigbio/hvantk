"""Phase D smoke test: ancestry pipeline is registered with @algorithm metadata."""
from __future__ import annotations


def test_ancestry_has_algorithm_metadata():
    from hvantk.algorithms.ancestry.pipeline import run_ancestry_inference
    from hvantk.core.models.backends import get_algorithm_meta, Backend

    meta = get_algorithm_meta(run_ancestry_inference)
    assert meta.name == "ancestry_inference"
    assert Backend.HAIL in meta.backends


def test_ancestry_has_no_skill_imports():
    """AST-walk the ancestry subpackage."""
    import ast
    from pathlib import Path

    root = Path(__file__).resolve().parents[1] / "algorithms" / "ancestry"
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
    assert not bad, f"algorithms/ancestry/ has skill imports: {bad}"
