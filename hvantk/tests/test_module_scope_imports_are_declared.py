"""Every third-party package imported at module scope must be declared in pyproject (#363).

`hvantk/skills/ucsc_cellbrowser/shared/ucsc.py` imported h5py at module scope while h5py
appeared in no dependency set -- masked only because anndata (a base dependency) happens
to require it. That is the pattern pyproject's own numpy comment warns about: "the day
pandas drops or re-pins numpy, the break surfaces somewhere unrelated." Nothing asserted
the rule, so this does. Stdlib-only, no Hail; runs in the default selection.
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "hvantk"

#: import name -> distribution name, where PyPI and the import disagree.
IMPORT_TO_DIST = {
    "sklearn": "scikit-learn",
    "yaml": "PyYAML",
    "sorted_nearest": "sorted-nearest",
}

#: Import names that are legitimately undeclared, each with the reason. Add here only
#: with a reason a reviewer can check; "it works on my machine" is what this test exists
#: to stop.
ALLOWED_UNDECLARED = {
    "hailtop": "part of the hail distribution, which is a base dependency",
    "pyspark": "installed by hail, which pins the compatible version",
}

_STDLIB = set(sys.stdlib_module_names)


def _normalize(dist: str) -> str:
    return re.sub(r"[-_.]+", "-", dist).lower()


def _declared_distributions() -> set[str]:
    try:
        import tomllib as toml  # Python >= 3.11
    except ModuleNotFoundError:  # Python 3.10
        import tomli as toml
    project = toml.loads((ROOT / "pyproject.toml").read_text())["project"]
    specs = list(project["dependencies"])
    for extra in project.get("optional-dependencies", {}).values():
        specs.extend(extra)
    return {_normalize(re.split(r"[<>=!~;@\[\s]", s, maxsplit=1)[0]) for s in specs}


def _source_files():
    for path in sorted(PACKAGE.rglob("*.py")):
        parts = path.relative_to(ROOT).parts
        if "tests" in parts or "testdata" in parts:
            continue
        yield path


def _is_type_checking_block(node: ast.AST) -> bool:
    if not isinstance(node, ast.If):
        return False
    test = node.test
    return (isinstance(test, ast.Name) and test.id == "TYPE_CHECKING") or (
        isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING"
    )


def _guards_import_error(node: ast.Try) -> bool:
    for handler in node.handlers:
        names = []
        if isinstance(handler.type, ast.Name):
            names = [handler.type.id]
        elif isinstance(handler.type, ast.Tuple):
            names = [e.id for e in handler.type.elts if isinstance(e, ast.Name)]
        if {"ImportError", "ModuleNotFoundError"} & set(names):
            return True
    return False


def module_scope_imports(path: Path) -> set[str]:
    """Top-level import names bound at module scope, excluding TYPE_CHECKING blocks and
    imports guarded by an except ImportError (those are optional by construction)."""
    tree = ast.parse(path.read_text(), filename=str(path))
    names: set[str] = set()

    def visit(body):
        for node in body:
            if isinstance(node, ast.Import):
                names.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                if node.level == 0 and node.module:
                    names.add(node.module.split(".")[0])
            elif _is_type_checking_block(node):
                visit(node.orelse)
            elif isinstance(node, ast.If):
                visit(node.body)
                visit(node.orelse)
            elif isinstance(node, ast.Try):
                if not _guards_import_error(node):
                    visit(node.body)
                for handler in node.handlers:
                    visit(handler.body)
                visit(node.orelse)
                visit(node.finalbody)

    visit(tree.body)
    return names


def _third_party(names: set[str]) -> set[str]:
    return {
        n for n in names if n not in _STDLIB and n != "hvantk" and not n.startswith("_")
    }


def test_the_scanner_sees_a_known_module_scope_import():
    """Guard the guard: if the scanner silently found nothing, the audit below would pass
    for the wrong reason."""
    ucsc = PACKAGE / "skills" / "ucsc_cellbrowser" / "shared" / "ucsc.py"
    assert "h5py" in module_scope_imports(ucsc)


def test_every_module_scope_third_party_import_is_declared():
    declared = _declared_distributions()
    offenders: dict[str, list[str]] = {}
    for path in _source_files():
        for name in _third_party(module_scope_imports(path)):
            if name in ALLOWED_UNDECLARED:
                continue
            if _normalize(IMPORT_TO_DIST.get(name, name)) in declared:
                continue
            offenders.setdefault(name, []).append(str(path.relative_to(ROOT)))
    assert not offenders, (
        "module-scope imports of packages declared in no dependency set (declare them in "
        "pyproject.toml [project.dependencies] or an extra, and mirror requirements.txt / "
        f"environment.yml):\n{offenders}"
    )
