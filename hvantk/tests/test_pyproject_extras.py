"""Guard the declared Poetry extras -- both what they contain and how they are documented.

#198: cptac imports pyranges, whose compiled dep `sorted_nearest` is not always
pulled by an extras install, breaking `import cptac`. The ptm extra must declare
`sorted-nearest` explicitly.

The docs half exists because the extras table is duplicated in THREE places -- the
`[tool.poetry.extras]` block, README.md, and docs_site/getting-started/installation.md --
and only the first is executable. Three separate hand-fixes to the two prose copies were
needed in as many sessions, two of them caught only by adversarial review, and a fourth
drift (psroc/ancestry/ml missing scipy, ptm missing sorted-nearest -- eight wrong cells)
survived a release. A reader following a wrong table installs an environment that cannot run
the command the table promises, which is exactly the failure `pip install hvantk[constraint]`
produced. Deliberately non-Hail so it runs in the default suite.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
# (path, the column holding the dependency list). Both tables are markdown pipe tables whose
# first column is the extra name in backticks; they differ in column order, so each doc names
# its own header rather than assuming a position.
DOC_TABLES = [
    (ROOT / "README.md", "Pulls in"),
    (ROOT / "docs_site" / "getting-started" / "installation.md", "Pulls in"),
]


def _load_pyproject() -> dict:
    try:
        import tomllib as toml  # Python >= 3.11
    except ModuleNotFoundError:  # Python 3.10
        import tomli as toml
    return toml.loads((ROOT / "pyproject.toml").read_text())


def _declared_extras() -> dict[str, set[str]]:
    return {k: set(v) for k, v in _load_pyproject()["tool"]["poetry"]["extras"].items()}


def _parse_doc_table(path: Path, dep_header: str) -> dict[str, set[str]]:
    """Extras -> dependency set, read from the first markdown table that has `dep_header`.

    Only rows whose first cell is a single backticked token are treated as extras rows, so
    prose tables elsewhere in the file cannot be picked up by accident.
    """
    rows: dict[str, set[str]] = {}
    header_cols: list[str] | None = None
    done = False
    for line in path.read_text().splitlines():
        if not line.strip().startswith("|"):
            # The table ENDS at the first non-table line. Without this the parser would run
            # to end of file and read any later markdown table as more extras rows, using the
            # extras table's column index -- so a future unrelated table with a backticked
            # first column would silently overwrite a real row and mask the very drift this
            # test exists to catch.
            if header_cols is not None:
                done = True
            continue
        if done:
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if header_cols is None:
            if dep_header in cells:
                header_cols = cells
            continue
        if set("".join(cells)) <= set("-: "):        # the |---|---| separator
            continue
        name = re.fullmatch(r"`([a-z0-9-]+)`", cells[0])
        if not name:
            continue
        deps = cells[header_cols.index(dep_header)]
        rows[name.group(1)] = {d.strip().strip("`") for d in deps.split(",") if d.strip()}
    assert header_cols is not None, f"{path.name}: no table with a {dep_header!r} column"
    return rows


def test_parser_stops_at_the_end_of_the_extras_table(tmp_path):
    """A later table must not be read as more extras rows.

    Regression guard: the first version of this parser never reset on a non-table line, so it
    consumed the rest of the file. A second table whose first column is a backticked token --
    an ordinary thing to add to these docs -- would overwrite a real row using the extras
    table's column index, turning this guard into a source of false passes and false failures.
    """
    doc = tmp_path / "doc.md"
    doc.write_text(
        "| Extra | Pulls in |\n"
        "|---|---|\n"
        "| `viz` | matplotlib |\n"
        "\n"
        "Some prose between the tables.\n"
        "\n"
        "| Command | Pulls in |\n"
        "|---|---|\n"
        "| `viz` | something-else |\n"
    )
    assert _parse_doc_table(doc, "Pulls in") == {"viz": {"matplotlib"}}


def test_ptm_extra_includes_sorted_nearest():
    poetry = _load_pyproject()["tool"]["poetry"]
    assert "sorted-nearest" in poetry["extras"]["ptm"]
    assert "sorted-nearest" in poetry["dependencies"]


def test_every_extra_dependency_is_declared_optional():
    """An extra naming a package that is not an optional dependency installs nothing."""
    poetry = _load_pyproject()["tool"]["poetry"]
    deps = poetry["dependencies"]
    for extra, packages in poetry["extras"].items():
        for pkg in packages:
            assert pkg in deps, f"extra {extra!r} names undeclared dependency {pkg!r}"
            spec = deps[pkg]
            assert isinstance(spec, dict) and spec.get("optional") is True, \
                f"extra {extra!r} names {pkg!r}, which is not optional = true"


@pytest.mark.parametrize("path,dep_header", DOC_TABLES, ids=lambda p: getattr(p, "name", p))
def test_documented_extras_match_pyproject(path: Path, dep_header: str):
    declared = _declared_extras()
    documented = _parse_doc_table(path, dep_header)

    missing = sorted(set(declared) - set(documented))
    assert not missing, f"{path.name}: extras declared but undocumented: {missing}"
    unknown = sorted(set(documented) - set(declared))
    assert not unknown, f"{path.name}: extras documented but not declared: {unknown}"

    wrong = {
        name: {"documented": sorted(documented[name]), "declared": sorted(declared[name])}
        for name in sorted(declared)
        if documented[name] != declared[name]
    }
    assert not wrong, f"{path.name}: extras table disagrees with pyproject: {wrong}"
