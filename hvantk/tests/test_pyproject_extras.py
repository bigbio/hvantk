"""Guard the declared Poetry extras -- both what they contain and how they are documented.

#198: cptac imports pyranges, whose compiled dep `sorted_nearest` is not always
pulled by an extras install, breaking `import cptac`. The ptm extra must declare
`sorted-nearest` explicitly.

The docs half exists because the extras table is duplicated in THREE places -- the
`[project.optional-dependencies]` table, README.md, and
docs_site/getting-started/installation.md --
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


def _optional_dependencies() -> dict[str, list[str]]:
    """Extras as PEP 621 declares them: name -> list of full requirement specifiers."""
    return _load_pyproject()["project"]["optional-dependencies"]


def _requirement_name(spec: str) -> str:
    """`scikit-learn>=1.4,<2.0` -> `scikit-learn`. Bare names pass through unchanged."""
    return re.split(r"[<>=!~;@\[\s]", spec, maxsplit=1)[0].strip()


def _declared_extras() -> dict[str, set[str]]:
    """Extras as bare package names, for comparison against the prose tables."""
    return {k: {_requirement_name(s) for s in v}
            for k, v in _optional_dependencies().items()}


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
    assert "sorted-nearest" in _declared_extras()["ptm"]


def test_no_extra_duplicates_a_base_dependency():
    """A package in [project.dependencies] is always installed; gating it behind an extra
    would advertise a choice that does not exist."""
    base = {_requirement_name(s) for s in _load_pyproject()["project"]["dependencies"]}
    for extra, packages in _declared_extras().items():
        overlap = sorted(packages & base)
        assert not overlap, f"extra {extra!r} re-declares base dependencies {overlap}"


def test_extras_agree_on_every_shared_constraint():
    """PEP 621 puts the full specifier in each extra, so `scipy>=1.8` is written six times
    and `scikit-learn>=1.4,<2.0` three times. Nothing stops one from being re-pinned and the
    rest left behind -- an inconsistency that would resolve differently depending on which
    extra a user installed. This is the drift the old [tool.poetry.dependencies] split could
    not have, and it arrived with the migration, so it is guarded here."""
    seen: dict[str, dict[str, str]] = {}
    for extra, specs in _optional_dependencies().items():
        for spec in specs:
            seen.setdefault(_requirement_name(spec), {})[extra] = spec
    for pkg, by_extra in sorted(seen.items()):
        distinct = set(by_extra.values())
        assert len(distinct) == 1, f"{pkg} is spelled inconsistently across extras: {by_extra}"


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
