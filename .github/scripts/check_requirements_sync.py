"""Check that requirements.txt and environment.yml stay in sync with pyproject.toml.

pyproject.toml's `[project.dependencies]` is the source of truth for hvantk's base
runtime dependencies. requirements.txt and environment.yml exist only because the
pip-path CI jobs (python-app.yml, python-package.yml) and the conda-path job
(python-package-conda.yml) install from them rather than from pyproject/poetry.lock,
and both files say in their own header/comments that they are kept in sync BY HAND --
nothing checked that until now (#408).

Three deterministic rules, run against the current tree:
  (a) every base dependency in pyproject's [project.dependencies] appears, by
      normalized (PEP 503) name, in requirements.txt and in environment.yml (its
      conda `dependencies:` list or the nested `pip:` list) -- a hard failure.
  (b) every package in requirements.txt is declared somewhere in pyproject.toml
      (base, any [project.optional-dependencies] extra, or the
      [tool.poetry.group.dev.dependencies] group) -- a hard failure.
  (c) base dependencies whose lower version bound differs between files are printed
      as warnings only. Some divergence here is intentional (environment.yml pins
      `hail` exactly for Spark compatibility, and a couple of entries there carry no
      floor at all), so this never fails the build -- it is a tripwire for someone to
      read, not a gate.

Stdlib only, including the TOML read (`tomllib`, Python >= 3.11 -- run this with such
an interpreter). environment.yml has no stdlib YAML reader available, so it is read
with a small line-based scan narrow enough for this file's actual shape (a flat conda
`dependencies:` list with exactly one nested `pip:` block) rather than a general
parser; see `_environment_yml_entries`.

Exit 1 on (a)/(b) violations, 0 otherwise -- (c) never changes the exit code.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    print(
        "FAIL check_requirements_sync.py needs Python >= 3.11 for tomllib, got "
        f"{sys.version.split()[0]} -- run it with python3.11+",
        file=sys.stderr,
    )
    sys.exit(1)

REPO = Path(__file__).resolve().parents[2]

#: Conda package names that differ from their PyPI distribution name, mapped to the
#: PyPI name. Empty today: every base dependency in environment.yml is installed
#: through its nested `pip:` list, under its PyPI name already. Add an entry here
#: (conda name -> PyPI name) the day a base dependency moves to the top-level conda
#: `dependencies:` list under a name that disagrees with PyPI's.
CONDA_TO_PYPI: dict[str, str] = {}

_NAME_RE = re.compile(r"^([A-Za-z0-9][A-Za-z0-9._-]*)")
_LOWER_BOUND_RE = re.compile(r"(?:>=|==)\s*([0-9]+(?:\.[0-9]+)*)")


def _normalize(name: str) -> str:
    """PEP 503: case-fold and collapse runs of -_. to a single '-'."""
    return re.sub(r"[-_.]+", "-", name).lower()


def _split_name_and_spec(entry: str) -> tuple[str, str] | None:
    """Split "pandas>=2.0,<3" into ("pandas", ">=2.0,<3"); None if unparseable."""
    entry = entry.strip()
    if not entry:
        return None
    m = _NAME_RE.match(entry)
    if not m:
        return None
    name = m.group(1)
    return name, entry[len(name) :]


def _lower_bound(spec: str) -> tuple[int, ...] | None:
    """The `>=` or `==` floor of a specifier, as an int tuple; None if unconstrained.

    `==` counts as a floor too -- environment.yml pins `hail` exactly (`==0.2.137`)
    where the other files say `>=0.2.137`, and that pair should read as agreeing, not
    as drift, since the pinned value matches the floor everywhere else.
    """
    m = _LOWER_BOUND_RE.search(spec)
    if not m:
        return None
    return tuple(int(p) for p in m.group(1).split("."))


def _bounds_equal(a: tuple[int, ...] | None, b: tuple[int, ...] | None) -> bool:
    if a is None or b is None:
        return a == b
    n = max(len(a), len(b))
    return a + (0,) * (n - len(a)) == b + (0,) * (n - len(b))


def _render_bound(b: tuple[int, ...] | None) -> str:
    return "unconstrained" if b is None else ".".join(str(p) for p in b)


def _load_pyproject() -> dict:
    return tomllib.loads((REPO / "pyproject.toml").read_text())


def _pyproject_base_deps(data: dict) -> dict[str, str]:
    """normalized name -> raw specifier, from [project.dependencies]."""
    out: dict[str, str] = {}
    for entry in data["project"]["dependencies"]:
        parsed = _split_name_and_spec(entry)
        if parsed:
            out[_normalize(parsed[0])] = parsed[1]
    return out


def _pyproject_extra_and_dev_names(data: dict) -> set[str]:
    """normalized names declared in any [project.optional-dependencies] extra or in
    [tool.poetry.group.dev.dependencies] (a Poetry-only table: keys are bare names)."""
    names: set[str] = set()
    for extra_entries in data["project"].get("optional-dependencies", {}).values():
        for entry in extra_entries:
            parsed = _split_name_and_spec(entry)
            if parsed:
                names.add(_normalize(parsed[0]))
    dev_group = (
        data.get("tool", {})
        .get("poetry", {})
        .get("group", {})
        .get("dev", {})
        .get("dependencies", {})
    )
    names.update(_normalize(name) for name in dev_group)
    return names


def _requirements_txt_entries() -> dict[str, str]:
    """normalized name -> raw specifier, from requirements.txt."""
    out: dict[str, str] = {}
    for line in (REPO / "requirements.txt").read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        parsed = _split_name_and_spec(line)
        if parsed:
            out[_normalize(parsed[0])] = parsed[1]
    return out


def _environment_yml_entries() -> dict[str, str]:
    """normalized name -> raw specifier, read from environment.yml's top-level
    `dependencies:` list and its nested `pip:` list.

    Not a general YAML reader. environment.yml's shape here is a flat list of scalar
    conda specs with exactly one nested mapping item (`- pip:` followed by a
    more-indented list), and this reads exactly that shape by indentation. A file
    restructured into something else (a second nested list, `pip:` spelled
    differently) needs an actual YAML parser instead of a wider version of this.
    """
    out: dict[str, str] = {}
    in_deps = False
    top_indent: int | None = None

    for raw in (REPO / "environment.yml").read_text().splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip(" "))
        stripped = line.strip()

        if not in_deps:
            in_deps = stripped == "dependencies:"
            continue

        if top_indent is None:
            top_indent = indent
        elif indent < top_indent:
            break  # dedented past the list -- a new top-level key started

        if indent == top_indent and stripped in ("- pip:", "pip:"):
            continue  # the "pip:" marker itself, not a package

        item = stripped[2:] if stripped.startswith("- ") else stripped
        parsed = _split_name_and_spec(item)
        if not parsed:
            continue
        pypi_name = CONDA_TO_PYPI.get(parsed[0], parsed[0])
        out[_normalize(pypi_name)] = parsed[1]

    return out


def _rule_a_and_b_failures(
    base: dict[str, str],
    requirements: dict[str, str],
    environment: dict[str, str],
    declared_anywhere: set[str],
) -> list[str]:
    """(a) every base dep is in both files; (b) every requirements.txt package is
    declared somewhere in pyproject.toml. Hard failures."""
    failures: list[str] = []
    for name in sorted(base):
        if name not in requirements:
            failures.append(
                f"{name!r}: a base dependency in pyproject.toml, missing from requirements.txt"
            )
        if name not in environment:
            failures.append(
                f"{name!r}: a base dependency in pyproject.toml, missing from environment.yml"
            )
    for name in sorted(requirements):
        if name not in declared_anywhere:
            failures.append(
                f"{name!r}: listed in requirements.txt but declared in no pyproject.toml "
                "dependency set (base, an extra, or the poetry dev group)"
            )
    return failures


def _rule_c_warnings(
    base: dict[str, str], requirements: dict[str, str], environment: dict[str, str]
) -> list[str]:
    """(c) base dependencies whose lower bound differs between files. Non-fatal."""
    warnings: list[str] = []
    for name in sorted(base):
        bounds = {
            "pyproject.toml": _lower_bound(base[name]),
            "requirements.txt": _lower_bound(requirements[name]),
            "environment.yml": _lower_bound(environment[name]),
        }
        reference = bounds["pyproject.toml"]
        if any(not _bounds_equal(reference, b) for b in bounds.values()):
            rendered = ", ".join(f"{f}={_render_bound(b)}" for f, b in bounds.items())
            warnings.append(f"{name!r}: lower bound differs -- {rendered}")
    return warnings


def main() -> int:
    data = _load_pyproject()
    base = _pyproject_base_deps(data)
    declared_anywhere = set(base) | _pyproject_extra_and_dev_names(data)
    requirements = _requirements_txt_entries()
    environment = _environment_yml_entries()

    failures = _rule_a_and_b_failures(
        base, requirements, environment, declared_anywhere
    )
    if failures:
        print("FAIL requirements-sync check:", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        return 1

    warnings = _rule_c_warnings(base, requirements, environment)
    if warnings:
        print("WARN requirements-sync lower-bound drift (non-fatal):", file=sys.stderr)
        for w in warnings:
            print(f"  - {w}", file=sys.stderr)

    print(
        "OK requirements.txt and environment.yml agree with pyproject.toml"
        + (" (see warnings above)" if warnings else "")
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
