"""Assert that the documentation's checkable claims still match the code.

Three classes of doc rot, each of which has actually shipped:

1. **Commands with options that do not exist.** `hvantk cohort validate
   --manifest ...` (real flag: `--cohort`), `hvantk enrichex burden --min-gq`
   (never existed -- genotype QC is upstream of EnrichEx), `hvantk download
   hgnc --output-dir` (that plugin takes `--output`). All three survived
   review and at least one release.

2. **A command that overwrites its own input.** A documented `annotate
   prepare --input X.ht --output X.ht` cannot run: Hail refuses a path that is
   both an input and an output of one query, and it fails only after full
   Spark startup.

3. **Structure trees naming paths that do not exist.** `docs_site/
   architecture.md` accumulated seven -- `core/constants.py`,
   `core/protocols.py`, `core/models/metadata.py`, `core/utils/catalog.py`,
   `core/utils/writers.py`, `algorithms/training_sets/`, `tools/build/` --
   while the parallel tree in README.md was corrected by hand. A contributor
   following the stale one is told to put new code in directories that are not
   there, in a layer that would violate the dependency rule.

Deliberately non-Hail so it runs in the default suite. It does import every
subcommand module (resolving a documented command imports the module that
implements it), which is the same cost `test_top_level_cli.py` already pays.
"""

from __future__ import annotations

import re
import shlex
import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner

from hvantk.hvantk import cli

REPO_ROOT = Path(__file__).resolve().parents[2]
PACKAGE_ROOT = REPO_ROOT / "hvantk"


def _tracked(*suffixes: str) -> list[Path]:
    """Repo files git actually tracks, with any of `suffixes`.

    Deliberately not `REPO_ROOT.rglob`: that also walks ignored and untracked
    files. The first version of this test did, and it failed in a working copy
    while passing in CI -- CI checks out a clean tree, so it never saw the
    gitignored agent-instruction files a developer has locally. A test that is
    green in CI and red on your machine trains people to ignore it. It also had
    no business reading files outside version control in the first place.
    """
    out = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "ls-files", "-z", *(f"*{s}" for s in suffixes)],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    return [REPO_ROOT / rel for rel in out.split("\0") if rel]


DOC_FILES = [
    p
    for p in _tracked(".md")
    if p.parts[len(REPO_ROOT.parts)] in {"docs_site", "README.md"}
]


# `hvantk <something>` inside a ```bash fence. Backslash continuations are
# joined first so a multi-line invocation is checked as one command.
_BASH_BLOCK = re.compile(r"```bash\n(.*?)```", re.S)
_INVOCATION = re.compile(r"^\s*(hvantk\s+.*?)(?=\n(?!\s)|\Z)", re.S | re.M)


def _documented_commands() -> list[tuple[Path, str]]:
    out: list[tuple[Path, str]] = []
    for path in DOC_FILES:
        for block in _BASH_BLOCK.findall(path.read_text()):
            for raw in _INVOCATION.findall(block.replace("\\\n", " ")):
                out.append((path, " ".join(raw.split())))
    return out


DOCUMENTED = _documented_commands()


def test_the_extractor_actually_finds_commands():
    """Guard the guard: a regex that silently matches nothing passes everything."""
    assert len(DOCUMENTED) > 100, (
        f"only {len(DOCUMENTED)} documented hvantk commands found -- the "
        "extractor is probably broken, which would make the tests below vacuous"
    )


@pytest.mark.parametrize(
    "path,command",
    [pytest.param(p, c, id=f"{p.name}::{c[:60]}") for p, c in DOCUMENTED],
)
def test_documented_command_options_exist(path: Path, command: str):
    """Every option in a documented command must exist on that command.

    Checked by appending `--help`, so click parses the option names without
    running anything: no file is read, no Hail context starts.
    """
    try:
        argv = shlex.split(command)[1:]
    except ValueError:  # an illustrative snippet with unbalanced quotes
        pytest.skip("not a parseable shell command")
    # Placeholders like `...` are prose, not arguments.
    argv = [a for a in argv if a != "..."]

    result = CliRunner().invoke(cli, argv + ["--help"], catch_exceptions=True)
    if result.exit_code != 0 and "No such option" in (result.output or ""):
        pytest.fail(
            f"{path.relative_to(REPO_ROOT)} documents a command with an option "
            f"that does not exist:\n    {command}\n"
            f"{result.output.strip().splitlines()[-1]}"
        )


@pytest.mark.parametrize(
    "path,command",
    [pytest.param(p, c, id=f"{p.name}::{c[:60]}") for p, c in DOCUMENTED],
)
def test_documented_command_does_not_overwrite_its_input(path: Path, command: str):
    """`--input X --output X` cannot run: Hail rejects a path that is both."""
    src = re.search(r"--input\s+(\S+)", command)
    dst = re.search(r"--output\s+(\S+)", command)
    if not (src and dst):
        return
    assert src.group(1) != dst.group(1), (
        f"{path.relative_to(REPO_ROOT)} documents a command that reads and "
        f"writes {src.group(1)}:\n    {command}\n"
        "Hail raises 'path is both an input and output source in a single "
        "query' -- after paying full Spark startup."
    )


def _tree_paths(markdown: str) -> set[str]:
    """Reconstruct the repo-relative path of every entry in an ASCII tree.

    Depth comes from the indentation before the `├──`/`└──` marker (4 columns
    per level), so `models/` nested under `core/` yields `core/models`, not a
    bare `models`. Matching on the basename alone is far too weak: a bogus
    `core/constants.py` would be "found" by any other constants.py in the tree.
    """
    paths: set[str] = set()
    # Any fenced block whose first line is a `hvantk/...` root. Trees are not
    # all rooted at the package root -- the Testing Strategy tree is rooted at
    # `hvantk/tests/` -- so each block carries its own prefix.
    for fence in re.findall(r"```[a-z]*\n(.*?)```", markdown, re.S):
        lines = fence.splitlines()
        if not lines or not re.fullmatch(r"hvantk/[A-Za-z0-9_/]*", lines[0].strip()):
            continue
        root = lines[0].strip()[len("hvantk/") :].rstrip("/")
        stack: dict[int, str] = {}
        for line in lines[1:]:
            match = re.search(r"([├└])── ([A-Za-z0-9_.*<>-]+/?)", line)
            if not match:
                continue
            depth = match.start(1) // 4
            name = match.group(2).rstrip("/")
            parent = stack.get(depth - 1, root)
            full = f"{parent}/{name}" if parent else name
            stack[depth] = full
            for deeper in [d for d in stack if d > depth]:
                del stack[deeper]
            paths.add(full)
    return paths


def _is_illustrative(path: str) -> bool:
    """Wildcards and <placeholders> stand for many real paths, not one."""
    return "*" in path or "<" in path or path.endswith("...")


@pytest.mark.parametrize("doc", ["README.md", "docs_site/architecture.md"])
def test_structure_tree_names_only_real_paths(doc: str):
    """Every entry drawn in a project-structure tree must exist at that path.

    architecture.md accumulated seven that did not: core/constants.py,
    core/protocols.py, core/models/metadata.py, core/utils/catalog.py,
    core/utils/writers.py, algorithms/training_sets/ and tools/build/. Two of
    those basenames DO exist elsewhere in the tree, which is why this checks
    the reconstructed path rather than the name.
    """
    drawn = _tree_paths((REPO_ROOT / doc).read_text())
    assert drawn, f"no hvantk/ structure tree found in {doc}"
    missing = sorted(
        path
        for path in drawn
        if not _is_illustrative(path) and not (PACKAGE_ROOT / path).exists()
    )
    assert not missing, (
        f"{doc} draws paths that do not exist:\n  "
        + "\n  ".join(f"hvantk/{m}" for m in missing)
        + "\nEither the tree is stale or the code moved."
    )


def test_the_tree_parser_reconstructs_nested_paths():
    """Guard the guard: if depth parsing broke, everything above goes vacuous."""
    drawn = _tree_paths((REPO_ROOT / "README.md").read_text())
    assert "core/models" in drawn, sorted(d for d in drawn if "models" in d)
    assert "core" in drawn and len(drawn) > 20


#: Names hvantk has presented as builder return types, current and historical.
#: Hardcoded rather than pattern-matched because Hail owns `Table` and
#: `MatrixTable`, which the HGC docs legitimately annotate with -- a shape-based
#: rule ("ends with Table/Matrix/Set") flags those and is unusable here.
#: `Artifact` is listed because it was deleted in #303 and restored as a union
#: later; that round trip is exactly the failure this test exists to catch.
HVANTK_ARTIFACT_NAMES = frozenset(
    {"Artifact", "AnnotationTable", "ExpressionMatrix", "VariantMatrix", "GeneSet"}
)


def test_documented_artifact_annotations_are_importable():
    """A `-> X` in the docs must be a real name in `hvantk.core.models`.

    #303 deleted the `Artifact` Protocol while four documents still wrote
    `-> Artifact` in copy-pasteable builder signatures, handing out-of-tree
    plugin authors an annotation that raised on import. The first version of
    this guard banned the string outright; `Artifact` is importable again -- as
    a union of the four concrete types -- so a ban would now forbid working
    code. What actually matters is that the name resolves.

    Concrete types remain the better style, since the manifest declares one and
    `run_builder_for_spec` checks against it. But style is not what breaks
    someone else's build.
    """
    import hvantk.core.models as models

    exported = set(models.__all__)
    stale = sorted(HVANTK_ARTIFACT_NAMES - exported)
    offenders = []
    for path in _tracked(".md", ".svg"):
        for match in re.finditer(r"(?:->|\u2192)\s*`?(\w+)`?", path.read_text()):
            name = match.group(1)
            if name in stale:
                offenders.append(f"{path.relative_to(REPO_ROOT)}: -> {name}")
    assert not offenders, (
        f"documentation annotates {stale}, which hvantk.core.models no longer "
        f"exports (it exports {sorted(exported)}): {offenders}"
    )
