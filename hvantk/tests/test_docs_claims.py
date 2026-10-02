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
   there, in a layer that would violate the dependency rule. That divergence is
   now structurally impossible: README.md's duplicate tree was replaced by a
   link, so architecture.md holds the only one and these tests guard it there.

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

    These run at module scope, so a failure here is a *collection* error for the
    whole file rather than a skip. A tree without `.git` is a real case -- a
    `git archive` export, an sdist, a container built from a copied source tree
    (this project's Apptainer image) -- and so is `git` not being on PATH.
    Neither says anything about the docs, so neither should read as a failure.
    """
    try:
        out = subprocess.run(
            [
                "git",
                "-C",
                str(REPO_ROOT),
                "ls-files",
                "-z",
                *(f"*{s}" for s in suffixes),
            ],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        pytest.skip(
            f"docs guards need a git checkout to enumerate tracked files: {exc}",
            allow_module_level=True,
        )
    return [REPO_ROOT / rel for rel in out.split("\0") if rel]


#: The prose whose `hvantk ...` invocations are checked. `docs_site/` and `README.md`
#: are read by humans; the per-provider `SKILL.md` files are read by coding agents that
#: then *write plugin code* from them, so a wrong command there is executed rather than
#: merely misleading. They were outside this filter until #350, and the promotion review
#: found ten invocations of `hvantk <name>-download` across five specs (plus thirty more
#: across seven downloader docstrings, which this filter still does not reach) -- a form
#: that has never existed, because the loader strips the `-download` suffix and binds the
#: command under the `download` group.
DOC_FILES = [
    p
    for p in _tracked(".md")
    if p.parts[len(REPO_ROOT.parts)] in {"docs_site", "README.md"}
    or (p.name == "SKILL.md" and "skills" in p.parts)
]


# `hvantk <something>` inside a ```bash fence. Backslash continuations are
# joined first so a multi-line invocation is checked as one command.
_BASH_BLOCK = re.compile(r"```bash\n(.*?)```", re.S)
_INVOCATION = re.compile(r"^\s*(hvantk\s+.*?)(?=\n(?!\s)|\Z)", re.S | re.M)


#: `hvantk ...` inside a single-backtick span. The fenced-block extractor above is not
#: enough for the per-provider specs: they name commands inline, in prose, so a fence-only
#: scan found 5 commands across 24 `SKILL.md` files and missed every one of the ten
#: `hvantk <name>-download` invocations #353 removed. Restricted to a single line so it
#: cannot swallow a fenced block, and to spans that start with `hvantk`.
_INLINE_INVOCATION = re.compile(r"`(hvantk\s+[^`\n]+)`")


#: Characters that mark a token as a stand-in rather than a literal subcommand:
#: `<source>`, `{list,show}`, `list/describe/errors`, and a bare ellipsis.
_PLACEHOLDER_CHARS = set("<>{}/,\u2026")


def _is_prose_shorthand(argv: list[str]) -> bool:
    """True when the *command path* is a placeholder rather than a real subcommand.

    Only the leading tokens are examined -- those before the first `-`-prefixed option
    -- because that is the part click has to resolve. A placeholder further right is an
    ordinary argument (`--raw-dir <dir>` is fine and common) and must stay checked, or
    widening this test to the specs would drop most of what it is meant to catch.

    Deliberately narrow: this skips documentation that never claimed to be literal
    (`hvantk download <source>`), not documentation that is simply wrong.
    """
    for token in argv:
        if token.startswith("-"):
            break
        if set(token) & _PLACEHOLDER_CHARS:
            return True
    return False


def _documented_commands_by_origin() -> dict[str, list[tuple[Path, str]]]:
    """Extracted commands, split by which extractor produced them.

    Split at the source rather than recounted afterwards so the floors below measure
    what actually reaches `DOCUMENTED`. A floor that re-runs the regexes independently
    guards the *patterns* but not the *wiring*: gutting this function's fenced loop left
    such a floor green, because it was still matching `_BASH_BLOCK` itself.
    """
    out: dict[str, list[tuple[Path, str]]] = {"fenced": [], "inline": []}
    for path in DOC_FILES:
        text = path.read_text()
        for block in _BASH_BLOCK.findall(text):
            for raw in _INVOCATION.findall(block.replace("\\\n", " ")):
                out["fenced"].append((path, " ".join(raw.split())))
        # Strip fenced blocks first so their contents are not counted twice.
        for raw in _INLINE_INVOCATION.findall(_BASH_BLOCK.sub("", text)):
            out["inline"].append((path, " ".join(raw.split())))
    return out


def _documented_commands() -> list[tuple[Path, str]]:
    by_origin = _documented_commands_by_origin()
    return by_origin["fenced"] + by_origin["inline"]


BY_ORIGIN = _documented_commands_by_origin()
DOCUMENTED = BY_ORIGIN["fenced"] + BY_ORIGIN["inline"]


def test_the_extractor_actually_finds_commands():
    """Guard the guard: a regex that silently matches nothing passes everything.

    Per-source floors, not one total. A single `> 100` bound was calibrated against a
    corpus a fifth this size, so by #356 it had stopped guarding anything: deleting the
    whole inline extractor still left 268, and dropping `SKILL.md` from `DOC_FILES`
    still left 373. Both are the very additions the bound was supposed to protect.

    The SKILL.md floor is the load-bearing one -- a fence-only scan found 5 commands
    across 24 specs, which is what let ten bogus `hvantk <name>-download` invocations sit
    in the provider docs until #353 removed them. (#348 reviewed the same promotion but
    did not touch them: all ten were still present at its tip.)

    Splitting by SOURCE FILE was not enough, because `docs_site/` is read by BOTH
    extractors. `from_docs > 100` claims to protect the fenced-block path, and does not:
    the inline path alone supplies 110 of those 373, so the fenced extractor can die
    completely and the bound still passes with margin. That is the same miscalibration
    this docstring says was fixed, one level down -- so the floors below are per
    EXTRACTOR, and each one is checked against what the other would leave behind.
    """
    from_specs = [c for p, c in DOCUMENTED if p.name == "SKILL.md"]
    from_docs = [c for p, c in DOCUMENTED if p.name != "SKILL.md"]

    # Read off what actually reached DOCUMENTED, not a recount from the regexes --
    # see _documented_commands_by_origin. Commands, not blocks: one ```bash fence
    # routinely holds several invocations.
    fenced = len(BY_ORIGIN["fenced"])
    inline = len(BY_ORIGIN["inline"])

    assert len(from_docs) > 100, (
        f"only {len(from_docs)} commands found in docs_site/README -- the extractors or "
        "the DOC_FILES filter are probably broken, which would make the tests below "
        "vacuous"
    )
    assert len(from_specs) > 100, (
        f"only {len(from_specs)} commands found in SKILL.md files -- the inline-span "
        "extractor or the DOC_FILES filter is probably broken. A fence-only scan finds "
        "about 5, so anything in that range means the inline path is gone."
    )
    assert fenced > 200, (
        f"only {fenced} commands from fenced ```bash blocks (was 268) -- the fenced "
        "extractor is broken. A floor on from_docs cannot catch this: the inline "
        "extractor alone supplies 110 of those 373."
    )
    assert inline > 180, (
        f"only {inline} commands from inline `hvantk ...` spans (was 242) -- the inline "
        "extractor is broken. SKILL.md uses backticks rather than fences, so this is "
        "the path that reaches the provider specs at all."
    )


def test_almost_every_extracted_command_is_actually_checked():
    """Guard the guard, one level in: extraction is not the same as assertion.

    `test_the_extractor_actually_finds_commands` bounds how many commands are FOUND.
    Nothing bounded how many are then `--help`-CHECKED, and `_is_prose_shorthand` sits
    between the two: widen it by accident -- a new placeholder token, a looser
    heuristic -- and every parametrised case below turns into a skip. The suite stays
    green, the count of collected tests does not move, and 488 assertions quietly stop
    running. A predicate that can disable the whole check is the one thing that most
    needs a bound of its own.

    The gap is deliberately narrow (22 of 510 today) because prose shorthand is rare:
    `<provider>:<dataset>` placeholders and `a/b/c`-style alternations. If a legitimate
    new shorthand pushes this over, raise the ratio and say which one -- do not delete
    the bound.
    """
    parseable = []
    for _, command in DOCUMENTED:
        try:
            argv = shlex.split(command)[1:]
        except ValueError:
            continue
        parseable.append([a for a in argv if a != "..."])

    checked = [a for a in parseable if not _is_prose_shorthand(a)]
    skipped = len(parseable) - len(checked)

    assert len(checked) > 450, (
        f"only {len(checked)} of {len(parseable)} documented commands survive "
        f"_is_prose_shorthand to reach a --help assertion ({skipped} skipped). The "
        "predicate has widened and the option checks below are now largely decorative."
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
    if _is_prose_shorthand(argv):
        pytest.skip("prose shorthand, not a literal invocation")

    result = CliRunner().invoke(cli, argv + ["--help"], catch_exceptions=True)
    # "No such command" matters as much as "No such option": a doc that says
    # `hvantk anotate ...` is exactly as broken for a reader as a bad flag, and
    # only the option case was checked -- so a typo'd or renamed *command* was
    # invisible to a test whose name claims to cover documented commands.
    output = result.output or ""
    broken = next(
        (m for m in ("No such option", "No such command") if m in output), None
    )
    if result.exit_code != 0 and broken is not None:
        noun = "an option" if broken == "No such option" else "a subcommand"
        pytest.fail(
            f"{path.relative_to(REPO_ROOT)} documents a command with {noun} "
            f"that does not exist:\n    {command}\n"
            f"{output.strip().splitlines()[-1]}"
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


# --- `reprocess <provider>:<dataset>` must name a dataset that actually exists -----------
#
# `test_documented_command_options_exist` appends `--help`, which click resolves
# WITHOUT ever looking at the `reprocess` command's positional argument -- `--help`
# short-circuits before click gets anywhere near resolving it against the plugin
# registry. So a stale or typo'd dataset name reads as a passing command there. That
# is exactly how `docs_site/getting-started/quickstart.md` and `docs_site/
# architecture.md` came to document `hvantk reprocess ucsc-cellbrowser:adultPancreas`,
# which is not a real dataset (the plugin ships `adult-ctx` / `dev-ctx`) and raises a
# bare `KeyError` if actually run.


def _documented_reprocess_invocations() -> list[tuple[Path, str, list[str], str]]:
    """(source file, full command, argv, dataset positional) for every documented,
    non-prose ``hvantk reprocess <dataset> ...`` invocation that names a dataset.

    A bare ``hvantk reprocess`` with no further token (prose naming the subcommand
    itself, e.g. "the `reprocess` command") is not an invocation of anything and is
    excluded, same as ``_is_prose_shorthand`` excludes a placeholder dataset like
    ``<provider>:<dataset>``.
    """
    out = []
    for path, command in DOCUMENTED:
        try:
            argv = [a for a in shlex.split(command)[1:] if a != "..."]
        except ValueError:
            continue
        if not argv or argv[0] != "reprocess" or _is_prose_shorthand(argv):
            continue
        dataset = next((a for a in argv[1:] if not a.startswith("-")), None)
        if dataset is None:
            continue
        out.append((path, command, argv, dataset))
    return out


def test_documented_reprocess_datasets_exist_in_the_registry():
    """Every documented ``hvantk reprocess <provider>:<dataset>`` must name a dataset
    the plugin registry actually knows about.

    Checked against ``list_manifests()`` -- the descriptive, pre-resolution view --
    rather than ``get_dataset()``: resolving callables can fail for reasons that have
    nothing to do with whether the documented NAME is real (a missing optional
    runtime, say), and that is not what this test is about.
    """
    from hvantk.core.plugin.loader import get_registry

    known = {m.name for m in get_registry().list_manifests()}
    offenders = [
        f"{path.relative_to(REPO_ROOT)}: {command}"
        for path, command, _argv, dataset in _documented_reprocess_invocations()
        if dataset not in known
    ]
    assert not offenders, (
        "documented `hvantk reprocess <provider>:<dataset>` names a dataset the "
        "plugin registry does not know (a --help check cannot catch this -- it "
        "short-circuits before the positional is ever resolved):\n  "
        + "\n  ".join(offenders)
    )


def test_documented_reprocess_without_a_downloader_passes_skip_download():
    """A dataset with no ``lifecycle.download`` and no ``acquisition.mode: byo`` has
    no downloader to run at all (tracked by issue #386). A documented invocation that
    omits ``--skip-download`` cannot actually run as written: `reprocess` would try
    to run a download step the plugin does not implement.

    Skips any dataset name the previous test already flags as unknown, so one stale
    name produces one failure here rather than two.
    """
    from hvantk.core.plugin.loader import get_registry

    manifests = {m.name: m for m in get_registry().list_manifests()}
    offenders = []
    for path, command, argv, dataset in _documented_reprocess_invocations():
        dm = manifests.get(dataset)
        if dm is None:
            continue
        needs_downloader_skip = not dm.has_download_fn and not dm.acquisition.is_byo
        if needs_downloader_skip and "--skip-download" not in argv:
            offenders.append(f"{path.relative_to(REPO_ROOT)}: {command}")
    assert not offenders, (
        "documented `hvantk reprocess <provider>:<dataset>` for a dataset with "
        "neither a downloader nor a BYO acquisition note must pass --skip-download, "
        "or the command cannot actually run as written (see issue #386):\n  "
        + "\n  ".join(offenders)
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


@pytest.mark.parametrize("doc", ["docs_site/architecture.md"])
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
    drawn = _tree_paths((REPO_ROOT / "docs_site" / "architecture.md").read_text())
    assert "core/models" in drawn, sorted(d for d in drawn if "models" in d)
    assert "core" in drawn and len(drawn) > 20


# --- backticked repo paths must exist ----------------------------------------------------
#
# `test_documented_command_options_exist` catches a wrong COMMAND, but a wrong PATH in prose
# rode through: `hvantk/resources/registry/...` (removed in the five-package refactor),
# `hvantk/qtlcascade/`, `hvantk/enrichex/`, `hvantk/algorithms/ptm/test.py` all survived
# in provider specs until #359. A path in backticks is a claim about the tree.

_REPO_PATH_SPAN = re.compile(r"`(hvantk/[^`\s]+)`")

#: Paths a document names precisely to say they do NOT exist. Keep this list short and
#: each entry justified; a stale entry here is the same rot this test exists to stop.
DELIBERATELY_ABSENT_PATHS = {
    # _conventions/SKILL.md: "There is no `hvantk/core/builders/table.py`" -- the sentence
    # explains why builders return artifacts rather than a shared base class.
    "hvantk/core/builders/table.py",
}


def _documented_repo_paths() -> list[tuple[str, str]]:
    """(source file, span) for every backticked hvantk/... path in the docs."""
    out = []
    files = sorted((REPO_ROOT / "hvantk" / "skills").rglob("SKILL.md")) + sorted(
        (REPO_ROOT / "docs_site").rglob("*.md")
    )
    for md in files:
        for span in _REPO_PATH_SPAN.findall(md.read_text()):
            out.append((str(md.relative_to(REPO_ROOT)), span))
    return out


def _path_part(span: str) -> str:
    """`hvantk/x/test_y.py::test_z` -> `hvantk/x/test_y.py`; `hvantk/x/c.py:SYMBOL` -> `hvantk/x/c.py`."""
    path = span.split("::", 1)[0]
    head, sep, tail = path.rpartition(":")
    if sep and tail and not tail.startswith("/") and "/" not in tail:
        path = head
    return path.rstrip("/")


def test_documented_repo_paths_exist():
    missing = sorted(
        {
            (src, span)
            for src, span in _documented_repo_paths()
            if not _is_illustrative(span)
            and _path_part(span) not in DELIBERATELY_ABSENT_PATHS
            and not (REPO_ROOT / _path_part(span)).exists()
        }
    )
    assert not missing, (
        "backticked hvantk/... paths that do not exist in the tree (fix the doc, or if the "
        "sentence is about the path NOT existing, add it to DELIBERATELY_ABSENT_PATHS with "
        f"a reason):\n" + "\n".join(f"  {src}: `{span}`" for src, span in missing)
    )


def test_deliberately_absent_paths_are_still_absent_and_still_cited():
    """The allowlist must not outlive the sentences it excuses."""
    cited = {_path_part(span) for _, span in _documented_repo_paths()}
    for path in DELIBERATELY_ABSENT_PATHS:
        assert not (REPO_ROOT / path).exists(), (
            f"{path} exists now; drop it from the allowlist"
        )
        assert path in cited, (
            f"{path} is no longer cited anywhere; drop it from the allowlist"
        )


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

    # Resolve each name for real. Comparing HVANTK_ARTIFACT_NAMES against
    # `__all__` alone compares two hand-maintained lists: both are literals in
    # this repo, so `stale` was permanently empty and the scan below was dead
    # code. `getattr` goes through the package's PEP 562 __getattr__ and so
    # actually imports the defining module -- deleting core/models/artifact.py
    # fails here, which is the whole point of a test named "are_importable".
    exported = set(models.__all__)
    unresolvable = []
    for name in sorted(HVANTK_ARTIFACT_NAMES & exported):
        try:
            getattr(models, name)
        except Exception as exc:  # noqa: BLE001 - report any import failure
            unresolvable.append(f"{name}: {type(exc).__name__}: {exc}")
    assert not unresolvable, (
        f"hvantk.core.models lists these in __all__ but cannot resolve them: {unresolvable}"
    )

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


def test_the_retired_download_command_form_is_gone_from_the_whole_tree():
    """`hvantk <name>-download` never existed and must not come back -- anywhere.

    `DOC_FILES` is built from `_tracked(".md")`, so no `.py` file is ever read. When
    #353 removed this form it found ten occurrences in specs and **thirty** more across
    seven `cli.py` module docstrings -- three times as many in exactly the files the
    markdown guard cannot see. Those docstrings are what `--help` prints, so they are
    read by more people than the specs are.

    The loader strips the `-download` suffix and binds the command under the `download`
    group, so the correct form is `hvantk download <name>`. This is a plain tree-wide
    grep rather than an extractor: the point is that no file type is exempt.
    """
    tracked = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split("\0")

    pattern = re.compile(r"hvantk\s+[a-z0-9_-]+-download\b")
    offenders: list[str] = []
    for rel in tracked:
        if not rel or not rel.endswith((".py", ".md", ".yaml", ".yml", ".rst")):
            continue
        path = REPO_ROOT / rel
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        # This test names the form in its own docstring; skip itself.
        if path.resolve() == Path(__file__).resolve():
            continue
        offenders.extend(f"{rel}: {m}" for m in pattern.findall(text))

    assert offenders == [], (
        "the retired `hvantk <name>-download` form is back in "
        f"{len(offenders)} place(s); use `hvantk download <name>`:\n  "
        + "\n  ".join(offenders)
    )
