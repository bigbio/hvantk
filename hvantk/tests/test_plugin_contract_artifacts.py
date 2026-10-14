"""Ratchet: every plugin.yaml artifact it declares must exist, or be listed as a known gap.

The plugin contract has each dataset declare a fixture, a schema snapshot, a row snapshot
and a drift fingerprint. The loader resolves those paths but never checks them
(``TestPaths`` is built from ``(plugin_dir / rel).resolve()``, and ``resolve()`` does not
require existence), so a manifest could promise a snapshot it did not ship and still load
clean -- which is how the tree reached 25 dataset declarations against 10 complete ones
with nothing failing.

This test makes the gap explicit and monotonically shrinking:

* a dataset that starts missing an artifact **fails immediately** -- new plugins cannot
  quietly land incomplete;
* a dataset listed here that becomes complete **also fails**, forcing the entry to be
  removed, so the ledger can never overstate the gap.

Deliberately reads the manifests with ``yaml`` instead of going through the registry, so
it needs no Hail and therefore runs in the default ``pytest`` selection. The equivalent
check behind the registry is ``TestPaths.missing_artifacts()``; keeping this one
Hail-free means the ratchet is enforced by every CI run, not only the ``hail``-marked job.

``test_dataset_round_trips_are_selected_and_unskipped`` below replaces an older static
AST scan of decorators on ``test_*round_trip*`` functions (issue #432). That scan read
source text, so it both missed real evasions -- ``pytest.skip()`` in the test body,
``pytest.importorskip``, a ``slow`` marker the default addopts deselects, a ``-k``/``-m``
filter buried in the manifest's own ``tests.command``, an aliased marker
(``skip_rt = pytest.mark.skip(...)``), a module-level ``pytestmark`` list, a skipped
``pytest.param``, ``collect_ignore`` in a conftest -- and raised false positives on valid
layouts it could not parse (a round trip written as a class method, ``skipif(False,
...)``, a round trip in a ``*_hail.py`` file, hail-marked only by the filename convention
``hvantk/tests/conftest.py`` applies). Worse, it never tied a round trip to *its*
dataset: deleting one test could stay green if a same-named unit test elsewhere in the
file happened to match the ``test_*round_trip*`` name pattern.

The replacement runs each dataset's own ``tests.command`` through ``pytest
--collect-only`` in a real subprocess (see ``hvantk/tests/_round_trip_introspection.py``)
and asks pytest's own collection and marker-evaluation machinery what it decided, instead
of re-deriving that from source:

* exactly one selected item must carry the dataset's ``round_trip`` marker -- a deleted
  marker, a misrouted one, and a ``-k``/``-m``/``collect_ignore``/``importorskip`` that
  keeps the item from ever being collected all surface the same way, as "0 matches" (or
  "N != 1" if a marker is duplicated across items, or stacked on one item);
* that item must not be something pytest's own skip/xfail evaluation would skip or
  xfail -- evaluating the actual condition, not just decorator presence, so
  ``skipif(False, ...)`` is correctly left alone;
* a ``backend: hail`` dataset's item must carry the ``hail`` marker, however it arrived
  (an explicit decorator or the ``*_hail.py`` filename convention).

What real collection still cannot see is a ``pytest.skip()`` call inside the test body,
because that only fires once the test actually runs, not at collection time. The
``--require-round-trips`` pytest hook in the repo-root ``conftest.py`` covers that half:
it fails the run if a ``round_trip``-marked test ends skipped or xfailed, and is wired
into the CI jobs that actually execute the round trips.

Nor does replaying ``tests.command`` faithfully prove CI can reach it: a command that
collects cleanly on its own can still run in no actual CI job, if its ``-m``/``-k``
differs from what CI's own invocations use, or its path sits outside ``hvantk/skills``
(``plugin-contract`` only walks that directory). ``test_tests_command_matches_how_ci_actually_invokes_it``
closes that with a fast, string-only check of the command's shape -- no subprocess, so
it runs in the default selection rather than behind ``round_trip_contract``.

Each subprocess costs a few seconds regardless of what it collects (pytest/Hail-adjacent
import time dominates), so across ~25 distinct ``tests.command`` strings the collection
check (not the fast shape check above) takes on the order of 30-40 seconds locally --
CI timing may differ, but either way it is too slow for the default ``pytest -q`` inner
loop. It is marked ``round_trip_contract``, a marker pytest.ini deselects by default
alongside slow/network/hail/integration; ``.github/workflows/python-package-conda.yml``'s
``plugin-contract`` job, the one CI job that actually runs the hail round trips, selects
it back in with ``-m``.
"""

from __future__ import annotations

import json
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

import pytest
import yaml

SKILLS_DIR = Path(__file__).resolve().parents[1] / "skills"

ARTIFACT_FIELDS = ("fixture", "schema_snapshot", "row_snapshot", "drift_fingerprint")

# Datasets that do NOT yet ship every declared artifact, with the fields still missing.
# Shrink this as snapshots land; never grow it for a new provider without a reason here.
#
# Empty since the plugin-completion work (#413, #414, #415, #417 and the alphagenome
# rewrite): every dataset now ships its fixture, both snapshots and its drift
# fingerprint. The dicts stay so the ratchet keeps guarding new plugins.
#
# Neither a licence that forbids redistributing real rows nor a source with no static
# file is a reason to sit here: such a constraint limits what a fixture may contain, never
# whether one exists. pqtl, cosmic-cgc and alphagenome used to sit here under that
# reasoning; pqtl and cosmic-cgc now ship synthetic, format-faithful fixtures, and
# alphagenome a small real subset under its NOTICE.md -- see "History of departures"
# below.
#
# Another cause used to dominate this list -- "simply never seeded, though the builder
# runs from a committed fixture". As of the #341 plugin follow-ups it is empty:
# expression-atlas and peptideatlas were the last two, and both now ship a fixture plus
# schema/row snapshots and a round-trip test, so both left this list. That matters
# beyond tidiness: a dataset with no executable contract cannot be graded, which is why
# the agent-authoring evaluation (#341) could say nothing at all about those five.
#
# History of departures, kept so the list is not re-grown by accident:
#   clingen, gencc, hgnc -- snapshots landed.
#   gevir -- drift fingerprint shipped with the gevir plugin-review work.
#   dbnsfp, gnomad-metrics -- real drift probes replaced their stub sentinels and their
#     baselines were captured from live probe runs. Both had been recorded under issue
#     #177 as having no probeable URL; re-checking showed the gnomAD constraint tables
#     sit in a public GCS bucket returning an MD5 ETag, and that the dbNSFP landing page
#     -- though its advertised S3 archives are all dead (issue #321) -- still exposes a
#     stable release list. (Both previously appeared here twice over: first for a fixture
#     dir never created, then for a missing drift fingerprint.)
#   uniprot-ptm, cptac:expression, cptac:phospho -- had no committed fixture at all,
#     their round-trip tests synthesized inputs into tmp_path. Small fixtures committed.
#   expression-atlas -- fixture derived by truncating the real 116k-transcript x
#     317-sample (320-column) inputs of E-MTAB-6798 (#341; the full copy once committed under hvantk/tests/testdata/raw/expression_atlas/ was later removed as unused).
#   peptideatlas -- fixture is a synthetic raw atlas_build_*.tsv.zip (real column
#     headers and real modification notation, confirmed against a live 202512/606
#     build; every row fabricated), so the round-trip test grades parse_raw_dir
#     together with the builder, not just the builder (#341; the raw-zip fixture
#     replaced an earlier parsed-intermediate one).
#   pqtl -- synthetic fixture (the Fang et al. preprint is cc_no); format checked
#     against the real allpairs files.
#   cosmic-cgc -- synthetic fixture (the COSMIC licence forbids real rows); header and
#     conventions checked against a licensed v103 export.
#   alphagenome -- builder rewritten to ingest SDK tidy_scores parquet; real-subset
#     fixture with NOTICE.md (AlphaGenome Output Terms, non-commercial).
#
# A future entry must be a dataset that cannot ship a *fixture* or *snapshot* at all;
# missing drift probes do not belong here (every dataset in the tree ships a live one).
KNOWN_INCOMPLETE: dict[str, tuple[str, ...]] = {}

# Why each remaining entry can never be completed. Every KNOWN_INCOMPLETE key must appear
# here (enforced below), so a future entry cannot join a list of permanent exemptions
# without someone writing down why it belongs there. If the honest answer is "not seeded
# yet", it does not belong in KNOWN_INCOMPLETE at all -- seed it, as #341 did for
# expression-atlas and peptideatlas.
PERMANENTLY_UNGRADABLE: dict[str, str] = {}


def _declared_artifact_gaps() -> dict[str, tuple[str, ...]]:
    """Map ``provider:dataset`` -> declared artifact fields whose file is absent."""
    gaps: dict[str, tuple[str, ...]] = {}
    for manifest_path in sorted(SKILLS_DIR.glob("*/plugin.yaml")):
        manifest = yaml.safe_load(manifest_path.read_text())
        provider = manifest.get("name") or manifest_path.parent.name
        for dataset in manifest.get("datasets", []):
            tests_block = dataset.get("tests") or {}
            missing = tuple(
                field
                for field in ARTIFACT_FIELDS
                if field in tests_block
                and not (manifest_path.parent / tests_block[field]).exists()
            )
            if missing:
                gaps[f"{provider}:{dataset['name']}"] = missing
    return gaps


def _manifest_datasets() -> list[dict[str, str]]:
    """Flatten every manifest dataset into its id, backend and ``tests.command``.

    Includes a dataset whose ``tests.command`` is empty (the schema requires the key
    but, before ``minLength: 1`` was added, not a non-empty value) instead of silently
    dropping it: an empty command used to make a dataset invisible to every check here.
    """
    datasets = []
    for manifest_path in sorted(SKILLS_DIR.glob("*/plugin.yaml")):
        manifest = yaml.safe_load(manifest_path.read_text())
        provider = manifest.get("name") or manifest_path.parent.name
        for dataset in manifest.get("datasets", []):
            datasets.append(
                {
                    "id": f"{provider}:{dataset['name']}",
                    "backend": dataset.get("backend", ""),
                    "command": (dataset.get("tests") or {}).get("command") or "",
                    "provider_dir": manifest_path.parent.name,
                }
            )
    return datasets


def _command_shape_problem(entry: dict[str, str]) -> str | None:
    """Flag a ``tests.command`` shaped differently from how CI actually invokes it.

    The collection check below replays a dataset's own ``tests.command`` faithfully --
    which means it cannot tell a command CI can reach from one that merely collects
    cleanly. ``plugin-contract`` (``.github/workflows/python-package-conda.yml``) only
    ever discovers a round trip by walking ``hvantk/skills`` with ``-m hail``; the
    default ``pytest`` selection walks the whole tree with no ``-m``/``-k`` override.
    Two commands pass this file's other checks while running in neither:

    * a non-hail dataset marked ``slow`` with ``command: ... -m 'slow or not slow'``
      collects its round trip cleanly, but the default addopts deselect ``slow`` and
      ``plugin-contract``'s ``-m hail`` was never going to select a non-hail test
      either, so it runs nowhere;
    * a command pointed at a path under ``hvantk/tests`` with ``-m hail`` collects
      cleanly on its own, but ``plugin-contract`` only walks ``hvantk/skills``, so it
      is never even considered there.

    Requiring the exact shape ``pytest hvantk/skills/<this dataset's own provider
    folder>/...`` plus, only when ``backend: hail``, exactly ``-m hail`` and nothing
    else, closes both: neither example matches it.
    """
    command = entry["command"]
    if not command:
        return "tests.command is empty"
    tokens = shlex.split(command)
    if not tokens or tokens[0] != "pytest":
        return f"must start with `pytest`, got {command!r}"
    # -m/-k each take the following token as their value, not a second path; everything
    # else that starts with "-" is a standalone flag (none of the current 26 commands
    # use one, but a bare boolean flag would otherwise misparse as a path).
    paths, flags, rest = [], [], tokens[1:]
    i = 0
    while i < len(rest):
        token = rest[i]
        if token in ("-m", "-k"):
            if i + 1 >= len(rest):
                return f"`{token}` with no value, got {command!r}"
            flags.extend([token, rest[i + 1]])
            i += 2
        elif token.startswith("-"):
            flags.append(token)
            i += 1
        else:
            paths.append(token)
            i += 1
    if len(paths) != 1:
        return f"must name exactly one path, got {command!r}"
    prefix = f"hvantk/skills/{entry['provider_dir']}/"
    if not paths[0].startswith(prefix):
        return f"path {paths[0]!r} is not under {prefix!r} (its own provider folder)"
    if entry["backend"] == "hail":
        if flags != ["-m", "hail"]:
            return (
                "a hail-backend dataset's command must end with exactly `-m hail`, "
                f"got {command!r}"
            )
    elif flags:
        return (
            f"a non-hail dataset's command must carry no extra flags, got {command!r}"
        )
    return None


def test_tests_command_matches_how_ci_actually_invokes_it():
    """Each dataset's ``tests.command`` must be shaped like a real CI invocation.

    Fast and string-only (no subprocess), so it runs in the default selection rather
    than behind ``round_trip_contract``. See ``_command_shape_problem`` for the two
    commands this specifically closes.
    """
    problems = []
    for entry in _manifest_datasets():
        problem = _command_shape_problem(entry)
        if problem is not None:
            problems.append(f"{entry['id']}: {problem}")
    assert not problems, (
        "Every dataset's tests.command must be `pytest hvantk/skills/<its own "
        "provider folder>/...`, with exactly `-m hail` appended when backend is hail "
        "and no other flags otherwise -- anything else can collect cleanly while "
        "running in no CI job:\n" + "\n".join(f"  {p}" for p in problems)
    )


# pytest's own exit codes (_pytest.config.ExitCode) that mean "collection happened
# without error", even though the second isn't 0: 5 (NO_TESTS_COLLECTED) is the normal
# outcome of e.g. a `-k`/`-m` filter matching nothing, which the "0 matches" check below
# already reports more specifically. Anything else -- 2 (INTERRUPTED), most often a
# collection error in a sibling file -- means the report this run produced cannot be
# trusted: a broken neighbour does not stop pytest from still finding and reporting the
# round_trip item in the file that imports fine, so relying on "N matches" alone missed
# it (verified: a broken sibling under hvantk/skills/pqtl/tests left pqtl:metrics's row
# intact in the JSON report while pytest exited 2).
_OK_COLLECT_EXIT_CODES = (0, 5)


def _collect_round_trip_state(command: str, repo_root: Path) -> dict:
    """Run a manifest's own ``tests.command`` with ``--collect-only`` in a subprocess.

    Returns ``{"rows": [...], "returncode": int, "diagnostics": str}``. ``rows`` has
    one dict per item carrying a ``round_trip`` marker that pytest would actually run
    (i.e. after that command's own ``-k``/``-m`` deselection), as reported by
    ``_round_trip_introspection`` (loaded with ``-p``, so nothing about this repo's
    normal conftest chain has to change to make this observable). ``diagnostics`` is
    the last lines of stdout and stderr combined: pytest's own "ERROR collecting ..."
    report for a collection error -- the main reason to look at this -- goes to
    stdout, not stderr, so stderr alone is usually empty (verified).
    """
    args = [token for token in shlex.split(command) if token != "pytest"]
    with tempfile.TemporaryDirectory() as tmp_dir:
        report_path = Path(tmp_dir) / "round_trip_report.json"
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                *args,
                "--collect-only",
                "-q",
                "-p",
                "hvantk.tests._round_trip_introspection",
                f"--round-trip-report={report_path}",
            ],
            cwd=repo_root,
            capture_output=True,
            text=True,
            timeout=120,
        )
        if not report_path.exists():
            raise AssertionError(
                f"`{command} --collect-only` produced no report (exit "
                f"{result.returncode}):\n{result.stdout}\n{result.stderr}"
            )
        tail_lines = (result.stdout + result.stderr).splitlines()[-20:]
        return {
            "rows": json.loads(report_path.read_text()),
            "returncode": result.returncode,
            "diagnostics": "\n".join(tail_lines),
        }


def test_no_new_dataset_declares_a_missing_artifact():
    """A dataset may not declare an artifact it does not ship unless it is a known gap."""
    unexpected = {
        name: fields
        for name, fields in _declared_artifact_gaps().items()
        if name not in KNOWN_INCOMPLETE
    }
    assert not unexpected, (
        "These datasets declare validation artifacts that do not exist on disk:\n"
        + "\n".join(
            f"  {n}: missing {', '.join(f)}" for n, f in sorted(unexpected.items())
        )
        + "\n\nCommit the artifacts (see hvantk/skills/_conventions/SKILL.md), or add an "
        "entry to KNOWN_INCOMPLETE with the reason."
    )


def test_known_incomplete_ledger_is_not_stale():
    """Once a dataset ships everything, it must be dropped from the ledger."""
    gaps = _declared_artifact_gaps()
    now_complete = sorted(name for name in KNOWN_INCOMPLETE if name not in gaps)
    assert not now_complete, (
        "These datasets now ship every declared artifact and must be removed from "
        f"KNOWN_INCOMPLETE: {', '.join(now_complete)}"
    )


@pytest.mark.parametrize("name", sorted(KNOWN_INCOMPLETE))
def test_known_gap_matches_reality(name: str):
    """The recorded fields must match what is actually absent, so the ledger stays honest."""
    actual = _declared_artifact_gaps().get(name, ())
    assert set(actual) == set(KNOWN_INCOMPLETE[name]), (
        f"{name}: ledger says {sorted(KNOWN_INCOMPLETE[name])}, "
        f"disk says {sorted(actual)}"
    )


def test_every_known_gap_is_classified_permanent():
    """A dataset may only sit on the ledger with a written reason it can never leave it.

    The ledger's original third cause was "simply never seeded", which is not an
    exemption but a TODO -- and a TODO parked in an exemption list is how five datasets
    came to have no executable contract, which is what made them ungradable by the
    agent-authoring evaluation (#341). Requiring a ``PERMANENTLY_UNGRADABLE`` entry
    forces that distinction to be made explicitly at the moment a dataset is added.
    """
    # A blank or whitespace reason is not a reason. Comparing keys alone would let
    # PERMANENTLY_UNGRADABLE[name] = "" satisfy the invariant while explaining nothing.
    unclassified = sorted(
        name
        for name in KNOWN_INCOMPLETE
        if not str(PERMANENTLY_UNGRADABLE.get(name, "")).strip()
    )
    assert not unclassified, (
        "These datasets are on KNOWN_INCOMPLETE without a recorded reason they can "
        f"never ship the artifacts: {', '.join(unclassified)}.\n"
        "If the real reason is 'not seeded yet', seed it instead of listing it."
    )
    stale = sorted(set(PERMANENTLY_UNGRADABLE) - set(KNOWN_INCOMPLETE))
    assert not stale, (
        "These datasets carry a permanent-exemption reason but are no longer on "
        f"KNOWN_INCOMPLETE: {', '.join(stale)}. Remove the reason too."
    )


@pytest.mark.round_trip_contract
def test_dataset_round_trips_are_selected_and_unskipped():
    """Each dataset's ``tests.command`` must select exactly one clean round-trip item.

    Runs real collection per distinct ``tests.command`` (grouping datasets that share
    one, e.g. onek-genomes's two datasets and ucsc-cellbrowser's three), so it sees
    precisely what pytest would run: markers resolved through decorators, ``pytestmark``
    lists and aliases alike, with that command's own ``-k``/``-m`` deselection already
    applied. See the module docstring for what this replaces and why.
    """
    repo_root = Path(__file__).resolve().parents[2]
    by_command: dict[str, list[dict[str, str]]] = {}
    problems = []
    for entry in _manifest_datasets():
        if not entry["command"]:
            # An empty command has nothing to run as a subprocess (shelling out with no
            # path args would instead collect pytest.ini's whole testpaths), so this is
            # reported directly rather than joining by_command below.
            problems.append(f"{entry['id']}: tests.command is empty")
            continue
        by_command.setdefault(entry["command"], []).append(entry)

    for command, entries in sorted(by_command.items()):
        try:
            collected = _collect_round_trip_state(command, repo_root)
        except subprocess.TimeoutExpired:
            problems.append(f"`{command}`: collection timed out")
            continue
        except AssertionError as exc:
            problems.append(str(exc))
            continue

        if collected["returncode"] not in _OK_COLLECT_EXIT_CODES:
            problems.append(
                f"`{command} --collect-only` exited {collected['returncode']}, which "
                "usually means a collection error (e.g. a broken sibling test file) "
                "that a clean round_trip match would not otherwise reveal:\n"
                f"{collected['diagnostics']}"
            )

        matches: dict[str, list[dict]] = {}
        for row in collected["rows"]:
            if len(row["round_trip"]) > 1:
                problems.append(
                    f"{row['nodeid']}: carries {len(row['round_trip'])} round_trip "
                    f"markers {row['round_trip']} (want at most 1 per item)"
                )
            for dataset_id in row["round_trip"]:
                matches.setdefault(dataset_id, []).append(row)

        for entry in entries:
            selected = matches.get(entry["id"], [])
            if len(selected) != 1:
                problems.append(
                    f"{entry['id']}: `{command}` selects {len(selected)} item(s) "
                    "carrying its round_trip marker (want exactly 1)"
                )
                continue
            (row,) = selected
            if row["skip_or_xfail"]:
                problems.append(f"{entry['id']}: {row['nodeid']} is skipped or xfailed")
            if entry["backend"] == "hail" and not row["hail"]:
                problems.append(f"{entry['id']}: {row['nodeid']} is not marked hail")

    assert not problems, (
        "Each dataset's tests.command must select exactly one round-trip item, clean of "
        "skip/xfail, and hail-marked when its backend is hail:\n"
        + "\n".join(f"  {problem}" for problem in problems)
    )
