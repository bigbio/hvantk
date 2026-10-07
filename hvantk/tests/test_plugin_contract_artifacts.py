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
"""

from __future__ import annotations

import ast
from pathlib import Path
import shlex

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


def _has_marker(expression: ast.expr, marker: str) -> bool:
    if isinstance(expression, (ast.List, ast.Tuple, ast.Set)):
        return any(_has_marker(item, marker) for item in expression.elts)
    rendered = ast.unparse(expression)
    return rendered == f"pytest.mark.{marker}" or rendered.startswith(
        f"pytest.mark.{marker}("
    )


def _round_trip_test_files(command: str) -> list[Path]:
    """Resolve pytest path arguments from a manifest's test command."""
    tokens = shlex.split(command)
    try:
        pytest_index = tokens.index("pytest")
    except ValueError:
        return []
    repo_root = Path(__file__).resolve().parents[2]
    test_paths = [
        repo_root / token.split("::", 1)[0]
        for token in tokens[pytest_index + 1 :]
        if not token.startswith("-") and (repo_root / token.split("::", 1)[0]).exists()
    ]
    files = []
    for path in test_paths:
        files.extend(sorted(path.rglob("test_*.py")) if path.is_dir() else [path])
    return files


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


def test_dataset_round_trips_exist_and_are_not_marked_skipped():
    """Each dataset's test command must reach a round-trip test not marked to skip.

    The files a dataset's ``tests.command`` points pytest at must define a module-level
    ``test_*round_trip*`` function; none may carry a ``skip``/``skipif`` marker (as a
    decorator or through a module-level ``pytestmark``); and a hail-backed dataset needs
    one marked ``hail``. This is a static check of the source, not proof the test runs: a
    ``pytest.skip()`` call in the body, ``pytest.importorskip`` or a failing import
    still passes it. A run-time check is tracked in #432.
    """
    problems = []
    for manifest_path in sorted(SKILLS_DIR.glob("*/plugin.yaml")):
        manifest = yaml.safe_load(manifest_path.read_text())
        provider = manifest.get("name") or manifest_path.parent.name
        for dataset in manifest.get("datasets", []):
            tests = dataset.get("tests") or {}
            candidates = []
            for test_path in _round_trip_test_files(tests.get("command", "")):
                module = ast.parse(test_path.read_text())
                pytestmark_values = [
                    statement.value
                    for statement in module.body
                    if isinstance(statement, ast.Assign)
                    and any(
                        isinstance(target, ast.Name) and target.id == "pytestmark"
                        for target in statement.targets
                    )
                ]
                module_skipped = any(
                    _has_marker(value, "skip") or _has_marker(value, "skipif")
                    for value in pytestmark_values
                )
                module_hail = any(
                    _has_marker(value, "hail") for value in pytestmark_values
                )
                for node in module.body:
                    if (
                        isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                        and node.name.startswith("test_")
                        and "round_trip" in node.name
                    ):
                        skipped = module_skipped or any(
                            _has_marker(decorator, "skip")
                            or _has_marker(decorator, "skipif")
                            for decorator in node.decorator_list
                        )
                        has_hail = module_hail or any(
                            _has_marker(decorator, "hail")
                            for decorator in node.decorator_list
                        )
                        candidates.append((test_path, node.name, skipped, has_hail))
            if not candidates:
                problems.append(f"{provider}:{dataset['name']}: no round-trip test")
                continue
            skipped = [
                name for _, name, marked_skipped, _ in candidates if marked_skipped
            ]
            if skipped:
                problems.append(
                    f"{provider}:{dataset['name']}: skipped round-trip marker on "
                    + ", ".join(skipped)
                )
            if dataset.get("backend") == "hail" and not any(
                has_hail for _, _, _, has_hail in candidates
            ):
                problems.append(
                    f"{provider}:{dataset['name']}: no round-trip test is marked for Hail CI"
                )

    assert not problems, (
        "Every dataset's tests.command must reach a test_*round_trip* function with no "
        "skip/skipif marker, and a hail-marked one for a hail backend:\n"
        + "\n".join(f"  {problem}" for problem in problems)
    )
