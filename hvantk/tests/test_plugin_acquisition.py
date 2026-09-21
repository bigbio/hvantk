"""Acquisition policy: the manifest must say whether a dataset CAN fetch its own data.

Before issue #118 the loader inferred this from whether ``lifecycle.download`` was
present, which conflated two structurally different states:

  * "a downloader has not been written yet" -- a TODO, and
  * "a downloader is impossible" -- size, credentials, licence, or publication-only
    supplementary data.

They were indistinguishable in a manifest, so ``--skip-download`` read as a category
error on every BYO invocation and no tooling could tell the two apart.

The check these tests care most about is the one that was NOT in the issue: that
``acquisition.mode: byo`` never gets read as "exempt from shipping a fixture". A BYO
dataset can be perfectly testable -- ``dbnsfp`` is BYO and ships a committed fixture
plus both snapshots -- and conflating acquisition with testability is what left five
datasets ungradable (#341).
"""

from __future__ import annotations

import json
from pathlib import Path

import jsonschema
import pytest
import yaml

from hvantk.core.plugin.api import Acquisition, PluginLoadError
from hvantk.core.plugin.loader import _acquisition_of

SKILLS_DIR = Path(__file__).resolve().parents[1] / "skills"
SCHEMA = json.loads(
    (
        Path(__file__).resolve().parents[1] / "core" / "plugin" / "manifest.schema.json"
    ).read_text()
)


def _manifest_with(acquisition: dict | None) -> dict:
    """A minimal VALID manifest, so only the acquisition block is under test.

    Validated whole rather than against ``properties.datasets.items``: the dataset
    sub-schema uses ``$ref: #/$defs/callable_ref``, which resolves against the document
    root, so validating the fragment alone raises PointerToNowhere and every case
    "fails" for a reason that has nothing to do with acquisition.
    """
    dataset = {
        "name": "d",
        "domain": "genomics",
        "backend": "hail",
        "builder": {"module": "m", "function": "f"},
        "drift_probe": {"module": "m", "function": "f"},
        "skill": "SKILL.md",
        "tests": {
            "command": "pytest",
            "fixture": "f",
            "schema_snapshot": "s",
            "row_snapshot": "r",
            "drift_fingerprint": "d",
        },
    }
    if acquisition is not None:
        dataset["acquisition"] = acquisition
    return {"api_version": 2, "name": "p", "version": "0.1.0", "datasets": [dataset]}


@pytest.mark.parametrize(
    "block, valid",
    [
        ({"mode": "download"}, True),
        ({"mode": "byo", "reason": "size"}, True),
        ({"mode": "byo", "reason": "license", "instructions": "SKILL.md#2"}, True),
        # A bare "byo" explains nothing, which is the whole point of the field.
        ({"mode": "byo"}, False),
        ({"mode": "maybe"}, False),
        ({"mode": "byo", "reason": "because-i-said-so"}, False),
        ({"reason": "size"}, False),  # mode is required
    ],
)
def test_acquisition_schema(block, valid):
    try:
        jsonschema.validate(_manifest_with(block), SCHEMA)
        ok = True
    except jsonschema.ValidationError:
        ok = False
    assert ok is valid


def test_absent_block_defaults_to_download():
    """Every manifest predating this field must keep loading unchanged."""
    acq = _acquisition_of({}, {}, "x:y")
    assert acq == Acquisition()
    assert acq.mode == "download"
    assert not acq.is_byo


def test_byo_and_lifecycle_download_is_a_load_error():
    """The one incoherent combination: it cannot fetch, and here is how it fetches.

    Whichever the runner honoured, the manifest would be lying about the other, so
    this is rejected rather than resolved by a precedence rule.
    """
    with pytest.raises(
        PluginLoadError, match="either fetches its own inputs or it does not"
    ):
        _acquisition_of(
            {"acquisition": {"mode": "byo", "reason": "size"}},
            {"download": {"module": "m", "function": "f"}},
            "onek-genomes:variants",
        )


def test_download_mode_alongside_a_downloader_is_fine():
    acq = _acquisition_of(
        {"acquisition": {"mode": "download"}},
        {"download": {"module": "m", "function": "f"}},
        "hgnc:lookup",
    )
    assert acq.mode == "download"


def _declared_acquisitions() -> dict[str, dict]:
    out = {}
    for manifest_path in sorted(SKILLS_DIR.glob("*/plugin.yaml")):
        manifest = yaml.safe_load(manifest_path.read_text())
        provider = manifest.get("name") or manifest_path.parent.name
        for ds in manifest.get("datasets", []):
            out[f"{provider}:{ds['name']}"] = ds
    return out


def test_byo_datasets_are_not_exempt_from_the_validation_contract():
    """BYO is about ACQUISITION, never about testability.

    dbnsfp is the worked example: nobody can download it (issue #321 -- every
    advertised archive 404s), yet it ships a committed fixture and both snapshots. If
    `acquisition.mode: byo` ever starts implying "no fixture needed", the ledger in
    test_plugin_contract_artifacts.py quietly becomes a dumping ground again.
    """
    from hvantk.tests.test_plugin_contract_artifacts import KNOWN_INCOMPLETE

    byo = {
        name: ds
        for name, ds in _declared_acquisitions().items()
        if (ds.get("acquisition") or {}).get("mode") == "byo"
    }
    # Every BYO dataset still declares a full tests: block ...
    for name, ds in byo.items():
        tests_block = ds.get("tests") or {}
        missing = [
            f
            for f in ("fixture", "schema_snapshot", "row_snapshot", "drift_fingerprint")
            if f not in tests_block
        ]
        assert not missing, f"{name} is BYO but stops declaring {missing}"

    # ... and being BYO does not put a dataset on the incomplete ledger.
    gradable_byo = sorted(set(byo) - set(KNOWN_INCOMPLETE))
    assert gradable_byo, (
        "no BYO dataset ships a full artifact set -- if that is ever true, the "
        "distinction this test defends has collapsed in practice"
    )


def test_every_dataset_without_a_downloader_has_been_classified():
    """A dataset with no downloader must say WHICH kind of no-downloader it is.

    Leaving it undeclared is what #118 is about: silence currently means both "not
    written yet" and "impossible", and the loader cannot tell them apart.
    """
    unclassified = []
    for name, ds in _declared_acquisitions().items():
        has_dl = bool((ds.get("lifecycle") or {}).get("download"))
        declared = (ds.get("acquisition") or {}).get("mode")
        if not has_dl and declared is None:
            unclassified.append(name)
    assert not unclassified, (
        "these datasets ship no lifecycle.download and declare no acquisition.mode, "
        "so 'not implemented yet' and 'impossible by design' are indistinguishable:\n  "
        + "\n  ".join(unclassified)
    )


def test_acquisition_vocabulary_matches_the_schema():
    """The Python enum and the JSON schema enum must be one vocabulary, not two.

    They had already drifted when this was written: the dataclass documented four
    reasons while the schema accepted five (``unstable-url`` was missing from the
    Python side), in the very commit that introduced both. A comment listing valid
    values is not a contract -- this is.

    Mirrors ``test_contract_matches_the_conventions_document``, which does the same
    job for the SKILL.md section list, and is bidirectional for the same reason: a
    value added to either side alone fails here.
    """
    from typing import get_args

    from hvantk.core.plugin.api import AcquisitionMode, AcquisitionReason

    props = SCHEMA["properties"]["datasets"]["items"]["properties"]["acquisition"][
        "properties"
    ]
    assert list(props["mode"]["enum"]) == list(get_args(AcquisitionMode))
    assert list(props["reason"]["enum"]) == list(get_args(AcquisitionReason))


@pytest.mark.parametrize(
    "acquisition, lifecycle, valid",
    [
        # The rule itself: a dataset either fetches its own inputs or it cannot.
        (
            {"mode": "byo", "reason": "size"},
            {"download": {"module": "m", "function": "f"}},
            False,
        ),
        # Each half alone is fine.
        ({"mode": "byo", "reason": "size"}, None, True),
        ({"mode": "download"}, {"download": {"module": "m", "function": "f"}}, True),
        # An omitted block defaults to "download", so a downloader must stay legal.
        (None, {"download": {"module": "m", "function": "f"}}, True),
        # `byo` beside a parse-only lifecycle is coherent: BYO is about acquisition.
        (
            {"mode": "byo", "reason": "license"},
            {"parse": {"module": "m", "function": "f"}},
            True,
        ),
    ],
)
def test_byo_and_lifecycle_download_are_mutually_exclusive_in_the_schema(
    acquisition, lifecycle, valid
):
    """`hvantk plugins validate` must reject what the loader rejects (#351).

    `_acquisition_of` raised `PluginLoadError` on this combination, but `plugins
    validate` runs jsonschema and not the loader, so a third-party author got `ok:
    plugin.yaml` followed by a plugin that would not load. Two validation surfaces
    disagreeing is worse than one, because the one people run before shipping was the
    wrong one.

    Encoded in the schema rather than duplicated into the CLI so the loader, the CLI
    and any future consumer inherit it together.
    """
    manifest = _manifest_with(acquisition)
    if lifecycle is not None:
        manifest["datasets"][0]["lifecycle"] = lifecycle

    try:
        jsonschema.validate(manifest, SCHEMA)
        accepted = True
    except jsonschema.ValidationError:
        accepted = False

    assert accepted is valid
