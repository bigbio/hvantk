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
        # An explicit mode 'download' with no sibling `lifecycle.download` (as added by
        # `_manifest_with` below) is the #360 hole: it used to validate here and only
        # fail once `hvantk reprocess` actually ran it.
        ({"mode": "download"}, False),
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


#: Datasets whose downloader belongs here but is not written yet (#386). That state is
#: spelled by OMITTING the acquisition block; this set pins which datasets are in it.
DOWNLOADER_NOT_WRITTEN_YET = {
    "cptac:expression",
    "gevir:metrics",
    "gwas-catalog:associations",
    "insider:interfaces",
}


def test_every_dataset_without_a_downloader_has_been_classified():
    """A dataset with no downloader must say WHICH kind of no-downloader it is.

    Silence used to mean both "not written yet" and "impossible" (#118). ``byo`` says
    impossible, and buys the implicit download skip, the ``--raw-dir`` pre-flight and
    the instructions. "Not written yet" is spelled by omitting the ``acquisition``
    block (#360), so omission is legitimate only for the datasets pinned in
    ``DOWNLOADER_NOT_WRITTEN_YET``: a new dataset left silently undeclared fails here,
    and so does a pinned one whose downloader has landed (drop it from the set).
    """
    no_downloader = {
        name: ds
        for name, ds in _declared_acquisitions().items()
        if not (ds.get("lifecycle") or {}).get("download")
        and (ds.get("acquisition") or {}).get("mode") != "byo"
    }
    assert set(no_downloader) == DOWNLOADER_NOT_WRITTEN_YET, (
        f"unclassified: {sorted(set(no_downloader) - DOWNLOADER_NOT_WRITTEN_YET)}; "
        f"stale: {sorted(DOWNLOADER_NOT_WRITTEN_YET - set(no_downloader))}"
    )
    # An explicit `mode: download` is a claim that the downloader exists (#360).
    assert not [name for name, ds in no_downloader.items() if "acquisition" in ds]


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
        # An omitted block with NO lifecycle at all is the pre-existing, intentional
        # "not written yet" TODO state (see the 'mode' enum description) -- the #360
        # rule below must not touch it, only the explicit claim.
        (None, None, True),
        # `byo` beside a parse-only lifecycle is coherent: BYO is about acquisition.
        (
            {"mode": "byo", "reason": "license"},
            {"parse": {"module": "m", "function": "f"}},
            True,
        ),
        # The converse (#360): explicit mode 'download' with no lifecycle.download is
        # the claim gevir:metrics and gwas-catalog:associations both made -- accepted
        # here at load time, and only failing once `hvantk reprocess` actually ran it.
        ({"mode": "download"}, None, False),
        ({"mode": "download"}, {"parse": {"module": "m", "function": "f"}}, False),
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

    Also covers the converse rule added for #360: an explicit ``mode: download`` with
    no ``lifecycle.download`` used to validate here and only fail once `hvantk
    reprocess` actually ran the dataset (`reprocess_cli.py`'s "has no lifecycle.download
    declared" `UsageError`) -- the same one-directional hole `byo`/`lifecycle.download`
    had before #351, just facing the other way.
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


def test_validate_explains_the_byo_rule_instead_of_quoting_jsonschema(tmp_path):
    """The translated sentence is the whole point of `_explain`, and was untested.

    The schema rule itself is covered above, but a conditional `if/then` reports as
    `... should not be valid under {'required': [...]}` -- which names the mechanism and
    not the problem. Reverting `_explain(exc)` to `exc.message` passed the entire suite,
    so the message a third-party plugin author actually reads was unguarded.
    """
    from click.testing import CliRunner

    from hvantk.tools.plugins.plugins_cli import plugins_group

    manifest = _manifest_with({"mode": "byo", "reason": "size"})
    manifest["datasets"][0]["lifecycle"] = {
        "download": {"module": "m", "function": "download_dataset"}
    }
    path = tmp_path / "plugin.yaml"
    path.write_text(yaml.safe_dump(manifest))

    result = CliRunner().invoke(plugins_group, ["validate", str(path)])

    assert result.exit_code != 0, result.output
    assert (
        "acquisition.mode is 'byo' but lifecycle.download is declared" in result.output
    ), result.output
    assert "should not be valid under" not in result.output, (
        "raw jsonschema text leaked to the user: " + result.output
    )

    # Converse (#360): explicit mode 'download' with no lifecycle.download gets its
    # own translated sentence too, not the raw `'lifecycle' is a required property`.
    manifest2 = _manifest_with({"mode": "download"})
    path2 = tmp_path / "plugin2.yaml"
    path2.write_text(yaml.safe_dump(manifest2))

    result2 = CliRunner().invoke(plugins_group, ["validate", str(path2)])

    assert result2.exit_code != 0, result2.output
    assert (
        "acquisition.mode is 'download' but no lifecycle.download is declared"
        in result2.output
    ), result2.output
    assert "is a required property" not in result2.output, (
        "raw jsonschema text leaked to the user: " + result2.output
    )
