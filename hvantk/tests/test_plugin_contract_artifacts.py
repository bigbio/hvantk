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

from pathlib import Path

import pytest
import yaml

SKILLS_DIR = Path(__file__).resolve().parents[1] / "skills"

ARTIFACT_FIELDS = ("fixture", "schema_snapshot", "row_snapshot", "drift_fingerprint")

# Datasets that do NOT yet ship every declared artifact, with the fields still missing.
# Shrink this as snapshots land; never grow it for a new provider without a reason here.
#
# Three causes, and they are not the same problem:
#   1. no static upstream artifact exists at all -- alphagenome is a credentialed live
#      prediction API, pqtl is publication-only supplementary data. These cannot be
#      snapshotted and are expected to stay listed.
#   2. licence forbids redistributing rows (cosmic-cgc).
#   3. simply never seeded, though the builder runs from a committed fixture -- the
#      majority, and the ones the follow-up work removes from this list.
#
# clingen, gencc and hgnc were removed from this list once their snapshots landed, and
# dbnsfp / ensembl-gene / gnomad-metrics now lack only a drift fingerprint. The
# remaining fingerprint gaps are a separate concern from snapshots: a probe has to be run
# against the live upstream, which the snapshot tests deliberately never touch.
# (gevir shipped its drift fingerprint as part of the gevir plugin-review work, so it
# left this list.)
#
# (dbnsfp / ensembl-gene / gnomad-metrics previously appeared here for a different
# reason -- they declared a plugin-local fixture dir that was never created while the
# tests read one under hvantk/tests/testdata/raw/. Their manifests now point at the real
# shared location, so only their snapshot/fingerprint files remain outstanding.)
#
# uniprot-ptm and both cptac datasets had no committed fixture at all -- their round-trip
# tests synthesized inputs into tmp_path. Small fixtures were committed for each, so all
# three now ship a full set and have left this list.
#
# expression-atlas is the remaining seedable one: its declared fixture dir holds only a
# .gitkeep, and the real inputs under hvantk/tests/testdata/raw/expression_atlas/ are too
# large to use directly (116k transcripts x 320 samples), so a fixture must be *derived*
# by truncation rather than copied.
KNOWN_INCOMPLETE: dict[str, tuple[str, ...]] = {
    "alphagenome:predictions": ARTIFACT_FIELDS,
    "cosmic-cgc:submissions": ARTIFACT_FIELDS,
    "dbnsfp:variants": ("drift_fingerprint",),
    "ensembl-gene:genes": ("drift_fingerprint",),
    "expression-atlas:dataset": ("schema_snapshot", "row_snapshot"),
    "gnomad-metrics:metrics": ("drift_fingerprint",),
    "peptideatlas:phospho": ("schema_snapshot", "row_snapshot"),
    "pqtl:metrics": ARTIFACT_FIELDS,
}


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


def test_no_new_dataset_declares_a_missing_artifact():
    """A dataset may not declare an artifact it does not ship unless it is a known gap."""
    unexpected = {
        name: fields
        for name, fields in _declared_artifact_gaps().items()
        if name not in KNOWN_INCOMPLETE
    }
    assert not unexpected, (
        "These datasets declare validation artifacts that do not exist on disk:\n"
        + "\n".join(f"  {n}: missing {', '.join(f)}" for n, f in sorted(unexpected.items()))
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
