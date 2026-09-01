"""Every committed drift baseline must be real probe output.

`test_plugin_contract_artifacts` checks only that the file EXISTS. That cannot
support the claim the ledger leans on -- that a baseline "was captured from a live
probe run" -- and it is exactly the gap `placeholder_baseline_reason` was written
to close after six hand-seeded baselines produced a permanent, meaningless
`drifted` verdict while looking maximally alive.

A hand-edit (a merge-conflict resolution, a redaction, a copy-paste of another
plugin's shape) re-creates those dead comparators, and the existence ratchet would
still report the dataset complete. These tests read the committed JSON directly, so
they need no network and run in the default selection.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from hvantk.core.plugin.api import (
    PROBE_STATUS_STUB,
    placeholder_baseline_reason,
)

SKILLS_DIR = Path(__file__).resolve().parents[1] / "skills"


def _committed_baselines() -> dict[str, Path]:
    """Map ``provider:dataset`` -> the drift baseline its manifest declares."""
    found: dict[str, Path] = {}
    for manifest_path in sorted(SKILLS_DIR.glob("*/plugin.yaml")):
        manifest = yaml.safe_load(manifest_path.read_text())
        provider = manifest.get("name") or manifest_path.parent.name
        for dataset in manifest.get("datasets", []):
            rel = (dataset.get("tests") or {}).get("drift_fingerprint")
            if not rel:
                continue
            path = manifest_path.parent / rel
            if path.exists():
                found[f"{provider}:{dataset['name']}"] = path
    return found


BASELINES = _committed_baselines()
assert BASELINES, "no committed drift baselines discovered"


@pytest.mark.parametrize("name", sorted(BASELINES))
def test_baseline_is_valid_json_object(name):
    payload = json.loads(BASELINES[name].read_text())
    assert isinstance(payload, dict), f"{name}: baseline is not a JSON object"


@pytest.mark.parametrize("name", sorted(BASELINES))
def test_baseline_is_not_hand_seeded(name):
    """The check drift_runner performs before diffing, applied at rest.

    Without this a placeholder or empty-digest baseline can be committed and only
    surfaces later as a permanent probe_failed that regenerating cannot clear.
    """
    payload = json.loads(BASELINES[name].read_text())
    reason = placeholder_baseline_reason(payload)
    assert reason is None, f"{name}: {reason}"


@pytest.mark.parametrize("name", sorted(BASELINES))
def test_baseline_records_a_real_fetch(name):
    """A captured baseline carries the timestamp of the run that produced it."""
    payload = json.loads(BASELINES[name].read_text())
    if payload.get("probe_status") == PROBE_STATUS_STUB:
        pytest.skip("stub sentinel: no fetch to record")
    fetched_at = payload.get("fetched_at")
    assert isinstance(fetched_at, str) and fetched_at, f"{name}: no fetched_at"
    assert not fetched_at.startswith("1970-01-01"), f"{name}: epoch fetched_at"


@pytest.mark.parametrize("name", sorted(BASELINES))
def test_baseline_carries_a_comparable_signal(name):
    """A baseline whose entire compared surface is empty compares equal forever.

    `placeholder_baseline_reason` deliberately does not flag an empty `checksums`
    map (peptideatlas legitimately ships one), so nothing else catches the case
    where `headers`, `checksums` AND `source_version` are all empty at once.
    """
    payload = json.loads(BASELINES[name].read_text())
    if payload.get("probe_status") == PROBE_STATUS_STUB:
        pytest.skip("stub sentinel: no comparable signal by design")
    signal = (
        payload.get("headers")
        or payload.get("checksums")
        or payload.get("extras")
        or payload.get("source_version")
    )
    assert signal, (
        f"{name}: baseline carries no compared signal at all; every future run "
        "would report clean regardless of what upstream does"
    )
