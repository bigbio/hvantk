"""A build must survive an unreachable drift probe.

The probe supplies provenance metadata, not build input. Every
manual-acquisition plugin is built from a hand-staged file, often on a node with
no egress; before the documentation-only plugins gained live probes their probes
were pure in-process calls, so those builds needed no network at all. Letting a
DriftProbeError propagate would make `hvantk reprocess --skip-download
--no-check-drift` fail offline -- and `--no-check-drift` gates only the
post-build check, not the probe call in run_builder.
"""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from hvantk.core.plugin.api import DriftProbeError, PROBE_UNAVAILABLE_TOKEN
from hvantk.core.plugin.run_builder import _coerce_fingerprint, run_builder_for_spec


class _Artifact:
    """Minimal stand-in for a typed artifact: saves, and carries provenance."""

    def __init__(self):
        self.saved_to = None
        self.provenance = SimpleNamespace(schema_id=None)

    def save(self, path):
        self.saved_to = path


def _spec(probe):
    spec = Mock()
    spec.name = "insider:variants"
    spec.artifact_type = _Artifact
    spec.drift_probe = probe
    spec.builder = lambda parsed_input, ctx, **params: _Artifact()
    # Mock auto-creates truthy attributes; the schema_id check would then
    # dereference .provenance on the stand-in artifact.
    spec.schema_id = None
    return spec


def test_unreachable_probe_does_not_abort_the_build(tmp_path, caplog):
    spec = _spec(Mock(side_effect=DriftProbeError("HTTP failure: no route to host")))

    provenance = run_builder_for_spec(
        spec,
        parsed_input="in",
        output_path=str(tmp_path / "out.ht"),
        plugin_version="0.1.0",
    )

    assert provenance is not None
    assert "drift probe could not reach its source" in caplog.text


def test_provenance_never_implies_a_probe_ran(tmp_path):
    """The fallback is self-describing, not a synthesized digest -- the same
    reasoning as STUB_FINGERPRINT_TOKEN."""
    captured = {}

    def builder(parsed_input, ctx, **params):
        captured["fingerprint"] = ctx.source_fingerprint
        return _Artifact()

    spec = _spec(Mock(side_effect=DriftProbeError("unreachable")))
    spec.builder = builder

    run_builder_for_spec(
        spec,
        parsed_input="in",
        output_path=str(tmp_path / "o.ht"),
        plugin_version="0.1.0",
    )

    assert captured["fingerprint"] == PROBE_UNAVAILABLE_TOKEN
    assert not captured["fingerprint"].startswith("sha256:")


def test_a_working_probe_still_stamps_a_real_fingerprint(tmp_path):
    """The degradation must not mask a healthy probe."""
    captured = {}

    def builder(parsed_input, ctx, **params):
        captured["fingerprint"] = ctx.source_fingerprint
        return _Artifact()

    probe_result = {"probe_version": 2, "headers": {"f": {"content_length": "1"}}}
    spec = _spec(Mock(return_value=probe_result))
    spec.builder = builder

    run_builder_for_spec(
        spec,
        parsed_input="in",
        output_path=str(tmp_path / "o.ht"),
        plugin_version="0.1.0",
    )

    assert captured["fingerprint"] == _coerce_fingerprint(probe_result, "x")
    assert captured["fingerprint"] != PROBE_UNAVAILABLE_TOKEN


def test_a_contract_violation_still_raises(tmp_path):
    """Only DriftProbeError degrades; a probe returning a non-dict is a bug and
    must stay loud."""
    spec = _spec(Mock(return_value="not-a-dict"))

    with pytest.raises(Exception):
        run_builder_for_spec(
            spec,
            parsed_input="in",
            output_path=str(tmp_path / "o.ht"),
            plugin_version="0.1.0",
        )
