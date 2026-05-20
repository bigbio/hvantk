"""Tests for Provenance dataclass: required fields, unknown() factory, immutability."""
from __future__ import annotations

from datetime import datetime, timezone

import pytest

from hvantk.core.models.provenance import Provenance


def test_provenance_is_frozen():
    p = Provenance(
        plugin="clinvar",
        dataset="clinvar:variants",
        plugin_version="0.2.0",
        source_fingerprint="sha256:abc",
        schema_id="clinvar-variants-v1",
        build_timestamp=datetime(2026, 5, 20, tzinfo=timezone.utc),
        builder_commit="deadbeef",
    )
    with pytest.raises((AttributeError, Exception)):
        p.plugin = "other"  # frozen


def test_provenance_parents_default_empty():
    p = Provenance(
        plugin="clinvar",
        dataset="clinvar:variants",
        plugin_version="0.2.0",
        source_fingerprint="sha256:abc",
        schema_id="clinvar-variants-v1",
        build_timestamp=datetime(2026, 5, 20, tzinfo=timezone.utc),
        builder_commit=None,
    )
    assert p.parents == ()


def test_provenance_unknown_marks_reason():
    p = Provenance.unknown(reason="legacy file, no manifest")
    assert p.plugin == "<unknown>"
    assert p.schema_id == "<unknown>"
    assert "legacy file" in p.source_fingerprint
    assert p.parents == ()


def test_provenance_unknown_requires_reason():
    with pytest.raises(TypeError):
        Provenance.unknown()  # reason is required
