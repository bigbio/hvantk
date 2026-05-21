"""Tests for GeneSet: membership, set algebra, immutability."""
from __future__ import annotations

from datetime import datetime, timezone

import pytest

from hvantk.core.models.gene_set import GeneSet
from hvantk.core.models.provenance import Provenance


def _prov():
    return Provenance(
        plugin="t",
        dataset="t:set",
        plugin_version="0.0",
        source_fingerprint="sha256:abc",
        schema_id="t-set-v1",
        build_timestamp=datetime(2026, 5, 20, tzinfo=timezone.utc),
        builder_commit=None,
    )


def test_membership_and_len():
    gs = GeneSet(name="brca", provenance=_prov(), _members=frozenset({"BRCA1", "BRCA2"}))
    assert "BRCA1" in gs
    assert "TP53" not in gs
    assert len(gs) == 2
    assert sorted(iter(gs)) == ["BRCA1", "BRCA2"]


def test_intersection_union_difference():
    a = GeneSet(name="a", provenance=_prov(), _members=frozenset({"X", "Y"}))
    b = GeneSet(name="b", provenance=_prov(), _members=frozenset({"Y", "Z"}))
    assert a.intersection(b).to_set() == {"Y"}
    assert a.union(b).to_set() == {"X", "Y", "Z"}
    assert a.difference(b).to_set() == {"X"}


def test_to_list_and_to_set():
    gs = GeneSet(name="g", provenance=_prov(), _members=frozenset({"B", "A"}))
    assert sorted(gs.to_list()) == ["A", "B"]
    assert gs.to_set() == {"A", "B"}


def test_load_classmethod_round_trip(tmp_path):
    gs = GeneSet(name="brca", provenance=_prov(), _members=frozenset({"BRCA1", "BRCA2"}))
    out = tmp_path / "brca.geneset.json"
    gs.save(out)

    loaded = GeneSet.load(out)
    assert isinstance(loaded, GeneSet)
    assert loaded.provenance == gs.provenance
