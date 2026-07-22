"""Identifier mapping onto the gene spine, with measurable loss."""
from __future__ import annotations

import pytest

from hvantk.algorithms.annotation.mapping import (
    GeneIdMapper,
    MappingRateError,
    MappingReport,
    enforce_rate,
)


class FakeHGNC:
    """Minimal stand-in for HGNCGeneCatalogStreamer.

    Real HGNC tables need Hail; the mapper's own logic -- alias resolution, spine
    membership, loss accounting -- does not, so it is tested against this fake and the
    real streamer is exercised in Task 5's integration test.
    """

    _CANONICAL = {"MYL7": "MYL7", "OLD1": "NEW1", "NEW1": "NEW1"}
    _TO_HGNC = {"MYL7": "HGNC:7592", "NEW1": "HGNC:00001"}
    _TO_ENSEMBL = {"HGNC:7592": "ENSG00000106631", "HGNC:00001": "ENSG00000000001"}

    def resolve_to_canonical(self, symbol):
        return self._CANONICAL.get(symbol, symbol)

    def map_to_hgnc(self, ids, source_type):
        return {i: self._TO_HGNC.get(i) for i in ids}

    def map_from_hgnc(self, hgnc_ids, target_type):
        assert target_type == "ensembl_gene_id"
        return {h: self._TO_ENSEMBL.get(h) for h in hgnc_ids}


SPINE = {"ENSG00000106631", "ENSG00000000001"}


def test_report_rate_is_mapped_over_input():
    r = MappingReport("x", "symbol", 10, 8, 2, 0, ("A", "B"))
    assert r.rate == 0.8


def test_empty_input_reports_full_rate():
    r = MappingReport("x", "symbol", 0, 0, 0, 0, ())
    assert r.rate == 1.0


def test_symbols_map_through_alias_resolution():
    mapper = GeneIdMapper(FakeHGNC(), SPINE)
    mapping, report = mapper.from_symbols(["MYL7", "OLD1"], source="test")

    assert mapping["MYL7"] == "ENSG00000106631"
    assert mapping["OLD1"] == "ENSG00000000001"   # OLD1 -> NEW1 -> HGNC -> ENSG
    assert report.n_mapped == 2
    assert report.rate == 1.0


def test_unmapped_symbols_are_named_not_just_counted():
    mapper = GeneIdMapper(FakeHGNC(), SPINE)
    mapping, report = mapper.from_symbols(["MYL7", "NOTAGENE"], source="test")

    assert mapping["NOTAGENE"] is None
    assert report.n_unmapped == 1
    assert "NOTAGENE" in report.unmapped


def test_a_gene_absent_from_the_spine_counts_as_unmapped():
    """Mapping to an Ensembl ID the spine does not contain is still a loss."""
    mapper = GeneIdMapper(FakeHGNC(), spine_gene_ids={"ENSG00000106631"})
    mapping, report = mapper.from_symbols(["MYL7", "NEW1"], source="test")

    assert mapping["NEW1"] is None
    assert report.n_unmapped == 1


def test_hgnc_ids_map_directly():
    mapper = GeneIdMapper(FakeHGNC(), SPINE)
    mapping, report = mapper.from_hgnc_ids(["HGNC:7592"], source="test")

    assert mapping["HGNC:7592"] == "ENSG00000106631"
    assert report.key_type == "hgnc_id"


def test_enforce_rate_passes_at_threshold():
    enforce_rate(MappingReport("s", "symbol", 10, 9, 1, 0, ("X",)), min_rate=0.9)


def test_enforce_rate_raises_below_threshold_and_names_examples():
    report = MappingReport("cardoso", "symbol", 10, 5, 5, 0, ("A", "B", "C", "D", "E"))
    with pytest.raises(MappingRateError) as exc:
        enforce_rate(report, min_rate=0.9)

    message = str(exc.value)
    assert "cardoso" in message
    assert "0.50" in message
    assert "A" in message
