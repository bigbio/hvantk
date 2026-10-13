"""Plugin registration and gene-mapping tests. The round trip lives in test_builder.py."""

from __future__ import annotations

from typing import NamedTuple

import pytest

from hvantk.core.models import AnnotationTable
from hvantk.core.plugin import loader as plugin_loader


def test_pqtl_metrics_registered():
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("pqtl:metrics")
    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "pqtl-v1"


class _Gene(NamedTuple):
    hgnc_id: str
    symbol: str
    ensembl_id: str
    previous: tuple = ()
    aliases: tuple = ()


class _FakeHGNCCatalog:
    """In-memory stand-in for the HGNC gene catalog.

    Skills may not import each other, even in tests, so the real
    ``HGNCGeneCatalogStreamer`` is out of reach here. Like the real catalog,
    ``resolve_alias`` merges previous symbols with aliases and returns the first gene
    in hgnc_id order; the builder must not rely on it.
    """

    def __init__(self, *genes):
        self._genes = sorted(genes)

    def map_ids(self, ids, source_type, target_type):
        assert target_type == "ensembl_gene_id"
        field = {"gene_symbol": "symbol", "hgnc_id": "hgnc_id"}[source_type]
        to_ensembl = {getattr(g, field): g.ensembl_id for g in self._genes}
        return {id_: to_ensembl.get(id_) for id_ in ids}

    def is_canonical(self, symbol):
        return any(g.symbol == symbol for g in self._genes)

    def resolve_alias(self, symbol):
        return next(
            (g.symbol for g in self._genes if symbol in g.aliases + g.previous), None
        )

    def get_genes_listing_symbols(self, symbols, field):
        attr = {"prev_symbols": "previous", "alias_symbols": "aliases"}[field]
        found = {}
        for g in self._genes:
            for symbol in set(symbols) & set(getattr(g, attr)):
                found.setdefault(symbol, set()).add(g.hgnc_id)
        return found


def test_gene_symbol_mapping_resolves_outdated_symbols_and_reports_unmapped(caplog):
    from hvantk.skills.pqtl.builder import _map_gene_symbols

    catalog = _FakeHGNCCatalog(
        _Gene("HGNC:1", "CURRENT", "ENSG00000000001", previous=("OUTDATED",)),
        _Gene("HGNC:2", "OTHER", "ENSG00000000002", aliases=("NICKNAME",)),
    )

    mapping = _map_gene_symbols(
        catalog, {"Liver": {"OUTDATED", "NICKNAME", "NOT_IN_CATALOG"}}
    )

    assert mapping == {
        ("Liver", "OUTDATED"): "ENSG00000000001",
        ("Liver", "NICKNAME"): "ENSG00000000002",
    }
    assert (
        "Mapped 2 of 3 pQTL gene symbol(s) to Ensembl gene IDs: 0 approved, "
        "1 via a previous symbol, 1 via an alias; 0 ambiguous, 1 unmapped"
    ) in caplog.text
    assert "NOT_IN_CATALOG" in caplog.text


def test_ambiguous_symbol_stays_unmapped_instead_of_taking_the_first_gene(caplog):
    """Regression: a symbol listed by two genes went to the first in hgnc_id order.

    SHARED is an alias of both genes, so it identifies neither. OLDB is an alias of
    GENEA but the previous symbol of GENEB: previous symbols are tried before aliases,
    so it belongs to GENEB, not to the first gene that lists it.
    """
    from hvantk.skills.pqtl.builder import _map_gene_symbols

    catalog = _FakeHGNCCatalog(
        _Gene("HGNC:1", "GENEA", "ENSG00000000001", aliases=("SHARED", "OLDB")),
        _Gene(
            "HGNC:2",
            "GENEB",
            "ENSG00000000002",
            previous=("OLDB",),
            aliases=("SHARED",),
        ),
    )

    mapping = _map_gene_symbols(catalog, {"Liver": {"GENEA", "SHARED", "OLDB"}})

    assert mapping == {
        ("Liver", "GENEA"): "ENSG00000000001",
        ("Liver", "OLDB"): "ENSG00000000002",
    }
    assert "1 ambiguous, 0 unmapped" in caplog.text
    assert "1 pQTL gene symbol(s) match more than one gene" in caplog.text
    assert "Examples: SHARED" in caplog.text


def test_colliding_symbol_defers_to_the_approved_symbol_only_in_its_tissue(caplog):
    """Regression: two source symbols resolving to the same gene_id duplicated
    (locus, alleles, gene_id) within a tissue.

    OLDNAME is GENEA's previous symbol. In Colon, CURRENT -- GENEA's approved symbol --
    labels rows too, so OLDNAME defers to it there. In Heart, OLDNAME is GENEA's only
    name, so it keeps the gene_id: rows of different tissues are told apart by
    ``tissue``, and a source can name one protein differently in different tissues.
    OLDB1 (a previous symbol) and OLDB2 (an alias) both resolve to GENEB in Colon,
    where no approved symbol arbitrates, so neither keeps the gene_id there.
    """
    from hvantk.skills.pqtl.builder import _map_gene_symbols

    catalog = _FakeHGNCCatalog(
        _Gene("HGNC:1", "CURRENT", "ENSG00000000001", previous=("OLDNAME",)),
        _Gene(
            "HGNC:2",
            "GENEB",
            "ENSG00000000002",
            previous=("OLDB1",),
            aliases=("OLDB2",),
        ),
    )

    mapping = _map_gene_symbols(
        catalog,
        {"Colon": {"CURRENT", "OLDNAME", "OLDB1", "OLDB2"}, "Heart": {"OLDNAME"}},
    )

    assert mapping == {
        ("Colon", "CURRENT"): "ENSG00000000001",
        ("Heart", "OLDNAME"): "ENSG00000000001",
    }
    assert (
        "In tissue Colon, pQTL gene symbol(s) OLDNAME resolve to the same gene "
        "(ENSG00000000001) as the approved symbol 'CURRENT'"
    ) in caplog.text
    assert (
        "In tissue Colon, pQTL gene symbols OLDB1, OLDB2 all resolve to the same gene "
        "(ENSG00000000002)"
    ) in caplog.text
    assert "In tissue Heart" not in caplog.text


def test_no_mapped_symbol_raises():
    from hvantk.skills.pqtl.builder import _map_gene_symbols

    catalog = _FakeHGNCCatalog(_Gene("HGNC:1", "GENEA", "ENSG00000000001"))

    with pytest.raises(ValueError, match="None of the 2 pQTL gene symbols mapped"):
        _map_gene_symbols(catalog, {"Liver": {"UNKNOWN1", "UNKNOWN2"}})


def test_missing_gene_name_does_not_break_the_report(caplog):
    """Regression: a missing gene_name reaches Python as None, which crashed the join."""
    from hvantk.skills.pqtl.builder import _map_gene_symbols

    catalog = _FakeHGNCCatalog(_Gene("HGNC:1", "GENEA", "ENSG00000000001"))

    mapping = _map_gene_symbols(catalog, {"Liver": {"GENEA", "UNKNOWN", None}})

    assert mapping == {("Liver", "GENEA"): "ENSG00000000001"}
    assert "Some pQTL rows have no gene name" in caplog.text
