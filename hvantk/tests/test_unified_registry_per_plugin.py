"""Verify HvantkRegistry sources all dataset metadata from per-plugin catalogs.

All omics buckets — transcriptomics, proteomics, and genomics — are now
populated exclusively from per-plugin ``catalog/datasets.json`` files
declared in each plugin's ``plugin.yaml``. There is no legacy file fallback.
A duplicate ``accession`` contributed by two different plugins is a hard
``ValueError``.
"""

from __future__ import annotations

from hvantk.core.plugin import loader as plugin_loader
from hvantk.resources.unified_registry import HvantkRegistry


def _fresh_registry() -> HvantkRegistry:
    # Drop any cached plugin registry so each test sees a freshly built
    # aggregation (mirrors how callers see the registry on first import).
    plugin_loader.reset_registry_for_tests()
    return HvantkRegistry()


def test_transcriptomics_loaded_from_per_plugin_catalogs():
    reg = _fresh_registry()
    entries = reg.list_transcriptomics_datasets()
    # After migration, transcriptomics entries come from the
    # expression-atlas and ucsc-cellbrowser per-plugin catalogs. Before
    # migration there were ~300.
    assert (
        len(entries) > 100
    ), f"expected substantial transcriptomics catalog, got {len(entries)}"


def test_proteomics_loaded_from_per_plugin_catalogs():
    reg = _fresh_registry()
    entries = reg.list_proteomics_datasets()
    # The two HPA E-PROT entries live in the expression-atlas catalog but
    # have proteomics shape (abundance_unit / protein_database / E-PROT-*
    # accession), so the registry routes them to the proteomics bucket.
    assert len(entries) >= 2, f"expected at least 2 proteomics entries, got {len(entries)}"


def test_genomics_loaded_from_per_plugin_catalogs():
    reg = _fresh_registry()
    accs = {e.get("accession") for e in reg.list_genomics_datasets()}
    # The 10 genomics datasets now come from per-plugin catalogs. Two of them
    # (Ensembl_v110, MSigDB_*) come from `mapping`-domain plugins routed into
    # the genomics bucket.
    expected = {
        "dbNSFP_v4.7", "ClinVar_latest", "gnomAD_v4.1", "INSIDER_v1.0",
        "Ensembl_v110", "GeVIR_v1.0", "ClinGen_GeneDisease",
        "GWAS_Catalog_v1.0_e115_r2026-04-27",
        "MSigDB_C2_CP_v2026.1.Hs.symbols", "GTEx_v11_eQTL_signif_pairs",
    }
    assert expected <= accs, f"missing genomics plugin entries: {expected - accs}"
    # CCR has no owning plugin and must NOT appear (it left the catalog).
    assert "CCR_v2.0" not in accs


def test_search_organism_filter_uses_substring_match():
    """`--organism Homo` must match "Homo sapiens" without requiring the canonical form."""
    reg = _fresh_registry()
    partial = reg.search(query="", organism="Homo")
    exact = reg.search(query="", organism="Homo sapiens")
    # Substring `Homo` must include the canonical `Homo sapiens` results
    # (and may also include any other organism whose name contains `Homo`).
    assert len(partial) >= len(exact) > 0
    assert {e.get("accession") for e in exact} <= {e.get("accession") for e in partial}


def test_search_data_source_filter_uses_substring_match():
    """`--data-source UCSC` should match data_source values containing `UCSC` (case-insensitive)."""
    reg = _fresh_registry()
    upper = reg.search(query="", data_source="UCSC")
    lower = reg.search(query="", data_source="ucsc")
    assert len(upper) > 0
    # Case-insensitive: lowercase needle yields the same matches.
    assert {e.get("accession") for e in upper} == {e.get("accession") for e in lower}


def test_search_filters_with_missing_field_are_skipped():
    """Entries with `organism: None` or missing the field must not crash the filter."""
    reg = _fresh_registry()
    # Filtering by a needle that no record has — must return empty without raising.
    results = reg.search(query="", organism="NotAnOrganismValue")
    assert results == []


def test_duplicate_accession_across_plugins_is_error(monkeypatch, tmp_path):
    """Two providers contributing the same accession must raise, not silently skip."""
    import json
    import pytest
    import hvantk.resources.unified_registry as ur

    class _FakeProvider:
        def __init__(self, name, catalog_path):
            self.name = name
            self.catalog_path = str(catalog_path)
            self.primary_domain = "genomics"
            self.datasets = ()

    def _entry(acc):
        return {
            "accession": acc, "title": acc, "description": "x",
            "data_source": "ClinGen", "organism": "Homo sapiens", "files": [],
        }

    cat_a = tmp_path / "a.json"; cat_a.write_text(json.dumps([_entry("DUP")]))
    cat_b = tmp_path / "b.json"; cat_b.write_text(json.dumps([_entry("DUP")]))

    class _FakeRegistry:
        def list_providers(self):
            return [_FakeProvider("prov-a", cat_a), _FakeProvider("prov-b", cat_b)]

    monkeypatch.setattr(
        "hvantk.core.plugin.loader.get_registry", lambda: _FakeRegistry()
    )
    with pytest.raises(ValueError, match="DUP"):
        ur.HvantkRegistry()


def test_unmapped_primary_domain_is_error(monkeypatch, tmp_path):
    """A catalog-owning provider whose primary domain isn't in _DOMAIN_TO_OMICS
    must fail fast, not silently drop its datasets (regression for #185)."""
    import json
    import pytest
    import hvantk.resources.unified_registry as ur

    class _FakeProvider:
        def __init__(self, name, catalog_path, domain):
            self.name = name
            self.catalog_path = str(catalog_path)
            self.primary_domain = domain
            self.datasets = ()

    cat = tmp_path / "c.json"
    cat.write_text(
        json.dumps([{
            "accession": "X1", "title": "X1", "description": "x",
            "data_source": "ClinGen", "organism": "Homo sapiens", "files": [],
        }])
    )

    class _FakeRegistry:
        def list_providers(self):
            return [_FakeProvider("prov-x", cat, "metabolomics")]

    monkeypatch.setattr(
        "hvantk.core.plugin.loader.get_registry", lambda: _FakeRegistry()
    )
    with pytest.raises(ValueError, match="metabolomics"):
        ur.HvantkRegistry()


