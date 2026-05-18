"""Verify HvantkRegistry sources data from per-plugin catalogs after migration.

Before the migration the registry read three flat files:
- registry/transcriptomics/datasets.json   (~300 entries)
- registry/proteomics/datasets.json        (2 entries)
- registry/epigenomics/datasets.json       (0 entries)

After the migration those entries live inside each owning plugin's
``catalog/datasets.json`` (expression-atlas, ucsc-cellbrowser). The legacy
``registry/genomics/datasets.json`` is intentionally retained until orphan
entries (dbNSFP, gnomad-metrics, ensembl-gene, gevir, cosmic-cgc) gain
owning plugins.
"""

from __future__ import annotations

from hvantk.core import plugin_loader
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


def test_genomics_still_loaded_from_legacy_path():
    reg = _fresh_registry()
    entries = reg.list_genomics_datasets()
    # The legacy registry/genomics/datasets.json was NOT migrated; should
    # still have its in-tree entries.
    assert len(entries) >= 5


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
