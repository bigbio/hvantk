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


def test_legacy_entry_skipped_when_plugin_catalog_owns_accession(tmp_path, caplog):
    """Per-plugin catalogs take precedence over the legacy registry.

    Regression guard for F13: prior to the fix, a legacy registry entry whose
    accession overlapped a plugin catalog entry would silently double-insert
    into the omics bucket and ``get_dataset(accession)`` would return whichever
    iteration order happened to land first. The fix logs a WARNING and skips
    the legacy duplicate.
    """
    import json
    import logging

    # Pick an accession that the live plugin layer owns exactly once. Some
    # per-plugin catalogs include intra-file duplicates (a separate
    # data-quality issue); choosing a unique accession isolates the
    # cross-source collision behavior under test.
    plugin_loader.reset_registry_for_tests()
    reg_baseline = HvantkRegistry()
    accession_counts: dict = {}
    for entry in reg_baseline.list_transcriptomics_datasets():
        acc = entry.get("accession")
        if acc:
            accession_counts[acc] = accession_counts.get(acc, 0) + 1
    target_accession = next(
        (acc for acc, n in accession_counts.items() if n == 1), None
    )
    assert target_accession is not None, "no uniquely-owned transcriptomics accession to test against"
    baseline_count = accession_counts[target_accession]
    assert baseline_count == 1

    # Build a fake registry root containing a transcriptomics datasets.json
    # whose entry duplicates the plugin accession.
    fake_root = tmp_path / "registry"
    (fake_root / "transcriptomics").mkdir(parents=True)
    (fake_root / "transcriptomics" / "datasets.json").write_text(
        json.dumps(
            [
                {
                    "accession": target_accession,
                    "title": "LEGACY-DUPLICATE",
                    "description": "this entry should be skipped",
                },
                {
                    "accession": "LEGACY-ONLY-ACCESSION",
                    "title": "Legacy-only entry",
                    "description": "no plugin owns this; should survive",
                },
            ]
        )
    )

    plugin_loader.reset_registry_for_tests()
    with caplog.at_level(logging.WARNING, logger="hvantk.resources.unified_registry"):
        reg = HvantkRegistry(registry_root=fake_root)

    # Per-plugin entry survives; legacy duplicate is dropped.
    matched = [
        e for e in reg.list_transcriptomics_datasets() if e["accession"] == target_accession
    ]
    assert len(matched) == 1, f"expected exactly one entry for {target_accession}, got {len(matched)}"
    assert matched[0].get("title") != "LEGACY-DUPLICATE", (
        "legacy entry leaked through; per-plugin precedence violated"
    )

    # Non-colliding legacy entry is still aggregated.
    assert any(
        e["accession"] == "LEGACY-ONLY-ACCESSION"
        for e in reg.list_transcriptomics_datasets()
    ), "non-colliding legacy entry was incorrectly dropped"

    # The skip emitted a WARNING that names the offending accession.
    assert any(
        target_accession in record.getMessage() and record.levelname == "WARNING"
        for record in caplog.records
    ), "expected WARNING log naming the colliding accession"
