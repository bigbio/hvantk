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
