"""Plugin registration test. The round-trip test lives in test_builder.py."""

from __future__ import annotations

from hvantk.core.models import AnnotationTable
from hvantk.core.plugin import loader as plugin_loader


def test_pqtl_metrics_registered():
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("pqtl:metrics")
    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "pqtl-v1"


def test_gene_symbol_mapping_resolves_aliases_and_reports_unmapped(caplog):
    from hvantk.skills.pqtl.builder import _map_gene_symbols

    class Catalog:
        def map_ids(self, ids, source_type, target_type):
            assert source_type == "gene_symbol"
            assert target_type == "ensembl_gene_id"
            return {
                symbol: {"CURRENT": "ENSG00000000001"}.get(symbol) for symbol in ids
            }

        def resolve_alias(self, symbol):
            return "CURRENT" if symbol == "OUTDATED" else None

    mapping, unmapped = _map_gene_symbols(Catalog(), {"OUTDATED", "NOT_IN_CATALOG"})

    assert mapping == {"OUTDATED": "ENSG00000000001"}
    assert unmapped == ["NOT_IN_CATALOG"]
    assert "Could not map 1 of 2" in caplog.text
    assert "NOT_IN_CATALOG" in caplog.text
