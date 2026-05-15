"""Verify that the plugin loader populates TABLE_BUILDERS / MATRIX_BUILDERS
alongside the legacy create_table_adapter registrations."""

from __future__ import annotations

from pathlib import Path

from hvantk.core import plugin_loader


def test_legacy_registrations_still_present():
    from hvantk.tables.registry import TABLE_BUILDERS
    # Sanity: an existing legacy builder is still registered.
    assert "clinvar" in TABLE_BUILDERS


def test_plugin_loader_populates_table_builders(monkeypatch, tmp_path):
    # Point the loader at the fake plugin fixture instead of the real skills dir.
    from hvantk.tables import registry as tables_registry

    fixture = (
        Path(__file__).parent
        / "testdata"
        / "raw"
        / "plugins"
        / "fake_plugin"
    )
    plugin_loader.reset_registry_for_tests()
    # Build a registry seeded only with the fake plugin.
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(fixture)
    # Apply that registry to the legacy dict (function added in this task).
    tables_registry._apply_plugin_registrations(reg)
    assert "fake:default" in tables_registry.TABLE_BUILDERS
    # Legacy entry still intact.
    assert "clinvar" in tables_registry.TABLE_BUILDERS
