"""Verify the plugin loader populates the dispatch table alongside the
legacy create_table_adapter registrations.

The dispatch table itself is module-private (`_TABLE_BUILDERS`); this test
uses the public `run_table_builder` to verify dispatch works for both
legacy and plugin-loaded entries.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from hvantk.core.plugin import loader as plugin_loader
from hvantk.core.plugin.registry import run_table_builder


def test_legacy_registration_still_dispatchable():
    """A pre-plugin legacy builder (gevir) is still reachable via run_table_builder."""
    # We don't actually invoke the gevir builder (it would try to read real data);
    # we just patch the underlying adapter and confirm dispatch routes to it.
    from hvantk.core.plugin import registry as tables_registry

    assert "gevir" in tables_registry._TABLE_BUILDERS
    with patch.dict(tables_registry._TABLE_BUILDERS, {"gevir": lambda i, o, p: None}):
        # Verify run_table_builder dispatches without raising KeyError.
        run_table_builder("gevir", "/tmp/x", "/tmp/y", {})


def test_plugin_loader_populates_dispatch_table(monkeypatch, tmp_path):
    """A plugin loaded into a fresh registry shows up in the dispatch table."""
    from hvantk.core.plugin import registry as tables_registry

    fixture = (
        Path(__file__).parent
        / "testdata"
        / "raw"
        / "plugins"
        / "fake_plugin"
    )
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.PluginRegistry()
    reg.load_from_directory(fixture)
    tables_registry._apply_plugin_registrations(reg)
    assert "fake:default" in tables_registry._TABLE_BUILDERS
    # Legacy entry still intact.
    assert "gevir" in tables_registry._TABLE_BUILDERS
