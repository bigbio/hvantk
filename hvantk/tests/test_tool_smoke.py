"""Smoke test: every hvantk/tools/**/*.tool.yaml validates and loads cleanly."""

from __future__ import annotations

from hvantk.core import tool_loader


def test_all_tool_manifests_load_without_errors():
    tool_loader.reset_registry_for_tests()
    reg = tool_loader.get_registry()
    errors = reg.load_errors()
    assert errors == [], (
        "Unexpected tool manifest load errors:\n"
        + "\n".join(f"  {p}: {e}" for p, e in errors)
    )


def test_all_expected_tool_domains_have_at_least_one_tool():
    """Smoke: every domain in the Literal enum should have at least one manifest."""
    tool_loader.reset_registry_for_tests()
    reg = tool_loader.get_registry()
    domains = {t.domain for t in reg.list_tools()}
    expected = {
        "plugins", "build", "expression", "annotation", "genesets",
        "ptm", "ancestry", "qtl", "enrichex", "hgc", "infra",
    }
    missing = expected - domains
    assert not missing, f"Domains with no manifested tool: {missing}"
