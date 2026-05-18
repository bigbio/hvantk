"""Smoke test: every in-tree plugin loads cleanly and registers expected datasets.

Marked @pytest.mark.hail because plugin loading transitively imports each
builder, which in turn imports hail. This test guards against:
- A plugin.yaml that points at a missing module/function.
- An entry-points + filesystem-scan collision.
- A regression where a builder module fails to import for non-hail reasons.

Locally without hail installed, this test fails at collection. CI's
conda env (environment.yml) installs hail and exercises this.
"""

from __future__ import annotations

import pytest

from hvantk.core import plugin_loader


EXPECTED_PROVIDERS = {
    "clingen",
    "clinvar",
    "cptac",
    "expression-atlas",
    "gencc",
    "gtex-eqtl",
    "gwas-catalog",
    "hgnc",
    "insider",
    "msigdb",
    "peptideatlas",
    "ucsc-cellbrowser",
    "uniprot-ptm",
}

EXPECTED_DATASETS = {
    "clingen:gene-disease",
    "clinvar:variants",
    "cptac:expression",
    "cptac:phospho",
    "expression-atlas:dataset",
    "gencc:submissions",
    "gtex-eqtl:eqtls",
    "gwas-catalog:associations",
    "hgnc:lookup",
    "insider:variants",
    "msigdb:genesets",
    "peptideatlas:phospho",
    "ucsc-cellbrowser:default",
    "ucsc-cellbrowser:adult-ctx",
    "ucsc-cellbrowser:dev-ctx",
    "uniprot-ptm:sites",
}


@pytest.mark.hail
def test_all_plugins_load_without_errors():
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    errors = reg.load_errors()
    assert errors == [], f"Unexpected plugin load errors:\n" + "\n".join(
        f"  {pid}: {exc}" for pid, exc in errors
    )


@pytest.mark.hail
def test_all_expected_providers_registered():
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    registered = {p.name for p in reg.list_providers()}
    missing = EXPECTED_PROVIDERS - registered
    extra = registered - EXPECTED_PROVIDERS
    assert not missing, f"Missing providers: {missing}"
    assert not extra, f"Unexpected extra providers: {extra}"


@pytest.mark.hail
def test_all_expected_datasets_registered():
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    registered = {ds.name for ds in reg.list_datasets()}
    missing = EXPECTED_DATASETS - registered
    extra = registered - EXPECTED_DATASETS
    assert not missing, f"Missing datasets: {missing}"
    assert not extra, f"Unexpected extra datasets: {extra}"
