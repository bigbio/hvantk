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

from hvantk.core.plugin import loader as plugin_loader


# Deliberately hard-coded rather than derived from the manifests on disk: this is a
# tripwire, so adding or removing a provider MUST be a conscious edit here. It is only
# useful if it actually runs -- these tests are `hail`-marked, which pytest.ini deselects
# by default, so the list silently rotted 8 providers behind the tree (alphagenome,
# cosmic-cgc, dbnsfp, ensembl-gene, gevir, gnomad-metrics, onek-genomes, pqtl) until the
# CI job added alongside this change started exercising the marker.
EXPECTED_PROVIDERS = {
    "alphagenome",
    "clingen",
    "clinvar",
    "cosmic-cgc",
    "cptac",
    "dbnsfp",
    "ensembl-gene",
    "expression-atlas",
    "gencc",
    "gevir",
    "gnomad-metrics",
    "gtex-eqtl",
    "gwas-catalog",
    "hgnc",
    "insider",
    "msigdb",
    "onek-genomes",
    "peptideatlas",
    "pqtl",
    "ucsc-cellbrowser",
    "uniprot-ptm",
}

EXPECTED_DATASETS = {
    "alphagenome:predictions",
    "clingen:gene-disease",
    "clinvar:variants",
    "cosmic-cgc:submissions",
    "cptac:expression",
    "cptac:phospho",
    "dbnsfp:variants",
    "ensembl-gene:genes",
    "ensembl-gene:structure",
    "expression-atlas:dataset",
    "gencc:submissions",
    "gevir:metrics",
    "gnomad-metrics:metrics",
    "gtex-eqtl:eqtls",
    "gwas-catalog:associations",
    "hgnc:lookup",
    "insider:variants",
    "msigdb:genesets",
    "onek-genomes:samples",
    "onek-genomes:variants",
    "peptideatlas:phospho",
    "pqtl:metrics",
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
