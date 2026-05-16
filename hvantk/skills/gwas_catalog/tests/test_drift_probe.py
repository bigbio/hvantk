"""Sanity test for the placeholder drift probe."""

from hvantk.skills.gwas_catalog.drift_probe import fetch_fingerprint


def test_placeholder_fingerprint_shape():
    fp = fetch_fingerprint()
    assert fp["probe_version"] == 0
    assert "not implemented" in fp["source_version"]
