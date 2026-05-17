"""Sanity test for the placeholder UCSC Cell Browser drift probe."""

from hvantk.skills.ucsc_cellbrowser.drift_probe import fetch_fingerprint


def test_placeholder_fingerprint_shape():
    fp = fetch_fingerprint()
    assert fp["probe_version"] == 0
    assert "per-collection" in fp["source_version"]
