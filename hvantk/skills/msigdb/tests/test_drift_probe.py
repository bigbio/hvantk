"""MSigDB drift probe: index-page version scrape.

Runs OFFLINE - the MSigDB index endpoint is stubbed via requests_mock.
"""

from __future__ import annotations

import hashlib
import json

import requests
import requests_mock
import pytest

from hvantk.core.plugin_api import DriftProbeError
from hvantk.skills.msigdb.drift_probe import fetch_fingerprint, MSIGDB_INDEX_URL


def test_fetch_fingerprint_shape():
    fake_page = (
        "<html>... C2 CP v2026.1.Hs.symbols.gmt ... and v2026.1.Mm.symbols.gmt"
        " also v2025.1.Hs and v2024.1.Hs ...</html>"
    )
    with requests_mock.Mocker() as m:
        m.get(MSIGDB_INDEX_URL, text=fake_page)
        fp = fetch_fingerprint()

    expected_versions = sorted({"2024.1.Hs", "2025.1.Hs", "2026.1.Hs", "2026.1.Mm"})
    expected_checksum = hashlib.sha256(
        json.dumps(expected_versions, separators=(",", ":")).encode("utf-8")
    ).hexdigest()

    assert fp["probe_version"] == 1
    assert fp["source_version"] == expected_versions[-1]
    assert fp["headers"]["index.jsp"] == ["release_version"]
    assert fp["checksums"]["index.jsp"] == expected_checksum
    assert fp["extras"]["versions_found"] == expected_versions
    assert "fetched_at" in fp


def test_fetch_fingerprint_raises_on_request_failure():
    with requests_mock.Mocker() as m:
        m.get(MSIGDB_INDEX_URL, exc=requests.ConnectionError("boom"))
        with pytest.raises(DriftProbeError):
            fetch_fingerprint()


def test_fetch_fingerprint_raises_when_no_version_tokens():
    with requests_mock.Mocker() as m:
        m.get(MSIGDB_INDEX_URL, text="<html>no version here</html>")
        with pytest.raises(DriftProbeError):
            fetch_fingerprint()
