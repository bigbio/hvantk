"""Sanity test for the PeptideAtlas phospho drift probe.

Network is mocked: the probe shape (probe_version + headers keyed by the
pinned build's zip filename) is what we assert, not the live HTTP response.
"""

from unittest import mock


def test_fingerprint_shape():
    from hvantk.skills.peptideatlas.phospho import drift_probe

    class _FakeResp:
        headers = {
            "Last-Modified": "Wed, 01 Jan 2025 00:00:00 GMT",
            "Content-Length": "12345",
        }

        def raise_for_status(self):
            return None

    with mock.patch("requests.head", return_value=_FakeResp()):
        fp = drift_probe.fetch_fingerprint()

    assert fp["probe_version"] == drift_probe.PROBE_VERSION
    assert fp["source_version"] == "Wed, 01 Jan 2025 00:00:00 GMT"
    # Exactly one filename entry; keyed by the pinned build's zip filename.
    assert len(fp["headers"]) == 1
    (entry,) = fp["headers"].values()
    assert entry["content_length"] == "12345"
    assert entry["build_date"]
    assert entry["build_id"]
