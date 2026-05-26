"""Stub drift probe returning a deterministic fingerprint."""

FIXED_FINGERPRINT = {
    "probe_version": 1,
    "source_version": "fake-v1",
    "headers": {"a.tsv": ["col1", "col2"]},
    "checksums": {"a.tsv": "deadbeef"},
    "fetched_at": "2026-01-01T00:00:00Z",
}


def fetch_fingerprint() -> dict:
    return dict(FIXED_FINGERPRINT)
