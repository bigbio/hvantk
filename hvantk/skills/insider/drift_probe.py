"""INSIDER drift probe placeholder.

The Interactome Insider portal (http://interactomeinsider.yulab.org) hosts
the ``Whole_Human_Interactome_Interface_hg38.bed`` file behind a manual
download page with no machine-readable manifest or release feed; a real
drift probe would need to track the page contents or the file's
Last-Modified/ETag header. Writing that probe is out of scope for the
Phase 1 migration; returns a stable placeholder until the real one lands.
"""

PROBE_VERSION = 0
NOT_IMPLEMENTED_REASON = (
    "insider drift probe not implemented; track upstream release manually"
)


def fetch_fingerprint() -> dict:
    return {
        "probe_version": PROBE_VERSION,
        "source_version": NOT_IMPLEMENTED_REASON,
        "headers": {},
        "checksums": {},
        "fetched_at": "1970-01-01T00:00:00Z",
    }
