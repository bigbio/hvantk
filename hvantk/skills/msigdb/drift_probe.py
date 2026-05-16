"""MSigDB drift probe placeholder.

MSigDB distributes per-release GMT files at versioned URLs (e.g.
https://www.gsea-msigdb.org/gsea/msigdb/human/download/release/c2.cp.v2024.1.Hs.symbols.gmt).
A meaningful drift probe needs to scrape the release index, which is out of
scope for the Phase 1 migration. Returns a stable placeholder until a real
probe is written.
"""

PROBE_VERSION = 0
NOT_IMPLEMENTED_REASON = (
    "msigdb drift probe not implemented; track upstream releases manually"
)


def fetch_fingerprint() -> dict:
    return {
        "probe_version": PROBE_VERSION,
        "source_version": NOT_IMPLEMENTED_REASON,
        "headers": {},
        "checksums": {},
        "fetched_at": "1970-01-01T00:00:00Z",
    }
