"""GWAS Catalog drift probe placeholder.

The NHGRI-EBI GWAS Catalog publishes a stable URL for the latest full
associations TSV
(https://www.ebi.ac.uk/gwas/api/search/downloads/full); a HEAD request
against that URL could yield an ETag / Last-Modified fingerprint, and the
catalog itself exposes a release-name string (e.g. ``e0_r2024-09-01``) on
the GWAS portal that a richer probe could surface. Writing a real probe is
out of scope for the Phase 1 migration; returns a stable placeholder until
the real one lands.
"""

PROBE_VERSION = 0
NOT_IMPLEMENTED_REASON = (
    "gwas-catalog drift probe not implemented; track upstream release manually"
)


def fetch_fingerprint() -> dict:
    return {
        "probe_version": PROBE_VERSION,
        "source_version": NOT_IMPLEMENTED_REASON,
        "headers": {},
        "checksums": {},
        "fetched_at": "1970-01-01T00:00:00Z",
    }
