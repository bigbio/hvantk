"""GTEx eQTL drift probe placeholder.

GTEx ships per-tissue v11 cis-eQTL parquet files behind the GTEx portal
(https://gtexportal.org), and there is no machine-readable manifest of the
release artifacts; a meaningful drift probe needs to track release-version
strings in portal HTML / Google Cloud bucket listings, which is out of scope
for the Phase 1 migration. Returns a stable placeholder until a real probe
is written.
"""

PROBE_VERSION = 0
NOT_IMPLEMENTED_REASON = (
    "gtex-eqtl drift probe not implemented; track upstream v11 release manually"
)


def fetch_fingerprint() -> dict:
    return {
        "probe_version": PROBE_VERSION,
        "source_version": NOT_IMPLEMENTED_REASON,
        "headers": {},
        "checksums": {},
        "fetched_at": "1970-01-01T00:00:00Z",
    }
