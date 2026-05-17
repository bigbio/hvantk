"""Expression Atlas drift probe placeholder.

Expression Atlas hosts many per-accession experiments (e.g. ``E-MTAB-5214``,
``E-GTEX-8``), each with its own release cadence and per-file checksums. A
single "fingerprint of Expression Atlas" doesn't carry meaningful semantics:
drift detection must happen at the accession level (per-experiment SDRF /
TPM file ``Last-Modified`` + ``Content-Length`` HEAD), and the configured
list of tracked accessions lives in
``hvantk/resources/registry/transcriptomics/datasets.json`` rather than in
this provider.

Building that per-accession probe is out of scope for the Phase 1 migration
- it requires the registry / catalog wiring to enumerate which accessions
are pinned. Returns a stable placeholder until a real probe is written.
"""

PROBE_VERSION = 0
NOT_IMPLEMENTED_REASON = (
    "expression-atlas drift probe not implemented; track per-accession "
    "releases manually via registry/transcriptomics/datasets.json"
)


def fetch_fingerprint() -> dict:
    return {
        "probe_version": PROBE_VERSION,
        "source_version": NOT_IMPLEMENTED_REASON,
        "headers": {},
        "checksums": {},
        "fetched_at": "1970-01-01T00:00:00Z",
    }
