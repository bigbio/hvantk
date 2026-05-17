"""UCSC Cell Browser drift probe placeholder.

The UCSC Cell Browser ships data per collection (``adultPancreas``,
``cortex-dev``, ``hoc/all-heart``, ...), each with its own release cadence
and per-file MD5s. A single "fingerprint of UCSC Cell Browser" does not
carry meaningful semantics: drift detection must happen at the
collection level (per-collection ``dataset.json`` ``md5`` plus the
expression matrix and metadata file ``Last-Modified``/``Content-Length``
HEAD), and the configured list of tracked collections lives in
``hvantk/resources/cells_ucsc_datasets.json`` rather than in this
provider.

Building that per-collection probe is out of scope for the Phase 1
migration -- it requires the registry/catalog wiring to enumerate which
collections are pinned. Returns a stable placeholder until a real probe
is written.
"""

PROBE_VERSION = 0
NOT_IMPLEMENTED_REASON = (
    "UCSC Cell Browser is per-collection; manual version tracking via "
    "cells_ucsc_datasets.json"
)


def fetch_fingerprint() -> dict:
    return {
        "probe_version": PROBE_VERSION,
        "source_version": NOT_IMPLEMENTED_REASON,
        "headers": {},
        "checksums": {},
        "fetched_at": "1970-01-01T00:00:00Z",
    }
