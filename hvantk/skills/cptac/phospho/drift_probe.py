"""CPTAC phosphoproteomics drift probe.

Thin wrapper around :mod:`hvantk.skills.cptac.shared.drift`. The CPTAC
package version + upstream GitHub-release tag together form the drift
fingerprint; see the shared helper for details.
"""

from __future__ import annotations

from hvantk.skills.cptac.shared.drift import (
    PROBE_VERSION,
    cptac_version as _cptac_version,
    fetch_fingerprint as _shared_fetch_fingerprint,
)

__all__ = ["PROBE_VERSION", "fetch_fingerprint", "_cptac_version"]

_DATASET_NAME = "phosphoproteomics"


def fetch_fingerprint() -> dict:
    """Fingerprint of the CPTAC phosphoproteomics dataset shipping."""
    return _shared_fetch_fingerprint(_DATASET_NAME)
