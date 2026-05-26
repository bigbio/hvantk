"""Drift probe stub for gevir.

Phase K plugin promotion lifts this source from the legacy hardcoded
_TABLE_BUILDERS dict. A real drift probe (which queries the upstream
source and returns a fingerprint) is a follow-up — for now this stub
returns a placeholder fingerprint.
"""
from __future__ import annotations


def fetch_fingerprint() -> dict:
    return {
        "fingerprint": "sha256:phase-k-stub-not-implemented",
        "probe_status": "stub",
        "comment": "Phase K plugin promotion stub; replace with real probe.",
    }
