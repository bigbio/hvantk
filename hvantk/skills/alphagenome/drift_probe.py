"""Drift probe for alphagenome — documentation-only source (stub).

AlphaGenome is consumed as a live prediction API requiring credentials; there is
no static data file with a stable, programmatically-probeable URL to fingerprint.
This probe returns a structured stub sentinel so ``hvantk drift`` reports a
visible WARNING (status="stub") rather than a silent false-green. Replace with a
real probe if a direct data URL becomes available. See issue #177.
"""
from __future__ import annotations

from hvantk.core.plugin.api import stub_fingerprint

_REASON = (
    "AlphaGenome is a live prediction API requiring credentials; "
    "no static data file to fingerprint"
)


def fetch_fingerprint() -> dict:
    return stub_fingerprint(_REASON)
