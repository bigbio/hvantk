"""Drift probe for cosmic-cgc — documentation-only source (stub).

The COSMIC Cancer Gene Census is behind a login/license gate
(cancer.sanger.ac.uk/census); there is no public direct URL to fingerprint.
This probe returns a structured stub sentinel so ``hvantk drift`` reports a
visible WARNING (status="stub") rather than a silent false-green. Replace with a
real probe if a direct data URL becomes available. See issue #177.
"""
from __future__ import annotations

from hvantk.core.plugin.api import stub_fingerprint

_REASON = (
    "COSMIC Cancer Gene Census is login/license-gated "
    "(cancer.sanger.ac.uk/census); no public direct URL to fingerprint"
)


def fetch_fingerprint() -> dict:
    return stub_fingerprint(_REASON)
