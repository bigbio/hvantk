"""Drift probe for ensembl-gene:structure.

Unlike the BioMart ``genes`` dataset -- acquired manually, hence a stub probe -- the GTF
has a stable, predictable release URL, so the pinned release is a real fingerprint. A
release bump changes CDS lengths, MANE assignments and coordinates, which is exactly the
kind of upstream movement drift detection exists to surface.
"""
from __future__ import annotations

from hvantk.resources.ensembl_release import (
    ENSEMBL_GTF_URL,
    ENSEMBL_RELEASE,
)


def fetch_fingerprint() -> dict:
    return {
        "release": ENSEMBL_RELEASE,
        "url": ENSEMBL_GTF_URL,
    }
