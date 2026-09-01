"""INSIDER drift probe for ``insider:variants``: HEAD against the genomic BED.

Issue #177 recorded the source as needing a page-scraping probe, because the
download page carries no links in its markup. The underlying path is stable and
directly addressable, so no scraping is required.

This probe covers ``/bed/all.bed`` only -- the >1 GB
``Whole_Human_Interactome_Interface_hg38.bed`` the builder reads. Confirmed to be
that product rather than a same-named sibling: its first bytes are the
``browser hide all`` directive followed by ``track name=A0A0A0MS80_ppi_P56705``
and ``chr11 700235``, matching the committed fixture and the excerpt in
SKILL.md § 4.

``insider:interfaces`` has its own probe and its own baseline
(``interfaces/drift_probe.py``); see ``shared/http_probe.py`` for why the two
must not share one.

Only a HEAD is issued, so the 1.17 GB body is never transferred.
"""

from __future__ import annotations

from hvantk.skills.insider.shared.http_probe import (
    INSIDER_BASE_URL,
    head_fingerprint,
)

INSIDER_BED_URL = f"{INSIDER_BASE_URL}/bed/all.bed"

# The documented product name, not the URL basename (`all.bed`), so the
# fingerprint reads against the name SKILL.md and the catalog use.
INSIDER_BED_FILENAME = "Whole_Human_Interactome_Interface_hg38.bed"


def fetch_fingerprint() -> dict:
    """Fingerprint the live genomic BED (HEAD only)."""
    return head_fingerprint(INSIDER_BED_FILENAME, INSIDER_BED_URL)
