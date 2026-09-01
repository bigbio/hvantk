"""INSIDER drift probe for ``insider:interfaces``: HEAD against the pair table.

Covers ``/downloads/interfacesALL/H_sapiens_interfacesALL.txt`` (~49 MB) only --
the protein-pair product this dataset builds from, which is a different file from
the genomic BED behind ``insider:variants``.

Deliberately separate from the sibling probe, with its own committed baseline.
The two products are versioned independently upstream (the BED has not moved
since 2018-03-05; this table moved 2024-05-15), and sharing one probe and one
baseline had three consequences: the drift bot writes a ledger row only for the
group's anchor dataset, so an interfaces-only change could never be reported
against the dataset that actually needed rebuilding; a transient fault on the BED
aborted the probe before this file was ever reached, masking real drift here; and
each artifact's recorded provenance covered the other dataset's file as well as
its own.
"""

from __future__ import annotations

from hvantk.skills.insider.shared.http_probe import (
    INSIDER_BASE_URL,
    head_fingerprint,
)

INSIDER_INTERFACES_URL = (
    f"{INSIDER_BASE_URL}/downloads/interfacesALL/H_sapiens_interfacesALL.txt"
)

INSIDER_INTERFACES_FILENAME = "H_sapiens_interfacesALL.txt"


def fetch_fingerprint() -> dict:
    """Fingerprint the live protein-pair interfaces table (HEAD only)."""
    return head_fingerprint(INSIDER_INTERFACES_FILENAME, INSIDER_INTERFACES_URL)
