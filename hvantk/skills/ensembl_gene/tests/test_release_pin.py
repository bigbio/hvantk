"""The Ensembl release must be defined once and agree everywhere it is used.

The cross-layer check -- that ``hvantk.algorithms.ptm.constants`` re-exports this pin --
lives in ``hvantk/tests/test_ensembl_release_pin.py``, NOT here: a file under ``skills/``
may not import from ``algorithms/`` (enforced by test_dependency_directions), and a
skills-layer test importing the ptm algorithm would itself break that rule.
"""
from __future__ import annotations

import json
from pathlib import Path

from hvantk.resources.ensembl_release import (
    ENSEMBL_GTF_FILENAME,
    ENSEMBL_GTF_URL,
    ENSEMBL_RELEASE,
)

CATALOG = Path("hvantk/skills/ensembl_gene/catalog/datasets.json")


def test_gtf_url_and_filename_carry_the_pinned_release():
    assert ENSEMBL_RELEASE in ENSEMBL_GTF_URL
    assert ENSEMBL_GTF_FILENAME == f"Homo_sapiens.GRCh38.{ENSEMBL_RELEASE}.gtf.gz"
    assert ENSEMBL_GTF_URL.endswith(ENSEMBL_GTF_FILENAME)


def test_catalog_declares_the_same_release():
    entries = json.loads(CATALOG.read_text())
    ensembl = [e for e in entries if e.get("data_source") == "Ensembl"]
    assert ensembl, "no Ensembl entry in the catalog"
    for entry in ensembl:
        assert entry["accession"] == f"Ensembl_v{ENSEMBL_RELEASE}"
        for f in entry.get("files", []):
            if f["path"].endswith(".gtf.gz"):
                assert f["path"] == ENSEMBL_GTF_FILENAME
