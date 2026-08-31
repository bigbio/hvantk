"""pqtl drift probe should fingerprint the upstream preprint's version metadata.

Runs OFFLINE via requests_mock. The pQTL statistics ship as supplementary
material to a preprint, so the publication is the upstream (issue #177 recorded
the absence of a direct data URL correctly).
"""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.pqtl.drift_probe import (
    MEDRXIV_API_URL,
    PQTL_SOURCE_DOI,
    fetch_fingerprint,
)


def _payload(**overrides):
    record = {
        "title": "Regulation of protein abundance in normal human tissues",
        "doi": PQTL_SOURCE_DOI,
        "date": "2025-01-13",
        "version": "1",
        "published": "NA",
    }
    record.update(overrides)
    return {"messages": [{"status": "ok"}], "collection": [record]}


def test_fetch_fingerprint_shape():
    with requests_mock.Mocker() as m:
        m.get(MEDRXIV_API_URL, json=_payload())
        fp = fetch_fingerprint()

    assert fp["source_version"] == "1"
    assert fp["extras"]["published"] == "NA"
    assert fp["extras"]["doi"] == PQTL_SOURCE_DOI


def test_journal_publication_moves_the_checksum():
    """`published` flipping from NA to a journal DOI is the event that most likely
    means the supplementary data no longer matches what an artifact was built from."""
    with requests_mock.Mocker() as m:
        m.get(MEDRXIV_API_URL, json=_payload())
        before = fetch_fingerprint()

    with requests_mock.Mocker() as m:
        m.get(MEDRXIV_API_URL, json=_payload(published="10.1038/s41588-025-0000-0"))
        after = fetch_fingerprint()

    assert before["checksums"] != after["checksums"]


def test_editorial_title_change_does_not_move_the_checksum():
    """Title and abstract are outside the compared surface, so a typo fix must not
    open a pull request carrying no information."""
    with requests_mock.Mocker() as m:
        m.get(MEDRXIV_API_URL, json=_payload())
        before = fetch_fingerprint()

    with requests_mock.Mocker() as m:
        m.get(MEDRXIV_API_URL, json=_payload(title="Regulation of protein abundance"))
        after = fetch_fingerprint()

    assert before["checksums"] == after["checksums"]


def test_newest_version_wins_when_several_are_listed():
    with requests_mock.Mocker() as m:
        m.get(
            MEDRXIV_API_URL,
            json={
                "collection": [
                    {"version": "1", "date": "2025-01-13", "published": "NA",
                     "doi": PQTL_SOURCE_DOI, "title": "t"},
                    {"version": "2", "date": "2025-06-01", "published": "NA",
                     "doi": PQTL_SOURCE_DOI, "title": "t"},
                ]
            },
        )
        fp = fetch_fingerprint()

    assert fp["source_version"] == "2"
    assert fp["informational"]["versions_listed"] == 2


def test_empty_collection_fails_closed():
    with requests_mock.Mocker() as m:
        m.get(MEDRXIV_API_URL, json={"collection": []})
        with pytest.raises(DriftProbeError, match="no records"):
            fetch_fingerprint()
