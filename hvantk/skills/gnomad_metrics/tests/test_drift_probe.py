"""gnomad-metrics drift probe should fingerprint every constraint table offline.

Runs OFFLINE via requests_mock, so CI never hits the gnomAD bucket. The probe
replaced a stub sentinel (issue #177, which rated the source "marginally
feasible") once the public GCS objects were confirmed to answer a HEAD with an
MD5 ETag.
"""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.gnomad_metrics.drift_probe import _object_paths, fetch_fingerprint
from hvantk.skills.gnomad_metrics.shared.constants import GNOMAD_RELEASE_BASE_URL

_BY_GENE = "2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"


def _mock_all(m, *, headers=None):
    for path in _object_paths():
        m.head(
            f"{GNOMAD_RELEASE_BASE_URL}/{path}",
            headers=headers
            if headers is not None
            else {
                "ETag": '"967fb59dc509d9ee351b9d35efacf52a"',
                "Content-Length": "4609488",
                "x-goog-generation": "1597936745580178",
                "Last-Modified": "Thu, 20 Aug 2020 15:19:05 GMT",
            },
        )


def test_every_declared_constraint_table_is_fingerprinted():
    """All declared tables are probed, not just the default, so the fingerprint
    covers whichever release a given build pinned."""
    with requests_mock.Mocker() as m:
        _mock_all(m)
        fp = fetch_fingerprint()

    assert set(fp["checksums"]) == set(_object_paths())
    assert len(_object_paths()) == 3
    assert fp["checksums"][_BY_GENE] == "967fb59dc509d9ee351b9d35efacf52a"


def test_generation_is_compared_so_a_silent_rewrite_is_visible():
    """x-goog-generation changes on every object rewrite even when the bytes are
    identical, so it belongs in the compared surface."""
    with requests_mock.Mocker() as m:
        _mock_all(m)
        fp = fetch_fingerprint()

    assert fp["extras"][_BY_GENE]["generation"] == "1597936745580178"
    assert fp["extras"][_BY_GENE]["content_length"] == "4609488"


def test_last_modified_is_demoted_out_of_the_compared_surface():
    with requests_mock.Mocker() as m:
        _mock_all(m)
        fp = fetch_fingerprint()

    assert fp["source_version"] is None
    assert (
        fp["informational"][_BY_GENE]["last_modified"]
        == "Thu, 20 Aug 2020 15:19:05 GMT"
    )


def test_no_content_signal_fails_closed():
    with requests_mock.Mocker() as m:
        _mock_all(m, headers={})
        with pytest.raises(DriftProbeError, match="no content signal"):
            fetch_fingerprint()
