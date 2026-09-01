"""dbnsfp drift probe should fingerprint the advertised release list, not the page."""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.dbnsfp.drift_probe import DBNSFP_LANDING_URL, fetch_fingerprint

_LINKS = [
    "dbNSFP4.8a.zip", "dbNSFP4.8c.zip", "dbNSFP4.9a.zip", "dbNSFP4.9c.zip",
    "dbNSFPv3.5a.zip", "dbNSFPv2.9.3.zip",
]
_PAGE = "<html><body>" + "".join(
    f'<a href="https://dbnsfp.s3.amazonaws.com/{n}">{n}</a>' for n in _LINKS
) + "</body></html>"


def _page(links):
    return "<html><body>" + "".join(
        f'<a href="https://dbnsfp.s3.amazonaws.com/{n}">{n}</a>' for n in links
    ) + "</body></html>"


def test_v_prefixed_releases_are_captured():
    """dbNSFP's 2.x and 3.x generations publish as `dbNSFPv3.5a.zip`. A pattern
    requiring a digit straight after "dbNSFP" matched none of them, so a next
    release named `dbNSFPv5.0a.zip` would have left the projection unchanged and
    reported clean on a real roll-over."""
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text=_PAGE)
        fp = fetch_fingerprint()

    found = fp["headers"]["dbNSFP-release-index"]
    assert "3.5a" in found
    assert "2.9.3" in found


def test_reports_latest_academic_release():
    """`a` is the academic build hvantk uses; `c` is commercial and must not win."""
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text=_PAGE)
        fp = fetch_fingerprint()

    assert fp["source_version"] == "4.9a"


def test_versions_compare_componentwise_not_as_text():
    """4.10a must beat 4.9a. The committed baseline already sits at 4.9a, so the
    very next minor release is the one a text sort gets wrong."""
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text=_page(["dbNSFP4.9a.zip", "dbNSFP4.10a.zip"]))
        fp = fetch_fingerprint()

    assert fp["source_version"] == "4.10a"


def test_projection_ignores_markup_and_link_order():
    """The page re-renders per request (352,830 vs 352,716 bytes on two consecutive
    live fetches), so the compared value must depend only on the release SET."""
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text=_PAGE)
        first = fetch_fingerprint()

    reordered = _page(list(reversed(_LINKS))).replace(
        "<body>", "<body><p>unrelated editorial edit</p>"
    )
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text=reordered)
        second = fetch_fingerprint()

    assert first["headers"] == second["headers"]


def test_case_variation_does_not_invent_a_release():
    """The extraction pattern is case-insensitive, so an uppercase link must fold
    onto the same release rather than entering the set twice."""
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text=_page(["dbNSFP4.9a.zip", "dbNSFP4.9A.zip"]))
        fp = fetch_fingerprint()

    assert fp["headers"]["dbNSFP-release-index"] == ["4.9a"]
    assert fp["source_version"] == "4.9a"


def test_new_release_moves_the_signal():
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text=_PAGE)
        before = fetch_fingerprint()

    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text=_page(_LINKS + ["dbNSFP5.0a.zip"]))
        after = fetch_fingerprint()

    assert before["headers"] != after["headers"]
    assert after["source_version"] == "5.0a"


def test_empty_release_list_fails_closed():
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text="<html><body>no releases here</body></html>")
        with pytest.raises(DriftProbeError, match="advertised no"):
            fetch_fingerprint()


def test_no_academic_release_fails_closed():
    """Silently recording source_version: null would let the bot commit that null,
    after which the probe reports clean having stopped identifying a release."""
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text=_page(["dbNSFP4.9c.zip"]))
        with pytest.raises(DriftProbeError, match="academic"):
            fetch_fingerprint()
