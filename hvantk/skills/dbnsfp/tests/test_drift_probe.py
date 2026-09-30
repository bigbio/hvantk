"""dbnsfp drift probe should fingerprint the advertised release list, not the page."""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.dbnsfp.drift_probe import DBNSFP_RELEASES_URL, fetch_fingerprint

# Shaped like the live https://www.dbnsfp.org/releases/ page (verified 2026-09-30):
# a shared version normally covers both branches (no suffix), and one past release,
# v5.1.1c, shipped commercial-only with the branch suffix glued onto the version and
# no companion "5.1.1a" entry.
_RELEASES = [
    ("5.4", "August 1, 2026"),
    ("5.3.1", "January 1, 2026"),
    ("5.3", "October 6, 2025"),
    ("5.2", "July 2, 2025"),
    ("5.1.1c", "April 24, 2025"),
    ("5.1", "March 21, 2025"),
    ("5.0", "January 1, 2025"),
]


def _page(releases):
    """Minimal page shaped like the real one: an h1/h2 shell plus one
    ``release-title`` span per entry. The probe does not distinguish the "Current
    Release" section from "Past Releases", so this helper does not bother either.
    """
    items = "".join(
        '<div class="release-item"><p class="release-heading">'
        f'<span class="release-title">dbNSFP v{version} ({date})</span>'
        "</p></div>"
        for version, date in releases
    )
    return f"<html><body><h1>Current Release</h1>{items}</body></html>"


_PAGE = _page(_RELEASES)


def test_reports_latest_academic_release():
    """The newest release carrying an academic build wins `source_version`."""
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_RELEASES_URL, text=_PAGE)
        fp = fetch_fingerprint()

    assert fp["source_version"] == "5.4"


def test_commercial_only_patch_does_not_win_latest_academic():
    """A commercial-only patch (the live v5.1.1c pattern) must not be reported as
    the latest academic release even when it is the newest entry on the page --
    there is no academic file behind it."""
    with requests_mock.Mocker() as m:
        m.get(
            DBNSFP_RELEASES_URL,
            text=_page([("5.4", "August 1, 2026"), ("5.4.1c", "September 1, 2026")]),
        )
        fp = fetch_fingerprint()

    assert fp["source_version"] == "5.4"
    assert any("5.4.1c" in h for h in fp["headers"]["dbNSFP-release-index"])


def test_commercial_suffix_case_variation_is_still_excluded():
    """The suffix check is case-insensitive, so an uppercase 'C' must not be
    mistaken for an unrecognised (and therefore eligible) token."""
    with requests_mock.Mocker() as m:
        m.get(
            DBNSFP_RELEASES_URL,
            text=_page([("5.4", "August 1, 2026"), ("5.4.1C", "September 1, 2026")]),
        )
        fp = fetch_fingerprint()

    assert fp["source_version"] == "5.4"


def test_versions_compare_componentwise_not_as_text():
    """5.10 must beat 5.4. The committed baseline sits in the 5.x series, so the
    first double-digit minor release is the one a text sort gets wrong."""
    with requests_mock.Mocker() as m:
        m.get(
            DBNSFP_RELEASES_URL,
            text=_page([("5.4", "d1"), ("5.10", "d2")]),
        )
        fp = fetch_fingerprint()

    assert fp["source_version"] == "5.10"


def test_projection_ignores_markup_and_link_order():
    """Reordering the release entries, or adding unrelated markup, must not move
    the fingerprint -- only the release SET matters."""
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_RELEASES_URL, text=_PAGE)
        first = fetch_fingerprint()

    reordered = _page(list(reversed(_RELEASES))).replace(
        "<body>", "<body><p>unrelated editorial edit</p>"
    )
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_RELEASES_URL, text=reordered)
        second = fetch_fingerprint()

    assert first["headers"] == second["headers"]


def test_new_release_moves_the_signal():
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_RELEASES_URL, text=_PAGE)
        before = fetch_fingerprint()

    with requests_mock.Mocker() as m:
        m.get(DBNSFP_RELEASES_URL, text=_page(_RELEASES + [("5.5", "d")]))
        after = fetch_fingerprint()

    assert before["headers"] != after["headers"]
    assert after["source_version"] == "5.5"


def test_empty_release_list_fails_closed():
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_RELEASES_URL, text="<html><body>no releases here</body></html>")
        with pytest.raises(DriftProbeError, match="no release-title"):
            fetch_fingerprint()


def test_unrecognised_release_format_fails_closed():
    with requests_mock.Mocker() as m:
        m.get(
            DBNSFP_RELEASES_URL,
            text=(
                '<html><body><span class="release-title">'
                "dbNSFP version five point four"
                "</span></body></html>"
            ),
        )
        with pytest.raises(DriftProbeError, match="did not match"):
            fetch_fingerprint()


def test_no_academic_release_fails_closed():
    """Silently recording source_version: null would let the bot commit that null,
    after which the probe reports clean having stopped identifying a release."""
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_RELEASES_URL, text=_page([("5.4.1c", "September 1, 2026")]))
        with pytest.raises(DriftProbeError, match="academic"):
            fetch_fingerprint()
