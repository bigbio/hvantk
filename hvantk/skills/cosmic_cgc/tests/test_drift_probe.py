"""cosmic-cgc drift probe should read the release index, never page prose."""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.cosmic_cgc.drift_probe import (
    COSMIC_RELEASE_NOTES_URL,
    fetch_fingerprint,
)

# Anchors are the real index; the prose deliberately names OTHER product versions,
# which is what the live page does ("COSMIC v20 of the Actionability data").
_PAGE = """
<html><body>
  <h2 id="v104">COSMIC v104</h2>
  <p>COSMIC v104 is combined with COSMIC v20 of the Actionability data.</p>
  <h2 id="v103">COSMIC v103</h2>
  <p>COSMIC v18 of the Actionability data are released.</p>
  <h2 id="v102">COSMIC v102</h2>
</body></html>
"""


def test_prose_versions_are_excluded_from_the_index():
    """v20/v18 belong to the Actionability product, a different version series.
    Matching them produced a non-contiguous 'release list' and would open a no-op
    PR on any editorial edit naming an old release."""
    with requests_mock.Mocker() as m:
        m.get(COSMIC_RELEASE_NOTES_URL, text=_PAGE)
        fp = fetch_fingerprint()

    found = fp["headers"]["cosmic-release-index"]
    assert found == ["v102", "v103", "v104"]
    assert "v20" not in found and "v18" not in found


def test_editorial_prose_edit_does_not_move_the_signal():
    with requests_mock.Mocker() as m:
        m.get(COSMIC_RELEASE_NOTES_URL, text=_PAGE)
        before = fetch_fingerprint()

    with requests_mock.Mocker() as m:
        m.get(
            COSMIC_RELEASE_NOTES_URL,
            text=_PAGE.replace("</body>", "<p>Unchanged since COSMIC v95.</p></body>"),
        )
        after = fetch_fingerprint()

    assert before["headers"] == after["headers"]


def test_forward_looking_prose_cannot_bump_source_version():
    """'coming in COSMIC v105' must not report a release that does not exist."""
    with requests_mock.Mocker() as m:
        m.get(
            COSMIC_RELEASE_NOTES_URL,
            text=_PAGE.replace("</body>", "<p>Coming soon: COSMIC v105.</p></body>"),
        )
        fp = fetch_fingerprint()

    assert fp["source_version"] == "v104"


def test_versions_sort_numerically():
    with requests_mock.Mocker() as m:
        m.get(COSMIC_RELEASE_NOTES_URL, text='<i id="v99"></i><i id="v104"></i>')
        fp = fetch_fingerprint()

    assert fp["source_version"] == "v104"
    assert fp["headers"]["cosmic-release-index"] == ["v99", "v104"]


def test_new_release_moves_the_signal():
    with requests_mock.Mocker() as m:
        m.get(COSMIC_RELEASE_NOTES_URL, text=_PAGE)
        before = fetch_fingerprint()

    with requests_mock.Mocker() as m:
        m.get(COSMIC_RELEASE_NOTES_URL, text=_PAGE + '<h2 id="v105">COSMIC v105</h2>')
        after = fetch_fingerprint()

    assert before["headers"] != after["headers"]
    assert after["source_version"] == "v105"


def test_login_redirect_fails_closed():
    """An actual redirect, not just a login-shaped body: the host is known to 302
    the trailing-slash path to /cosmic/login, and a login page answers 200."""
    with requests_mock.Mocker() as m:
        m.get(
            COSMIC_RELEASE_NOTES_URL,
            status_code=302,
            headers={"Location": "https://cancer.sanger.ac.uk/cosmic/login"},
        )
        m.get("https://cancer.sanger.ac.uk/cosmic/login", text="<p>Please log in</p>")
        with pytest.raises(DriftProbeError, match="redirected"):
            fetch_fingerprint()


def test_missing_anchors_fail_closed():
    with requests_mock.Mocker() as m:
        m.get(COSMIC_RELEASE_NOTES_URL, text="<html><body>Please log in</body></html>")
        with pytest.raises(DriftProbeError, match="no 'id="):
            fetch_fingerprint()
