"""cosmic-cgc drift probe should fingerprint the public release list.

Runs OFFLINE via requests_mock. The data stays licence-gated (issue #177 was
right about that); what this probe reaches is the release notes page, which
COSMIC serves without a login.
"""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.cosmic_cgc.drift_probe import (
    COSMIC_RELEASE_NOTES_URL,
    fetch_fingerprint,
)

_PAGE = """
<html><body>
  <h2>COSMIC v104</h2><p>Release v104 notes</p>
  <h2>COSMIC v103</h2>
  <h2>COSMIC v102</h2>
</body></html>
"""


def test_reports_the_newest_release():
    with requests_mock.Mocker() as m:
        m.get(COSMIC_RELEASE_NOTES_URL, text=_PAGE)
        fp = fetch_fingerprint()

    assert fp["source_version"] == "v104"
    assert fp["extras"]["releases_found"] == ["v102", "v103", "v104"]


def test_versions_sort_numerically_not_lexicographically():
    """v104 must beat v99; a string sort would pick the wrong newest release."""
    with requests_mock.Mocker() as m:
        m.get(
            COSMIC_RELEASE_NOTES_URL,
            text="<p>COSMIC v99</p><p>COSMIC v104</p>",
        )
        fp = fetch_fingerprint()

    assert fp["source_version"] == "v104"


def test_new_release_moves_the_checksum():
    with requests_mock.Mocker() as m:
        m.get(COSMIC_RELEASE_NOTES_URL, text=_PAGE)
        before = fetch_fingerprint()

    with requests_mock.Mocker() as m:
        m.get(COSMIC_RELEASE_NOTES_URL, text=_PAGE + "<h2>COSMIC v105</h2>")
        after = fetch_fingerprint()

    assert before["checksums"] != after["checksums"]
    assert after["source_version"] == "v105"


def test_login_redirect_fails_closed():
    """If the page starts redirecting to the login form there are no version
    tokens, which must surface as probe_failed rather than an empty baseline."""
    with requests_mock.Mocker() as m:
        m.get(COSMIC_RELEASE_NOTES_URL, text="<html><body>Please log in</body></html>")
        with pytest.raises(DriftProbeError, match="named no"):
            fetch_fingerprint()
