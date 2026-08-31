"""dbnsfp drift probe should fingerprint the advertised release list, not the page.

Runs OFFLINE via requests_mock. The probe replaced a stub sentinel (issue #177)
once the landing page was confirmed to expose a stable release list, even though
every archive it links is dead.
"""

from __future__ import annotations

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.dbnsfp.drift_probe import DBNSFP_LANDING_URL, fetch_fingerprint

_PAGE = """
<html><body>
  <a href="https://dbnsfp.s3.amazonaws.com/dbNSFP4.8a.zip">4.8a</a>
  <a href="https://dbnsfp.s3.amazonaws.com/dbNSFP4.8c.zip">4.8c</a>
  <a href="https://dbnsfp.s3.amazonaws.com/dbNSFP4.9a.zip">4.9a</a>
  <a href="https://dbnsfp.s3.amazonaws.com/dbNSFP4.9c.zip">4.9c</a>
  <a href="https://dbnsfp.s3.amazonaws.com/dbNSFP2.0b1_variant.zip">2.0b1</a>
</body></html>
"""


def test_fetch_fingerprint_reports_latest_academic_release():
    """`a` is the academic build hvantk uses; `c` is commercial and must not win."""
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text=_PAGE)
        fp = fetch_fingerprint()

    assert fp["source_version"] == "4.9a"
    assert fp["extras"]["releases_found"] == [
        "2.0b1_variant",
        "4.8a",
        "4.8c",
        "4.9a",
        "4.9c",
    ]


def test_checksum_covers_the_release_list_not_the_page_body():
    """Google Sites re-renders per request -- two consecutive live fetches returned
    352,830 and 352,716 bytes -- so hashing the body would flag drift on every run.
    The digest must depend only on the extracted release list."""
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text=_PAGE)
        first = fetch_fingerprint()

    # Same releases, different surrounding markup and ordering.
    shuffled = _PAGE.replace("<html><body>", "<html><body><p>unrelated edit</p>")
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text=shuffled)
        second = fetch_fingerprint()

    assert first["checksums"] == second["checksums"]


def test_new_release_moves_the_checksum():
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text=_PAGE)
        before = fetch_fingerprint()

    with requests_mock.Mocker() as m:
        m.get(
            DBNSFP_LANDING_URL,
            text=_PAGE.replace(
                "</body>",
                '<a href="https://dbnsfp.s3.amazonaws.com/dbNSFP5.0a.zip">5.0a</a></body>',
            ),
        )
        after = fetch_fingerprint()

    assert before["checksums"] != after["checksums"]
    assert after["source_version"] == "5.0a"


def test_empty_release_list_fails_closed():
    """No matches means the page layout changed, not that dbNSFP has no releases."""
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_LANDING_URL, text="<html><body>no releases here</body></html>")
        with pytest.raises(DriftProbeError, match="advertised no"):
            fetch_fingerprint()
