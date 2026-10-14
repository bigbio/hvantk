"""dbnsfp drift probe: the academic release's column list, found via the releases page."""

from __future__ import annotations

import hashlib
import re

import pytest
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.core.plugin.drift_ledger import classify_risk
from hvantk.core.plugin.drift_runner import _compare_fingerprints
from hvantk.core.utils.http import USER_AGENT
from hvantk.skills.dbnsfp.drift_probe import DBNSFP_RELEASES_URL, fetch_fingerprint

_DIST = "https://dist.genos.us/release"
_ANY_COLUMNS_FILE = re.compile(r".*_variant\.columns\.txt$")

# Shaped like the live https://www.dbnsfp.org/releases/ page (verified 2026-10-13):
# (title version, date, branches with a README link). A shared release links an
# academic and a commercial README; v5.1.1c shipped commercial-only and links only its
# own README, with the branch suffix glued onto the title's version.
_RELEASES = [
    ("5.4", "August 1, 2026", "ac"),
    ("5.3.1", "January 1, 2026", "ac"),
    ("5.3", "October 6, 2025", "ac"),
    ("5.2", "July 2, 2025", "ac"),
    ("5.1.1c", "April 24, 2025", "c"),
    ("5.1", "March 21, 2025", "ac"),
    ("5.0", "January 1, 2025", "ac"),
]

_COLUMNS = ["#chr", "pos(1-based)", "ref", "alt", "REVEL_score"]


def _columns_tsv(names):
    """A column list shaped like the live dbNSFP5.4a_variant.columns.txt."""
    rows = ["index\tname\tdescription"]
    rows += [f"{i}\t{name}\tdescription of {name}" for i, name in enumerate(names, 1)]
    return "\n".join(rows) + "\n"


def _page(releases, columns=True):
    """Minimal page shaped like the real one: one ``release-title`` span per entry, then
    its README links, plus the academic column-list link when ``columns`` is set (the
    live page links it on the current release only; the probe reads only the latest
    academic one, so linking it everywhere changes nothing it sees).
    """
    items = []
    for version, date, branches in releases:
        number = version.rstrip("acAC")
        links = [
            f'<a href="{_DIST}/dbNSFP{number}{b}.readme.txt">README v{number}{b}</a>'
            for b in branches
        ]
        if columns and "a" in branches:
            links.append(
                f'<a href="{_DIST}/dbNSFP{number}a_variant.columns.txt">'
                f"Variant Columns v{number}a</a>"
            )
        items.append(
            '<div class="release-item"><p class="release-heading">'
            f'<span class="release-title">dbNSFP v{version} ({date})</span></p>'
            f'<p class="release-files">{"".join(links)}</p></div>'
        )
    return f"<html><body><h1>Current Release</h1>{''.join(items)}</body></html>"


def _mock(m, page, columns=_COLUMNS):
    m.get(DBNSFP_RELEASES_URL, text=page)
    m.get(_ANY_COLUMNS_FILE, text=_columns_tsv(columns))


def test_reports_latest_academic_release():
    """The newest release linking an academic README wins `source_version`, and its own
    column list -- not the page -- is the schema signal."""
    with requests_mock.Mocker() as m:
        _mock(m, _page(_RELEASES))
        fp = fetch_fingerprint()

    assert fp["source_version"] == "5.4"
    assert fp["informational"]["variant_columns_url"] == (
        f"{_DIST}/dbNSFP5.4a_variant.columns.txt"
    )
    assert fp["headers"] == {"dbNSFP_variant": _COLUMNS}
    assert fp["checksums"] == {
        "dbNSFP_variant": hashlib.sha256("\t".join(_COLUMNS).encode()).hexdigest()
    }
    # Both requests identify as hvantk, not as python-requests (#408).
    assert [r.headers["User-Agent"] for r in m.request_history] == [USER_AGENT] * 2


def test_commercial_only_patch_does_not_win_latest_academic():
    """A commercial-only patch (the live v5.1.1c pattern) must not be reported as
    the latest academic release even when it is the newest entry on the page --
    there is no academic file behind it."""
    with requests_mock.Mocker() as m:
        _mock(
            m,
            _page(
                [("5.4", "August 1, 2026", "ac"), ("5.4.1c", "September 1, 2026", "c")]
            ),
        )
        fp = fetch_fingerprint()

    assert fp["source_version"] == "5.4"
    assert any("5.4.1c" in title for title in fp["extras"]["releases"])


def test_a_bare_title_with_only_a_commercial_readme_is_not_academic():
    """The README link decides, not the title: a commercial-only patch printed without
    its suffix (here with an uppercase `C` README) is still not academic. Reading any
    bare title as academic, as probe_version 3 did, would report it as the latest."""
    with requests_mock.Mocker() as m:
        _mock(
            m,
            _page(
                [("5.4", "August 1, 2026", "ac"), ("5.4.1", "September 1, 2026", "C")]
            ),
        )
        fp = fetch_fingerprint()

    assert fp["source_version"] == "5.4"


def test_versions_compare_componentwise_not_as_text():
    """5.10 must beat 5.4. The committed baseline sits in the 5.x series, so the
    first double-digit minor release is the one a text sort gets wrong."""
    with requests_mock.Mocker() as m:
        _mock(m, _page([("5.4", "d1", "ac"), ("5.10", "d2", "ac")]))
        fp = fetch_fingerprint()

    assert fp["source_version"] == "5.10"


def test_projection_ignores_markup_and_link_order():
    """Reordering the release entries, or adding unrelated markup, must not move
    the fingerprint -- only the release SET and the column list matter."""
    with requests_mock.Mocker() as m:
        _mock(m, _page(_RELEASES))
        first = fetch_fingerprint()

    reordered = _page(list(reversed(_RELEASES))).replace(
        "<body>", "<body><p>unrelated editorial edit</p>"
    )
    with requests_mock.Mocker() as m:
        _mock(m, reordered)
        second = fetch_fingerprint()

    for key in ("source_version", "headers", "checksums", "extras"):
        assert first[key] == second[key], key


def test_new_release_moves_the_signal():
    """A new release moves `source_version` and `extras`; its tier then follows the
    columns. The same column list is routine (batched), a changed one is a schema PR --
    probe_version 3 tiered every release as schema."""
    with requests_mock.Mocker() as m:
        _mock(m, _page(_RELEASES))
        before = fetch_fingerprint()

    newer = [("5.5", "d", "ac"), *_RELEASES]
    with requests_mock.Mocker() as m:
        _mock(m, _page(newer))
        same_columns = fetch_fingerprint()
    with requests_mock.Mocker() as m:
        _mock(m, _page(newer), columns=[*_COLUMNS, "NewTool_score"])
        new_column = fetch_fingerprint()

    assert same_columns["source_version"] == "5.5"
    assert before["extras"] != same_columns["extras"]
    assert classify_risk(_compare_fingerprints(before, same_columns)) == "routine"
    assert classify_risk(_compare_fingerprints(before, new_column)) == "schema"


def test_empty_release_list_fails_closed():
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_RELEASES_URL, text="<html><body>no releases here</body></html>")
        with pytest.raises(DriftProbeError, match="no release-title"):
            fetch_fingerprint()


@pytest.mark.parametrize(
    "page,columns_tsv,match",
    [
        pytest.param(
            '<html><body><span class="release-title">'
            "dbNSFP version five point four"
            "</span></body></html>",
            None,
            "did not match",
            id="release-title",
        ),
        pytest.param(
            _page(_RELEASES, columns=False),
            None,
            "links no variant columns file",
            id="no-columns-link",
        ),
        pytest.param(
            _page(_RELEASES),
            "name\tdescription\n#chr\tchromosome\n",
            "does not start with",
            id="columns-header",
        ),
        pytest.param(
            _page(_RELEASES),
            "index\tname\tdescription\n1\t#chr\tx\n3\tref\tx\n",
            "row 2",
            id="columns-row",
        ),
    ],
)
def test_unrecognised_release_format_fails_closed(page, columns_tsv, match):
    """A layout change on the page or in the column list raises instead of recording
    a partial fingerprint that every later run would compare equal to."""
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_RELEASES_URL, text=page)
        if columns_tsv is not None:
            m.get(_ANY_COLUMNS_FILE, text=columns_tsv)
        with pytest.raises(DriftProbeError, match=match):
            fetch_fingerprint()


def test_no_academic_release_fails_closed():
    """Silently recording source_version: null would let the bot commit that null,
    after which the probe reports clean having stopped identifying a release."""
    with requests_mock.Mocker() as m:
        m.get(DBNSFP_RELEASES_URL, text=_page([("5.4.1c", "September 1, 2026", "c")]))
        with pytest.raises(DriftProbeError, match="academic"):
            fetch_fingerprint()
