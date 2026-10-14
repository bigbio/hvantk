"""dbnsfp drift probe: the academic release's column list, found via the releases page.

probe_version 1 and 2 watched the legacy Google Sites landing page
(``https://sites.google.com/site/jpopgen/dbNSFP``), which is frozen at v4.9 and tells
readers to use dbNSFP.org, so they reported clean through v5.1-v5.4 (issue #321).
probe_version 3 scraped the release titles from ``https://www.dbnsfp.org/releases/`` and,
being unable to see any schema, kept that list under ``headers``: per `_conventions`
section 12 a pure content signal belongs in ``extras`` only once the probe can see the
schema (issue #271), so every release opened its own ``drift:schema`` PR.

**What it fingerprints now (probe_version 4).** Since v5.4 the releases page links each
branch's column list as a plain-text file beside its README -- for the academic branch,
``https://dist.genos.us/release/dbNSFP5.4a_variant.columns.txt``. Verified 2026-10-13: a
public 82 KB file, served without the registration the data itself needs (SKILL.md
section 2), a tab-separated ``index``/``name``/``description`` table with one row per
column of the variant TSV (508 for v5.4a). The probe follows that link for the latest
academic release and records:

* ``headers`` -- the column names, in file order (``#chr`` first, with the hash the TSV
  header carries). A real schema signal: a release that adds columns, as v5.4 added
  the GPN-MSA scores, now opens a ``drift:schema`` PR because its columns moved.
* ``checksums`` -- sha256 over those names tab-joined, i.e. the header row the release
  declares for its TSV, the shape the header-hashing probes (clingen, gencc, hgnc) use.
  The ``description`` column is outside both: it is prose that changes between releases
  with no column changing (of the 503 columns v5.3.1a, whose README still carried the
  descriptions, shares with v5.4a, 2 were reworded), and a schema tier that moves on a
  reworded sentence is the constant tier of issue #333.
* ``source_version`` (the latest academic release) and ``extras.releases`` (every release
  title on the page) -- compared, but not schema keys
  (``SCHEMA_KEYS`` in ``hvantk/core/plugin/drift_ledger.py``). A new release with an
  unchanged column list is therefore ``drift:routine``, and a commercial-only patch still
  moves ``extras``.
* ``informational.variant_columns_url`` -- the file that was read, for the reviewer.

The key under ``headers``/``checksums`` is version-free (``dbNSFP_variant``). Keyed by
the versioned file name instead, every release would rename it and read as a schema
change regardless of its columns.

**Academic or commercial: the README links, not the version title.** Each release lists
one README per branch that shipped: ``dbNSFP5.4a.readme.txt`` and
``dbNSFP5.4c.readme.txt`` for v5.4, but only ``dbNSFP5.1.1c.readme.txt`` for the
commercial-only v5.1.1c patch (verified 2026-10-13). An academic README is direct
evidence that an academic build exists. probe_version 3 inferred it from the title
instead, treating any version printed without a ``c`` suffix as academic, so a
commercial-only patch printed bare would have been taken for an academic release.

**The raw page body is still never hashed.** Its footer carries a copyright year range
that rolls forward every January with no release involved.

**Identifies as hvantk.** dbnsfp.org and dist.genos.us both sit behind Cloudflare, where
a bot challenge served to the scheduled runner would fail the probe closed on every run.
Both requests send ``USER_AGENT``, as hvantk's downloaders do, instead of requests'
default ``python-requests/<version>``.

What it cannot detect: the columns file is the publisher's declaration of the TSV
header, not the registration-gated archive (~50 GB) itself, so a release whose files
diverge from their own column list, or an in-place change to file contents, goes unseen.

Fails closed (``DriftProbeError``) rather than record a fingerprint a later run would
compare equal to: no release-title entries, a title not shaped
``dbNSFP v<version> (<date>)``, no academic README on the page, no variant columns link
for the latest academic release, or a columns file that is not the
``index``/``name``/``description`` table it has been.
"""

from __future__ import annotations

import hashlib
import posixpath
import re
from datetime import datetime, timezone
from urllib.parse import unquote, urljoin, urlsplit

import requests

from hvantk.core.plugin.api import DriftProbeError
from hvantk.core.utils.http import USER_AGENT, request_with_retry

PROBE_VERSION = 4
DBNSFP_RELEASES_URL = "https://www.dbnsfp.org/releases/"

# Every release -- current or past -- renders its title in one of these spans (verified
# live 2026-09-30 and 2026-10-13).
_RELEASE_TITLE_REGEX = re.compile(r'<span class="release-title">([^<]+)</span>')

# "dbNSFP v5.4 (August 1, 2026)", "dbNSFP v5.1.1c (April 24, 2025)". Only validated: a
# title that does not fit means the page's release format changed.
_VERSION_DATE_REGEX = re.compile(r"^dbNSFP\s+v[0-9][\w.]*\s*\([^)]+\)$", re.IGNORECASE)

_HREF_REGEX = re.compile(r'href="([^"]+)"')

# One README per branch that shipped: dbNSFP5.4a.readme.txt, dbNSFP5.1.1c.readme.txt.
_README_REGEX = re.compile(r"^dbNSFP(\d+(?:\.\d+)*)([ac])\.readme\.txt$", re.IGNORECASE)

# The academic variant table's column list, linked beside the README since v5.4.
_COLUMNS_FILE = "dbNSFP{version}a_variant.columns.txt"
_COLUMNS_HEADER = ("index", "name", "description")
_COLUMNS_KEY = "dbNSFP_variant"

_HEADERS = {"User-Agent": USER_AGENT}
_TIMEOUT_S = (5.0, 15.0)


def _get_text(url: str) -> str:
    """GET ``url`` as hvantk and return its body as UTF-8 text.

    Decoded explicitly rather than via ``resp.text``: neither host sends a charset in
    Content-Type, so requests would fall back to ISO-8859-1 for text and ignore the
    page's own ``<meta charset="utf-8">``.
    """
    try:
        resp = request_with_retry(
            "GET",
            url,
            timeout=_TIMEOUT_S,
            allow_redirects=True,
            headers=_HEADERS,
            retry_on_empty_body=True,
        )
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc
    return resp.content.decode("utf-8", errors="replace")


def _file_name(href: str) -> str:
    return posixpath.basename(unquote(urlsplit(href).path))


def _numeric_parts(version: str) -> tuple[int, ...]:
    return tuple(int(p) for p in version.split("."))


def _academic_versions(hrefs: list[str]) -> set[str]:
    """Versions whose release links an academic (``a``) README."""
    versions = set()
    for href in hrefs:
        m = _README_REGEX.match(_file_name(href))
        if m and m.group(2).lower() == "a":
            versions.add(m.group(1))
    return versions


def _columns_href(hrefs: list[str], version: str) -> str | None:
    """The link to ``version``'s academic variant column list, if the page has one."""
    wanted = _COLUMNS_FILE.format(version=version).lower()
    return next((h for h in hrefs if _file_name(h).lower() == wanted), None)


def _parse_columns(text: str, url: str) -> list[str]:
    """Column names from an ``index``/``name``/``description`` table, in index order.

    Fails closed on anything else -- a different header, a row whose index does not
    follow its predecessor (a truncated or merged row), a blank name, or no rows at all
    -- since a partial list would be committed as the baseline.
    """
    lines = [line for line in text.splitlines() if line.strip()]
    header = tuple(f.strip().lower() for f in lines[0].split("\t")) if lines else ()
    if header != _COLUMNS_HEADER:
        raise DriftProbeError(
            f"dbNSFP column list at {url} does not start with an "
            f"'index<TAB>name<TAB>description' header (got {header!r}); its format "
            "has probably changed."
        )
    names: list[str] = []
    for expected, line in enumerate(lines[1:], start=1):
        fields = line.split("\t")
        name = fields[1].strip() if len(fields) > 1 else ""
        if fields[0].strip() != str(expected) or not name:
            raise DriftProbeError(
                f"dbNSFP column list at {url}: row {expected} reads {line[:80]!r}, not "
                f"'{expected}<TAB><name><TAB>...'; its format has probably changed."
            )
        names.append(name)
    if not names:
        raise DriftProbeError(f"dbNSFP column list at {url} lists no columns.")
    return names


def fetch_fingerprint() -> dict:
    """Fingerprint the latest academic dbNSFP release's declared variant columns."""
    body = _get_text(DBNSFP_RELEASES_URL)

    titles = [t.strip() for t in _RELEASE_TITLE_REGEX.findall(body)]
    # No release-title spans means the page moved or its markup changed, not that
    # dbNSFP shipped zero releases.
    if not titles:
        raise DriftProbeError(
            "dbNSFP releases page contained no release-title entries; "
            "the page layout has probably changed."
        )
    for title in titles:
        if not _VERSION_DATE_REGEX.match(title):
            raise DriftProbeError(
                f"dbNSFP release entry {title!r} did not match the expected "
                "'dbNSFP v<version> (<date>)' shape; the page layout has probably "
                "changed."
            )

    hrefs = _HREF_REGEX.findall(body)
    academic = _academic_versions(hrefs)
    if not academic:
        raise DriftProbeError(
            f"dbNSFP advertised {len(titles)} release(s) but none links an academic "
            "README (dbNSFP<version>a.readme.txt); the release scheme or the page "
            "layout has probably changed."
        )
    latest = max(academic, key=_numeric_parts)  # componentwise: 5.10 beats 5.4

    href = _columns_href(hrefs, latest)
    if href is None:
        raise DriftProbeError(
            f"the latest academic dbNSFP release, v{latest}a, links no variant columns "
            f"file ({_COLUMNS_FILE.format(version=latest)}), so its column list "
            "cannot be read."
        )
    columns_url = urljoin(DBNSFP_RELEASES_URL, href)
    names = _parse_columns(_get_text(columns_url), columns_url)

    return {
        "probe_version": PROBE_VERSION,
        "source_version": latest,
        "headers": {_COLUMNS_KEY: names},
        "checksums": {
            _COLUMNS_KEY: hashlib.sha256("\t".join(names).encode("utf-8")).hexdigest()
        },
        "extras": {"releases": sorted(set(titles))},
        "informational": {"variant_columns_url": columns_url},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
