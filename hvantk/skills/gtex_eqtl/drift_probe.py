"""GTEx eQTL drift probe: the pinned v11 release, read from the bucket GTEx serves it from.

``probe_version`` 1 fingerprinted the GTEx Portal download *page*: it hashed the page's
``gtex-version`` and ``gtex-build-date`` meta tags plus its ``ETag`` and filed the hash
under ``checksums``, a schema key to the drift bot (``SCHEMA_KEYS`` in
``.github/scripts/drift_to_pr.py``). Those tags identify a build of the portal web
application, not a data release. On 2026-10-12 the page read ``cl363``, built 2026-09-22,
while every v11 file this plugin reads was uploaded on 2026-01-15 and has not changed since;
the redeploy from ``cl362`` was nonetheless tiered as a schema change (bot PR #403).

The GTEx API's dataset list (``/api/v2/metadata/dataset``) is no better a release signal
here: it lists ``gtex_v7``, ``gtex_v8``, ``kids_first_harmonization`` and ``gtex_v10`` --
the releases the portal serves through its API -- and not v11, which has been downloadable
since January 2026.

Version 2 reads the source the portal's download links point at: the public ``adult-gtex``
Google Cloud Storage bucket, which allows anonymous listing through the JSON API. Three
small GETs (two listings and a 7 KB README):

* ``headers``

  - ``gs://adult-gtex/bulk-qtl/`` -- the bulk-QTL release folders (``v10``, ``v11``, ...).
    A new GTEx release appears here. Adopting one is a builder change (SKILL.md s 8: a new
    source constant, possibly a new import branch -- v8 to v11 went from TSV to parquet),
    which is what the schema tier asks a reviewer to check.
  - ``gs://adult-gtex/bulk-qtl/v11/single-tissue-cis-qtl/`` -- the object names in the
    pinned release's cis-QTL folder, i.e. its file layout. The builder's
    ``*.v11.eQTLs.signif_pairs.parquet`` inputs ship inside ``GTEx_Analysis_v11_eQTL.tar``
    (the Liver member is exactly the 27,912,129 bytes the catalog records).
  - ``README_eQTL_v11.txt:signif_pairs`` -- the columns the release README documents for
    ``*.signif_pairs.parquet``, in documented order. This is GTEx's *declared* schema, not
    the parquet's own: the real files already differ from it in two known ways (SKILL.md
    s 4). Reading a real parquet footer would mean walking the 3.3 GB archive's tar
    headers first, one ranged request per member -- 100 requests and ~40 s, measured on
    2026-10-12 -- which no scheduled probe should do.

* ``checksums`` -- sha256 of that README section's text, so a changed column description
  (the ``af`` direction the builder depends on, say) or file extension moves it even when
  the column names do not. The rest of the README (eGenes, allpairs, fine-mapping) is out
  of scope for this builder and deliberately not hashed.
* ``extras`` -- ``md5Hash`` and ``size`` of ``GTEx_Analysis_v11_eQTL.tar``, the content
  identity of the archive holding the builder's inputs. A re-upload of v11 data that
  keeps the layout is ``routine``. Recorded as ``null`` when the archive is missing from
  the folder; the folder listing under ``headers`` reports that layout change.
* ``source_version`` -- ``v11``, the release this probe is pinned to. A constant: moving
  to a newer release is a change to this probe, not an observation of it.

Portal build tags, build dates and page ``ETag``s are not recorded at all, not even under
``informational``: they describe deployments of the web application, and fetching the
page for them would add a request, and an HTML scrape that can fail the probe, for no
drift signal.

Fails closed (``DriftProbeError``) when a listing is empty, paginated or unreadable, or
when the README no longer documents the files this builder reads.
"""

from __future__ import annotations

import hashlib
import re
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError

PROBE_VERSION = 2

GTEX_RELEASE = "v11"
GTEX_BUCKET = "adult-gtex"
#: Google Cloud Storage JSON API object listing for the bucket (anonymous access).
GTEX_LIST_URL = f"https://storage.googleapis.com/storage/v1/b/{GTEX_BUCKET}/o"
_RELEASES_PREFIX = "bulk-qtl/"
_RELEASE_PREFIX = f"{_RELEASES_PREFIX}{GTEX_RELEASE}/single-tissue-cis-qtl/"
_ARCHIVE = f"GTEx_Analysis_{GTEX_RELEASE}_eQTL.tar"
_README = f"README_eQTL_{GTEX_RELEASE}.txt"
GTEX_README_URL = (
    f"https://storage.googleapis.com/{GTEX_BUCKET}/{_RELEASE_PREFIX}{_README}"
)

_RELEASES_KEY = f"gs://{GTEX_BUCKET}/{_RELEASES_PREFIX}"
_RELEASE_KEY = f"gs://{GTEX_BUCKET}/{_RELEASE_PREFIX}"
_README_KEY = f"{_README}:signif_pairs"

_TIMEOUT_S = 30

#: The README's section rule (``#----...``). Sections alternate title / body between rules.
_SECTION_RULE = re.compile(r"^#-{10,}[ \t]*$", re.M)
#: The body that documents the builder's input: its "File extension(s):" line names it.
_SIGNIF_PAIRS_EXTENSION = re.compile(r"^\s*File extensions?:.*signif_pairs", re.M)
#: A documented column: an indented ``name:`` line. "File extensions:" and
#: "Column headers:" never match, because their labels contain a space.
_COLUMN_LINE = re.compile(r"^\s+([A-Za-z0-9_.]+):\s", re.M)


def fetch_fingerprint() -> dict:
    """Fingerprint the pinned GTEx v11 cis-eQTL release; see the module docstring."""
    prefixes = _list(_RELEASES_PREFIX, fields="prefixes,nextPageToken").get("prefixes")
    releases = sorted(
        p[len(_RELEASES_PREFIX) :].rstrip("/")
        for p in prefixes or []
        if isinstance(p, str) and p.startswith(_RELEASES_PREFIX)
    )
    if not releases:
        raise DriftProbeError(
            f"{_RELEASES_KEY} lists no release folders; the bucket layout has changed"
        )

    items = _list(_RELEASE_PREFIX, fields="items(name,size,md5Hash),nextPageToken").get(
        "items"
    )
    objects = {
        item["name"][len(_RELEASE_PREFIX) :]: item
        for item in items or []
        if isinstance(item, dict)
        and isinstance(item.get("name"), str)
        and item["name"].startswith(_RELEASE_PREFIX)
    }
    if not objects:
        raise DriftProbeError(
            f"{_RELEASE_KEY} lists no objects; the pinned {GTEX_RELEASE} release has "
            "moved or been withdrawn"
        )
    archive = objects.get(_ARCHIVE)

    section = _signif_pairs_section(_get_text(GTEX_README_URL))
    columns = _COLUMN_LINE.findall(section)
    if not columns:
        raise DriftProbeError(
            f"{_README}: the signif_pairs section documents no columns; its format has "
            "changed"
        )

    return {
        "probe_version": PROBE_VERSION,
        "source_version": GTEX_RELEASE,
        "headers": {
            _RELEASES_KEY: releases,
            _RELEASE_KEY: sorted(objects),
            _README_KEY: columns,
        },
        "checksums": {
            _README_KEY: hashlib.sha256(section.encode("utf-8")).hexdigest(),
        },
        "extras": {
            _ARCHIVE: (
                {"md5Hash": archive.get("md5Hash"), "size": archive.get("size")}
                if archive is not None
                else None
            ),
        },
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }


def _list(prefix: str, *, fields: str) -> dict:
    """One page of the bucket's object listing under ``prefix`` (direct children only)."""
    try:
        resp = requests.get(
            GTEX_LIST_URL,
            params={"prefix": prefix, "delimiter": "/", "fields": fields},
            timeout=_TIMEOUT_S,
        )
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(
            f"HTTP failure listing gs://{GTEX_BUCKET}/{prefix}: {exc}"
        ) from exc
    try:
        payload = resp.json()
    except ValueError as exc:
        raise DriftProbeError(
            f"listing of gs://{GTEX_BUCKET}/{prefix} was not JSON: {exc}"
        ) from exc
    if not isinstance(payload, dict):
        raise DriftProbeError(
            f"listing of gs://{GTEX_BUCKET}/{prefix} returned a JSON "
            f"{type(payload).__name__}, not an object"
        )
    if payload.get("nextPageToken"):
        # One page holds 1,000 entries; these folders hold about a dozen. A second page
        # would mean the layout changed beyond recognition, and fingerprinting the first
        # page alone would silently drop the rest.
        raise DriftProbeError(
            f"listing of gs://{GTEX_BUCKET}/{prefix} spans more than one page"
        )
    return payload


def _get_text(url: str) -> str:
    try:
        resp = requests.get(url, timeout=_TIMEOUT_S)
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure on {url}: {exc}") from exc
    # Decoded from the bytes, not `resp.text`: the object is served as text/plain with no
    # charset, for which requests would guess ISO-8859-1.
    return resp.content.decode("utf-8", errors="replace")


def _signif_pairs_section(readme: str) -> str:
    """Body of the README section that documents ``*.signif_pairs.*`` files."""
    for block in _SECTION_RULE.split(readme):
        if _SIGNIF_PAIRS_EXTENSION.search(block):
            return block.strip()
    raise DriftProbeError(
        f"{_README} no longer documents *.signif_pairs.* files, the builder's input; "
        "the release format has changed"
    )
