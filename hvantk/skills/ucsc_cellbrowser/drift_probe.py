"""UCSC Cell Browser drift probes: one per dataset, each fingerprinting what it builds from.

``probe_version`` 1 was a single zero-argument probe shared by all three datasets. It
hashed the ``ETag`` and ``Content-Length`` of the browser's root index,
``https://cells.ucsc.edu/dataset.json``, and filed the hash under ``checksums`` -- a schema
key to the drift bot (``SCHEMA_KEYS`` in ``.github/scripts/drift_to_pr.py``). The root
index changes whenever any of the browser's ~320 collections changes, so every run that saw
upstream activity opened a ``drift:schema`` PR that said nothing about the collections this
plugin builds from (issue #396): between 2026-08-25 and 2026-09-29 the index gained 56
collections and 38 others changed, while the two parent collections below were untouched.

Moving that hash into ``extras`` alone would have turned uninformative positives into
false negatives, because the probe still could not see a schema (#271). Version 2 instead
gives each dataset the signal it actually has.

``adult-ctx`` and ``dev-ctx`` (``fetch_adult_ctx_fingerprint``, ``fetch_dev_ctx_fingerprint``)
    Per-class summaries derived from one pinned parent collection each:
    ``adult-ctx-meta-atlas`` grouped by ``Class`` and ``dev-ctx-meta-atlas`` grouped by
    ``Type.v2`` (SKILL.md s 4). Two requests per probe: a ranged GET of the parent's
    ``meta.tsv`` first line, and a GET of the parent's own ``dataset.json`` descriptor
    (62 KB and 236 KB respectively).

    * ``headers`` -- the parent's ``meta.tsv`` header row, i.e. the real column list the
      summary is grouped from; and, sorted, the grouping column's labels that the summary
      keeps as rows. Each kept label is one ``obs`` row of the derived dataset and one
      column of its expression TSV, so a class appearing, vanishing or being renamed is a
      layout change to the dataset, not just new content. The labels come from the
      descriptor's ``metaFields`` and pass through the same two filters as
      ``summarize_ucsc_streaming``, which ``hvantk expression summarize-ucsc`` runs. Cells
      whose label ``load_ucsc_metadata``'s ``pd.read_csv`` defaults read as missing (``NA``,
      ``n/a``, an empty field, ...) are dropped by its ``notna()``; pandas itself decides
      which labels those are, so the probe cannot disagree with the build. Groups below
      ``SUMMARIZE_MIN_CELLS_PER_GROUP`` cells (the ``--min-cells`` default, imported from
      ``shared/constants.py``) are dropped too. The cortex derivations are not scripted in
      this repository and record no cut-off of their own, so the default is the one that
      applies. Checked on 2026-10-12: this yields the 19 ``Class`` labels and the 15
      ``Type.v2`` labels of the committed fixtures, ``NA`` (305 cells) excluded.
    * ``checksums`` -- sha256 of the ``meta.tsv`` header line bytes.
    * ``extras`` -- the ``md5`` and ``size`` the descriptor's ``fileVersions`` publishes for
      ``exprMatrix.tsv.gz`` and ``meta.tsv``: the content identity of exactly the two files
      the summary is derived from. Also the cell count of every published label, the
      dropped ones included, so a count that moves, even below the cut-off, shows as
      content drift; only a label crossing the cut-off changes the rows. A re-release that
      keeps the layout is ``routine``.
    * ``informational`` -- the collection-level ``md5`` (the value the root index lists)
      and ``lastBuildTime``. Not compared: the collection ``md5`` also covers coordinates,
      descriptions and the config, none of which reach the summary, and a rebuild bumps the
      build time without changing any data.

``default`` (``fetch_default_fingerprint``)
    Builds whichever collection the user names, so it has no fixed parent and no schema a
    probe could see. It keeps tracking the root index, now as a content signal under
    ``extras`` (the normalised ``ETag`` plus ``Content-Length``), which makes index churn
    ``routine``. ``headers`` and ``checksums`` are empty on purpose: an empty map is the
    honest statement that this dataset has no schema signal, where any hash placed there
    would tier every browser update as a schema change again. ``Last-Modified`` is kept
    under ``informational``; the ``ETag`` already moves with it.

The pinned local catalog (``data/cells_ucsc_datasets.json``) is no longer hashed. It is a
tracked repository file, versioned by git rather than observed upstream, and under
``checksums`` an in-repo edit to it surfaced as an upstream schema change.

Every probe fails closed (``DriftProbeError``) when a response lacks a signal it records,
rather than writing a fingerprint that a later run would compare equal to.
"""

from __future__ import annotations

import csv
import hashlib
import io
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError, normalize_etag
from hvantk.skills.ucsc_cellbrowser.shared.constants import (
    EXPRESSION_MATRIX_FILE_NAME,
    METADATA_FILE_NAME,
    SUMMARIZE_MIN_CELLS_PER_GROUP,
    UCSC_CELL_BROWSER_BASE_URL,
)

PROBE_VERSION = 2
UCSC_DATASET_INDEX_URL = f"{UCSC_CELL_BROWSER_BASE_URL}/dataset.json"

#: Parent collection and grouping column of each derived dataset (SKILL.md s 4, and the
#: ``*_summary`` entries in ``catalog/datasets.json``). The column is the raw ``meta.tsv``
#: header; ``load_ucsc_metadata`` later rewrites ``Type.v2`` to the ``Type_v2`` obs column.
ADULT_CTX_PARENT = "adult-ctx-meta-atlas"
ADULT_CTX_GROUP_COLUMN = "Class"
DEV_CTX_PARENT = "dev-ctx-meta-atlas"
DEV_CTX_GROUP_COLUMN = "Type.v2"

_TIMEOUT_S = 30

#: Bytes requested for the ``meta.tsv`` header row (100 and 129 bytes on 2026-10-12). The
#: read stops at this cap even if the server ignores ``Range``, so a 200 costs one chunk,
#: not the whole 63-98 MB file.
_RANGE_BYTES = 65536

#: ``requests`` asks for gzip by default, and Apache then answers the index HEAD with the
#: compressed length (20,151 rather than 151,816 bytes on 2026-10-12) and an ``ETag`` whose
#: ``-gzip`` suffix sits inside the quotes, where ``normalize_etag`` does not strip it.
#: Asking for the identity encoding records the file's real size and plain ``ETag``, and
#: keeps the ranged GET's byte offsets about the file rather than a compressed stream.
_IDENTITY = {"Accept-Encoding": "identity"}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def fetch_default_fingerprint() -> dict:
    """``default``: the root collection index as a content signal; no schema signal."""
    try:
        head = requests.head(
            UCSC_DATASET_INDEX_URL,
            headers=_IDENTITY,
            timeout=_TIMEOUT_S,
            allow_redirects=True,
        )
        head.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(
            f"HTTP failure on {UCSC_DATASET_INDEX_URL}: {exc}"
        ) from exc

    etag = normalize_etag(head.headers.get("ETag"))
    content_length = head.headers.get("Content-Length")
    if etag is None or not content_length:
        # Content-Length alone cannot stand in for the ETag: a collection update rewrites
        # its fixed-width md5 in the index, which can leave the length unchanged.
        raise DriftProbeError(
            f"{UCSC_DATASET_INDEX_URL} answered without an ETag or Content-Length "
            f"(ETag={head.headers.get('ETag')!r}, Content-Length={content_length!r}); "
            "refusing to record a fingerprint with no content signal."
        )

    return {
        "probe_version": PROBE_VERSION,
        "source_version": None,
        "headers": {},
        "checksums": {},
        "extras": {
            "dataset.json": {"etag": etag, "content_length": content_length},
        },
        "informational": {"last_modified": head.headers.get("Last-Modified")},
        "fetched_at": _now(),
    }


def fetch_adult_ctx_fingerprint() -> dict:
    """``adult-ctx``: the ``Class`` summary of ``adult-ctx-meta-atlas``."""
    return _parent_fingerprint(ADULT_CTX_PARENT, ADULT_CTX_GROUP_COLUMN)


def fetch_dev_ctx_fingerprint() -> dict:
    """``dev-ctx``: the ``Type.v2`` summary of ``dev-ctx-meta-atlas``."""
    return _parent_fingerprint(DEV_CTX_PARENT, DEV_CTX_GROUP_COLUMN)


def _parent_fingerprint(parent: str, group_column: str) -> dict:
    """Fingerprint a pinned parent collection; see the module docstring for the layout."""
    meta_key = f"{parent}/{METADATA_FILE_NAME}"
    header_line = _first_line(f"{UCSC_CELL_BROWSER_BASE_URL}/{meta_key}")
    try:
        columns = header_line.decode("utf-8").split("\t")
    except UnicodeDecodeError as exc:
        raise DriftProbeError(f"{meta_key} header row is not UTF-8: {exc}") from exc

    descriptor = _get_descriptor(parent)
    labels_key = f"{meta_key}:{group_column}"
    rows, counts = _group_rows(descriptor, parent, group_column)
    return {
        "probe_version": PROBE_VERSION,
        "source_version": None,
        "headers": {meta_key: columns, labels_key: rows},
        "checksums": {meta_key: hashlib.sha256(header_line).hexdigest()},
        "extras": {
            f"{parent}/{EXPRESSION_MATRIX_FILE_NAME}": _published_file(
                descriptor, parent, "outMatrix"
            ),
            meta_key: _published_file(descriptor, parent, "outMeta"),
            labels_key: counts,
        },
        "informational": {
            "collection_md5": descriptor.get("md5"),
            "last_build_time": descriptor.get("lastBuildTime"),
        },
        "fetched_at": _now(),
    }


def _first_line(url: str) -> bytes:
    """The first line of ``url``, read with a capped ranged GET (no trailing newline)."""
    buf = bytearray()
    try:
        with requests.get(
            url,
            headers={"Range": f"bytes=0-{_RANGE_BYTES - 1}", **_IDENTITY},
            timeout=_TIMEOUT_S,
            stream=True,
            allow_redirects=True,
        ) as resp:
            # Inside the `with`, so a 4xx/5xx still releases the pooled connection.
            resp.raise_for_status()
            for chunk in resp.iter_content(8192):
                buf += chunk
                if b"\n" in buf or len(buf) >= _RANGE_BYTES:
                    break
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure on {url}: {exc}") from exc

    line, newline, _ = bytes(buf[:_RANGE_BYTES]).partition(b"\n")
    if not newline:
        # Without the terminator the row may be cut short, and a truncated column list
        # would bake into the baseline as if it were the real one.
        raise DriftProbeError(
            f"no complete header row within the first {_RANGE_BYTES} bytes of {url} "
            f"({len(buf)} bytes read)"
        )
    return line.rstrip(b"\r")


def _get_descriptor(parent: str) -> dict:
    """The collection's own ``dataset.json`` (the Cell Browser's per-collection descriptor)."""
    url = f"{UCSC_CELL_BROWSER_BASE_URL}/{parent}/dataset.json"
    try:
        resp = requests.get(url, timeout=_TIMEOUT_S, allow_redirects=True)
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure on {url}: {exc}") from exc
    try:
        descriptor = resp.json()
    except ValueError as exc:
        raise DriftProbeError(f"{url} did not return JSON: {exc}") from exc
    if not isinstance(descriptor, dict):
        raise DriftProbeError(
            f"{url} returned a JSON {type(descriptor).__name__}, not an object"
        )
    return descriptor


def _group_rows(
    descriptor: dict, parent: str, column: str
) -> tuple[list[str], dict[str, int]]:
    """The summary rows ``column`` yields, and the cell count of every published label.

    Rows: the labels ``summarize_ucsc_streaming`` keeps -- not read as missing, and at
    least ``SUMMARIZE_MIN_CELLS_PER_GROUP`` cells -- sorted, because ``valCounts`` is
    ordered by cell count, which a re-annotation can reshuffle without adding or
    removing a class. The labels are matched on the column's raw ``meta.tsv`` header,
    which the descriptor calls ``label``.

    A column missing from ``metaFields`` yields ``([], {})`` rather than an error: that is
    the upstream schema change this entry exists to report, and the ``meta.tsv`` header
    row recorded beside it moves too.
    """
    fields = descriptor.get("metaFields")
    if not isinstance(fields, list):
        raise DriftProbeError(
            f"{parent}/dataset.json has no metaFields list; the Cell Browser descriptor "
            "format has changed"
        )
    field = next(
        (f for f in fields if isinstance(f, dict) and f.get("label") == column), None
    )
    if field is None:
        return [], {}
    try:
        counts = {
            str(label): int(count) for label, count in field.get("valCounts") or []
        }
    except (TypeError, ValueError) as exc:
        raise DriftProbeError(
            f"{parent}/dataset.json: unreadable valCounts for {column!r}: {exc}"
        ) from exc
    missing = _read_as_missing(list(counts))
    rows = sorted(
        label
        for label, count in counts.items()
        if label not in missing and count >= SUMMARIZE_MIN_CELLS_PER_GROUP
    )
    return rows, dict(sorted(counts.items()))


def _read_as_missing(labels: list[str]) -> set[str]:
    """The labels pandas' ``read_csv`` defaults turn into NaN.

    That is how ``load_ucsc_metadata`` reads ``meta.tsv`` before ``summarize_ucsc_streaming``
    drops NaN-labelled cells with ``notna()``. Asking pandas itself, rather than copying
    its NA list, keeps the probe and the build from disagreeing. Each label gets its own
    indexed line, so an empty label stays a row. ``QUOTE_NONE`` pins one row per label;
    it departs from the real read only for a label containing a double quote, and no
    ``Class`` or ``Type.v2`` label contains one. pandas is imported here, not at module
    level, so loading the plugin registry does not pay for it.
    """
    import pandas as pd

    if not labels:
        return set()
    text = "".join(f"{i}\t{label}\n" for i, label in enumerate(labels))
    try:
        parsed = pd.read_csv(
            io.StringIO(text),
            sep="\t",
            header=None,
            index_col=0,
            dtype=str,
            quoting=csv.QUOTE_NONE,
        ).iloc[:, 0]
    except ValueError as exc:  # pandas' ParserError and EmptyDataError included
        raise DriftProbeError(f"could not parse the published labels: {exc}") from exc
    if len(parsed) != len(labels):
        raise DriftProbeError(
            f"parsed {len(parsed)} labels from {len(labels)}; refusing to guess which "
            "ones pandas reads as missing"
        )
    return {
        label for label, value in zip(labels, parsed, strict=True) if pd.isna(value)
    }


def _published_file(descriptor: dict, parent: str, kind: str) -> dict:
    """``{"md5", "size"}`` the descriptor's ``fileVersions`` publishes for ``kind``."""
    versions = descriptor.get("fileVersions")
    entry = versions.get(kind) if isinstance(versions, dict) else None
    md5 = entry.get("md5") if isinstance(entry, dict) else None
    size = entry.get("size") if isinstance(entry, dict) else None
    if not md5 or size is None:
        raise DriftProbeError(
            f"{parent}/dataset.json publishes no md5/size for {kind!r}; refusing to record "
            "a fingerprint with no content signal for it."
        )
    return {"md5": md5, "size": size}
