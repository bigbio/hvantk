"""UCSC Cell Browser drift probes: one per dataset (issue #396).

Runs OFFLINE - every upstream request is stubbed via requests_mock.

**Expected URLs are pinned as literals here**, not imported from the module under test
(issue #341, item 2). Importing a URL constant and feeding it to the mock made the URL
assertion circular: a probe pointing at ``https://cells.ucsc.edu/datasets.json`` (plural -
a real mistake made while re-authoring this probe) registered its own wrong URL with the
mocker and passed green. Registering against the literals below makes a probe aimed
anywhere else raise ``NoMockAddress`` instead.
"""

from __future__ import annotations

import hashlib

import pytest
import requests
import requests_mock

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.ucsc_cellbrowser.drift_probe import (
    UCSC_DATASET_INDEX_URL,
    fetch_adult_ctx_fingerprint,
    fetch_default_fingerprint,
    fetch_dev_ctx_fingerprint,
)
from hvantk.skills.ucsc_cellbrowser.shared.constants import (
    SUMMARIZE_MIN_CELLS_PER_GROUP as CUT_OFF,
)

# Note the SINGULAR "dataset.json"; "datasets.json" does not exist upstream.
EXPECTED_INDEX_URL = "https://cells.ucsc.edu/dataset.json"

_INDEX_HEADERS = {
    "Last-Modified": "Wed, 30 Sep 2026 23:39:35 GMT",
    "ETag": '"25108-65cbbcfd8767d"',
    "Content-Length": "151816",
}


def test_probe_targets_the_documented_index_url():
    """The probe must point at the real upstream index, not a plausible-looking one."""
    assert UCSC_DATASET_INDEX_URL == EXPECTED_INDEX_URL


def test_default_tracks_the_index_as_content_with_no_schema_signal():
    """`default` has no fixed collection: index churn must tier `routine`, never `schema`."""
    with requests_mock.Mocker() as m:
        m.head(EXPECTED_INDEX_URL, headers=_INDEX_HEADERS)
        fp = fetch_default_fingerprint()
        # Identity encoding, or Apache reports the gzip length and a "-gzip" ETag.
        assert m.last_request.headers["Accept-Encoding"] == "identity"

    assert fp["probe_version"] == 2
    assert fp["source_version"] is None
    # Empty schema keys are the point of #396: nothing here may move `headers`/`checksums`.
    assert fp["headers"] == {} and fp["checksums"] == {}
    assert fp["extras"] == {
        "dataset.json": {"etag": "25108-65cbbcfd8767d", "content_length": "151816"}
    }
    assert fp["informational"] == {"last_modified": "Wed, 30 Sep 2026 23:39:35 GMT"}
    assert "fetched_at" in fp


def test_default_fails_closed_without_an_etag():
    headers = {k: v for k, v in _INDEX_HEADERS.items() if k != "ETag"}
    with requests_mock.Mocker() as m:
        m.head(EXPECTED_INDEX_URL, headers=headers)
        with pytest.raises(DriftProbeError, match="no content signal"):
            fetch_default_fingerprint()


def test_default_raises_on_request_failure():
    with requests_mock.Mocker() as m:
        m.head(EXPECTED_INDEX_URL, exc=requests.ConnectionError("boom"))
        with pytest.raises(DriftProbeError):
            fetch_default_fingerprint()


# --- adult-ctx / dev-ctx: the pinned parent collection of each derived dataset ----------

# (probe, parent collection, raw meta.tsv grouping column), pinned as literals.
_PARENTS = [
    (fetch_adult_ctx_fingerprint, "adult-ctx-meta-atlas", "Class"),
    (fetch_dev_ctx_fingerprint, "dev-ctx-meta-atlas", "Type.v2"),
]


def _descriptor(group_column: str, *, with_group: bool = True) -> dict:
    """Minimal per-collection dataset.json: just the keys the probe reads."""
    fields = [{"name": "cellId", "label": "cellId", "type": "uniqueString"}]
    if with_group:
        fields.append(
            {
                "name": group_column.replace(".", ""),
                "label": group_column,
                "type": "enum",
                # Ordered by cell count upstream; the probe must sort. Built rows
                # follow summarize-ucsc: pandas reads NA, n/a and the empty field
                # as missing ("na" it keeps), and groups below the cut-off drop.
                "valCounts": [
                    ["Neuron", 437],
                    ["NA", 305],
                    ["n/a", 40],
                    ["na", 30],
                    ["Astrocyte", 29],
                    ["", 12],
                    ["Edge", CUT_OFF],
                    ["Rare", CUT_OFF - 1],
                ],
            }
        )
    return {
        "md5": "7cb6b243df",
        "lastBuildTime": "2025-01-08T14:59:57.970232",
        "fileVersions": {
            "outMatrix": {"md5": "1f2fbb3ff7", "size": 9321940510},
            "outMeta": {"md5": "8381a139e9", "size": 63416376},
        },
        "metaFields": fields,
    }


def _mock_parent(m, parent, group_column, *, meta_body=None, descriptor=None):
    header = f"cellId\tnCount_RNA\t{group_column}"
    m.get(
        f"https://cells.ucsc.edu/{parent}/meta.tsv",
        status_code=206,
        content=(
            meta_body if meta_body is not None else f"{header}\nc1\t5\tNeuron\n"
        ).encode("utf-8"),
    )
    m.get(
        f"https://cells.ucsc.edu/{parent}/dataset.json",
        json=descriptor if descriptor is not None else _descriptor(group_column),
    )
    return header


@pytest.mark.parametrize("probe, parent, group_column", _PARENTS)
def test_parent_fingerprint_layout(probe, parent, group_column):
    """Schema keys: the column list and the rows summarize-ucsc keeps. Extras: content."""
    with requests_mock.Mocker() as m:
        header = _mock_parent(m, parent, group_column)
        fp = probe()
        meta_request = next(r for r in m.request_history if r.path.endswith("meta.tsv"))
        assert meta_request.headers["Range"] == "bytes=0-65535"

    meta = f"{parent}/meta.tsv"
    labels = f"{meta}:{group_column}"
    assert fp["probe_version"] == 2
    assert fp["headers"] == {
        meta: ["cellId", "nCount_RNA", group_column],
        labels: ["Astrocyte", "Edge", "Neuron", "na"],
    }
    assert fp["checksums"] == {meta: hashlib.sha256(header.encode("utf-8")).hexdigest()}
    assert fp["extras"] == {
        f"{parent}/exprMatrix.tsv.gz": {"md5": "1f2fbb3ff7", "size": 9321940510},
        meta: {"md5": "8381a139e9", "size": 63416376},
        # Every published label's count, dropped ones included: routine content.
        labels: {
            "": 12,
            "Astrocyte": 29,
            "Edge": CUT_OFF,
            "NA": 305,
            "Neuron": 437,
            "Rare": CUT_OFF - 1,
            "n/a": 40,
            "na": 30,
        },
    }
    assert fp["informational"] == {
        "collection_md5": "7cb6b243df",
        "last_build_time": "2025-01-08T14:59:57.970232",
    }


def test_parent_reports_a_vanished_group_column_as_drift_not_failure():
    """A renamed grouping column is the schema change the probe exists to report."""
    probe, parent, group_column = _PARENTS[0]
    with requests_mock.Mocker() as m:
        _mock_parent(
            m,
            parent,
            group_column,
            descriptor=_descriptor(group_column, with_group=False),
        )
        fp = probe()
    assert fp["headers"][f"{parent}/meta.tsv:{group_column}"] == []
    assert fp["extras"][f"{parent}/meta.tsv:{group_column}"] == {}


def test_parent_fails_closed_on_a_truncated_header_row():
    """No newline within the range means the column list may be cut short."""
    probe, parent, group_column = _PARENTS[1]
    with requests_mock.Mocker() as m:
        _mock_parent(m, parent, group_column, meta_body="cellId\tnCount_RNA\tTy")
        with pytest.raises(DriftProbeError, match="no complete header row"):
            probe()
