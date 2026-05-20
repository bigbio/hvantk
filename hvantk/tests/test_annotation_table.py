"""Tests for AnnotationTable construction, conversion, identity."""
from __future__ import annotations

from datetime import datetime, timezone

import pandas as pd
import pytest

from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.provenance import Provenance


def _prov():
    return Provenance(
        plugin="t",
        dataset="t:rows",
        plugin_version="0.0",
        source_fingerprint="sha256:abc",
        schema_id="t-rows-v1",
        build_timestamp=datetime(2026, 5, 20, tzinfo=timezone.utc),
        builder_commit=None,
    )


@pytest.fixture()
def df():
    return pd.DataFrame(
        {"gene": ["BRCA1", "BRCA2"], "score": [0.7, 0.4]}
    )


def test_from_pandas_round_trip(df):
    ann = AnnotationTable.from_pandas(df, provenance=_prov())
    assert ann.backend == "pandas"
    assert ann.provenance == _prov()
    assert ann.schema == {"gene": "str", "score": "float"}
    out = ann.to_pandas()
    pd.testing.assert_frame_equal(out.reset_index(drop=True), df.reset_index(drop=True))


def test_from_pandas_requires_provenance(df):
    with pytest.raises(TypeError):
        AnnotationTable.from_pandas(df)  # provenance kwarg required


@pytest.mark.hail
def test_from_hail_round_trip():
    import hail as hl

    ht = hl.Table.parallelize(
        [{"gene": "BRCA1", "score": 0.7}],
        hl.tstruct(gene=hl.tstr, score=hl.tfloat64),
    )
    ann = AnnotationTable.from_hail(ht, provenance=_prov())
    assert ann.backend == "hail"
    assert ann.provenance == _prov()
    out = ann.to_hail()
    assert out.count() == 1


@pytest.mark.hail
def test_to_pandas_works_on_hail_backend():
    import hail as hl

    ht = hl.Table.parallelize(
        [{"gene": "BRCA1", "score": 0.7}],
        hl.tstruct(gene=hl.tstr, score=hl.tfloat64),
    )
    ann = AnnotationTable.from_hail(ht, provenance=_prov())
    df = ann.to_pandas()
    assert df.iloc[0]["gene"] == "BRCA1"
