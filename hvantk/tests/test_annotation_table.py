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


# --- filter tests ---

from hvantk.core.models._expr import col


def test_filter_pandas(df):
    ann = AnnotationTable.from_pandas(df, provenance=_prov())
    filtered = ann.filter(col("score") > 0.5)
    out = filtered.to_pandas()
    assert out["gene"].tolist() == ["BRCA1"]
    # filter is immutable
    assert ann.to_pandas().shape[0] == 2


@pytest.mark.hail
def test_filter_hail():
    import hail as hl

    ht = hl.Table.parallelize(
        [{"gene": "BRCA1", "score": 0.7}, {"gene": "BRCA2", "score": 0.4}],
        hl.tstruct(gene=hl.tstr, score=hl.tfloat64),
    )
    ann = AnnotationTable.from_hail(ht, provenance=_prov())
    filtered = ann.filter(col("score") > 0.5)
    out = filtered.to_pandas()
    assert out["gene"].tolist() == ["BRCA1"]


@pytest.mark.hail
def test_filter_parity(df):
    """Same predicate on both backends -> same rows."""
    import hail as hl

    pandas_ann = AnnotationTable.from_pandas(df, provenance=_prov())
    hail_ann = AnnotationTable.from_hail(
        hl.Table.from_pandas(df), provenance=_prov()
    )

    predicate = (col("score") > 0.4) & (col("gene") == "BRCA1")
    p = pandas_ann.filter(predicate).to_pandas()
    h = hail_ann.filter(predicate).to_pandas()
    assert sorted(p["gene"].tolist()) == sorted(h["gene"].tolist())


# --- select tests ---


def test_select_pandas(df):
    ann = AnnotationTable.from_pandas(df, provenance=_prov())
    out = ann.select("gene").to_pandas()
    assert list(out.columns) == ["gene"]


def test_with_columns_pandas(df):
    ann = AnnotationTable.from_pandas(df, provenance=_prov())
    out = ann.with_columns(score_sq=col("score") * col("score")).to_pandas()
    assert out["score_sq"].tolist() == pytest.approx([0.49, 0.16])


def test_rename_pandas(df):
    ann = AnnotationTable.from_pandas(df, provenance=_prov())
    out = ann.rename(gene="gene_symbol").to_pandas()
    assert "gene_symbol" in out.columns
    assert "gene" not in out.columns


@pytest.mark.hail
def test_select_hail():
    import hail as hl

    ht = hl.Table.parallelize(
        [{"gene": "BRCA1", "score": 0.7, "extra": 1}],
        hl.tstruct(gene=hl.tstr, score=hl.tfloat64, extra=hl.tint32),
    )
    ann = AnnotationTable.from_hail(ht, provenance=_prov())
    out = ann.select("gene", "score").to_pandas()
    assert set(out.columns) == {"gene", "score"}


@pytest.mark.hail
def test_with_columns_hail():
    import hail as hl

    ht = hl.Table.parallelize(
        [{"gene": "BRCA1", "score": 0.7}],
        hl.tstruct(gene=hl.tstr, score=hl.tfloat64),
    )
    ann = AnnotationTable.from_hail(ht, provenance=_prov())
    out = ann.with_columns(score_sq=col("score") * col("score")).to_pandas()
    assert out["score_sq"].iloc[0] == pytest.approx(0.49)


# --- count, collect, head, distinct, join tests ---


def test_count_collect_head_pandas(df):
    ann = AnnotationTable.from_pandas(df, provenance=_prov())
    assert ann.count() == 2
    rows = ann.collect()
    assert rows == [
        {"gene": "BRCA1", "score": 0.7},
        {"gene": "BRCA2", "score": 0.4},
    ]
    h = ann.head(1).collect()
    assert h == [{"gene": "BRCA1", "score": 0.7}]


def test_distinct_pandas():
    df = pd.DataFrame({"x": [1, 1, 2], "y": ["a", "a", "b"]})
    ann = AnnotationTable.from_pandas(df, provenance=_prov())
    assert ann.distinct().count() == 2
    assert ann.distinct(subset=["x"]).count() == 2


def test_join_pandas():
    a = AnnotationTable.from_pandas(
        pd.DataFrame({"gene": ["BRCA1", "BRCA2"], "score": [0.7, 0.4]}),
        provenance=_prov(),
    )
    b = AnnotationTable.from_pandas(
        pd.DataFrame({"gene": ["BRCA1", "TP53"], "chrom": ["chr17", "chr17"]}),
        provenance=_prov(),
    )
    joined = a.join(b, on="gene", how="inner").collect()
    assert joined == [{"gene": "BRCA1", "score": 0.7, "chrom": "chr17"}]


@pytest.mark.hail
def test_count_collect_head_hail():
    import hail as hl

    ht = hl.Table.parallelize(
        [{"gene": "BRCA1", "score": 0.7}, {"gene": "BRCA2", "score": 0.4}],
        hl.tstruct(gene=hl.tstr, score=hl.tfloat64),
    )
    ann = AnnotationTable.from_hail(ht, provenance=_prov())
    assert ann.count() == 2
    rows = ann.collect()
    assert {r["gene"] for r in rows} == {"BRCA1", "BRCA2"}
    assert len(ann.head(1).collect()) == 1


@pytest.mark.hail
def test_join_hail():
    import hail as hl

    a = AnnotationTable.from_hail(
        hl.Table.parallelize(
            [{"gene": "BRCA1", "score": 0.7}, {"gene": "BRCA2", "score": 0.4}],
            hl.tstruct(gene=hl.tstr, score=hl.tfloat64),
        ).key_by("gene"),
        provenance=_prov(),
    )
    b = AnnotationTable.from_hail(
        hl.Table.parallelize(
            [{"gene": "BRCA1", "chrom": "chr17"}],
            hl.tstruct(gene=hl.tstr, chrom=hl.tstr),
        ).key_by("gene"),
        provenance=_prov(),
    )
    joined = a.join(b, on="gene", how="inner").collect()
    assert len(joined) == 1
    assert joined[0]["gene"] == "BRCA1"
