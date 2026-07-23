"""Tests for shared Hail helpers.

Regression cover for ``agg_max_str``: three call sites (the GeneDiseaseTableStreamer base,
plus the ClinGen and GenCC streamers) previously used ``hl.agg.max`` on ISO-date *string*
columns, which raises ``TypeError: max: parameter 'expr': expected expression of type
int32 or int64 or float32 or float64``. Hail has no string-max aggregator, so the helper
builds one; these tests pin the behaviour it has to keep, in both the grouped and
ungrouped aggregation contexts the call sites use.
"""
import hail as hl
import pytest

from hvantk.core.utils.hail_helpers import agg_max_str

pytestmark = pytest.mark.hail


def _dates_table():
    """Rows covering: a plain group, a group with one missing date, an all-missing group."""
    rows = [
        {"g": "A", "d": "2025-06-15"},
        {"g": "A", "d": "2024-01-02"},
        {"g": "B", "d": "2026-01-15"},
        {"g": "B", "d": None},
        {"g": "C", "d": None},
    ]
    return hl.Table.parallelize(rows, hl.tstruct(g=hl.tstr, d=hl.tstr))


def test_agg_max_str_ungrouped_returns_latest_iso_date():
    """The GeneDiseaseTableStreamer.compute_stats context: one fused ht.aggregate."""
    ht = _dates_table()
    assert ht.aggregate(agg_max_str(ht.d)) == "2026-01-15"


def test_agg_max_str_grouped_returns_per_group_max():
    """The ClinGen/GenCC summary context: group_by(...).aggregate(...)."""
    ht = _dates_table()
    got = {r.g: r.v for r in ht.group_by(ht.g).aggregate(v=agg_max_str(ht.d)).collect()}
    assert got["A"] == "2025-06-15"
    # The max must ignore the missing value rather than be poisoned by it.
    assert got["B"] == "2026-01-15"


def test_agg_max_str_all_missing_group_is_missing_not_an_error():
    """A group with no defined dates must yield missing, not an array-index error."""
    ht = _dates_table()
    got = {r.g: r.v for r in ht.group_by(ht.g).aggregate(v=agg_max_str(ht.d)).collect()}
    assert got["C"] is None


def test_agg_max_str_empty_table_is_missing():
    ht = hl.Table.parallelize([], hl.tstruct(g=hl.tstr, d=hl.tstr))
    assert ht.aggregate(agg_max_str(ht.d)) is None


def test_hl_agg_max_still_rejects_strings():
    """Guards the premise: if Hail ever accepted strings here, the helper could be retired."""
    ht = _dates_table()
    with pytest.raises(TypeError):
        ht.aggregate(hl.agg.max(ht.d))
