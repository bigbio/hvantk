import json
from pathlib import Path

import pytest

from hvantk.tests._snapshot_utils import (
    hail_schema_to_dict,
    collect_sample_rows,
    load_snapshot,
)


@pytest.mark.hail
def test_hail_schema_to_dict_roundtrips_table(hail_session):
    import hail as hl

    ht = hl.utils.range_table(3).annotate(value=hl.str("hello"))
    schema = hail_schema_to_dict(ht)
    assert schema["key"] == ["idx"]
    assert schema["row"]["idx"] == "int32"
    assert schema["row"]["value"] == "str"
    # JSON-roundtrippable
    json.loads(json.dumps(schema))


@pytest.mark.hail
def test_collect_sample_rows_matches_by_key(hail_session):
    import hail as hl

    ht = hl.utils.range_table(5).annotate(value=hl.str("v"))
    rows = collect_sample_rows(ht, keys=[{"idx": 0}, {"idx": 2}])
    assert len(rows) == 2
    assert rows[0]["key"] == {"idx": 0}
    assert rows[1]["key"] == {"idx": 2}


def test_load_snapshot_reads_json(tmp_path):
    p = tmp_path / "schema.json"
    p.write_text(json.dumps({"key": ["idx"], "row": {"idx": "int32"}}))
    data = load_snapshot(p)
    assert data == {"key": ["idx"], "row": {"idx": "int32"}}


def test_to_jsonable_handles_primitives():
    from hvantk.tests._snapshot_utils import _to_jsonable

    assert _to_jsonable(None) is None
    assert _to_jsonable("hello") == "hello"
    assert _to_jsonable(42) == 42
    assert _to_jsonable(3.14) == 3.14
    assert _to_jsonable(True) is True


def test_to_jsonable_handles_collections():
    from hvantk.tests._snapshot_utils import _to_jsonable

    assert _to_jsonable([1, 2, 3]) == [1, 2, 3]
    assert _to_jsonable((1, 2)) == [1, 2]
    assert _to_jsonable({"b": 2, "a": 1}) == {"b": 2, "a": 1}
    assert _to_jsonable({3, 1, 2}) == [1, 2, 3]
    assert _to_jsonable([{1, 2}, {3, 4}]) == [[1, 2], [3, 4]]


def test_to_jsonable_falls_back_to_str_for_unknown():
    from hvantk.tests._snapshot_utils import _to_jsonable

    class Weird:
        def __str__(self):
            return "weird-repr"

    assert _to_jsonable(Weird()) == "weird-repr"


@pytest.mark.hail
def test_to_jsonable_renders_hail_locus(hail_session):
    import hail as hl
    from hvantk.tests._snapshot_utils import _to_jsonable

    locus = hl.Locus("chr1", 12345, reference_genome="GRCh38")
    assert _to_jsonable(locus) == "chr1:12345"


@pytest.mark.hail
def test_collect_sample_rows_raises_on_missing_key(hail_session):
    import hail as hl
    from hvantk.tests._snapshot_utils import collect_sample_rows

    ht = hl.utils.range_table(3)
    with pytest.raises(KeyError, match="not found"):
        collect_sample_rows(ht, keys=[{"idx": 999}])


@pytest.mark.hail
def test_collect_sample_rows_handles_array_keys(hail_session):
    """Variant tables key on (locus, alleles) where alleles is array<str>.

    Regression: list values inside the key tuple must be coerced to tuples
    so the dict-key path is hashable.
    """
    import hail as hl
    from hvantk.tests._snapshot_utils import collect_sample_rows

    ht = hl.Table.parallelize(
        [
            {"locus": hl.Locus("chr1", 100, reference_genome="GRCh38"), "alleles": ["A", "G"], "score": 1},
            {"locus": hl.Locus("chr1", 200, reference_genome="GRCh38"), "alleles": ["C", "T"], "score": 2},
        ],
        schema=hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            score=hl.tint32,
        ),
        key=["locus", "alleles"],
    )
    rows = collect_sample_rows(ht, keys=[{"locus": "chr1:100", "alleles": ["A", "G"]}])
    assert len(rows) == 1
    assert rows[0]["key"] == {"locus": "chr1:100", "alleles": ["A", "G"]}
    assert rows[0]["row"] == {"score": 1}
