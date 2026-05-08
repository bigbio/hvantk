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
def test_to_jsonable_renders_hail_struct_as_nested_dict(hail_session):
    import hail as hl
    from hvantk.tests._snapshot_utils import _to_jsonable

    s = hl.Struct(a=1, b="x", nested=hl.Struct(c=[1, 2]))
    assert _to_jsonable(s) == {"a": 1, "b": "x", "nested": {"c": [1, 2]}}


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


def test_anndata_schema_to_dict():
    import anndata as ad
    import numpy as np
    import pandas as pd
    from hvantk.tests._snapshot_utils import anndata_schema_to_dict

    X = np.zeros((3, 4), dtype=np.float32)
    obs = pd.DataFrame({"sample_id": ["s1", "s2", "s3"], "celltype": ["A", "B", "A"]})
    obs.index = obs["sample_id"]
    var = pd.DataFrame({"gene_id": ["g1", "g2", "g3", "g4"]})
    var.index = var["gene_id"]
    a = ad.AnnData(X=X, obs=obs, var=var)

    schema = anndata_schema_to_dict(a)
    assert schema["n_obs"] == 3
    assert schema["n_vars"] == 4
    assert "sample_id" in schema["obs_columns"]
    assert "celltype" in schema["obs_columns"]
    assert "gene_id" in schema["var_columns"]
    assert schema["X_dtype"] == "float32"
    assert "ndarray" in schema["X_format"] or "array" in schema["X_format"].lower()
    assert schema["layers"] == []


def test_anndata_sample_rows():
    import anndata as ad
    import numpy as np
    import pandas as pd
    from hvantk.tests._snapshot_utils import anndata_sample_rows

    X = np.arange(12, dtype=np.float32).reshape(3, 4)
    obs = pd.DataFrame({"celltype": ["A", "B", "C"]}, index=["c1", "c2", "c3"])
    var = pd.DataFrame({"gene_id": ["g1", "g2", "g3", "g4"]}, index=["g1", "g2", "g3", "g4"])
    a = ad.AnnData(X=X, obs=obs, var=var)

    rows = anndata_sample_rows(a, n=2)
    assert len(rows["obs_head"]) == 2
    assert rows["obs_head"][0]["celltype"] == "A"
    assert len(rows["var_head"]) == 2
    assert rows["var_head"][0]["gene_id"] == "g1"
    assert rows["X_corner"] == [[0.0, 1.0], [4.0, 5.0]]


def test_regenerate_snapshots_dispatches_on_anndata(tmp_path):
    """Builders returning AnnData write anndata-shape snapshots, not Hail-shape."""
    import json
    import anndata as ad
    import numpy as np
    import pandas as pd
    from hvantk.tests._snapshot_utils import regenerate_snapshots

    def fake_builder(expr_path: str, output_path: str):
        X = np.zeros((2, 3), dtype=np.float32)
        obs = pd.DataFrame({"celltype": ["A", "B"]}, index=["c1", "c2"])
        var = pd.DataFrame({"gene_id": ["g1", "g2", "g3"]}, index=["g1", "g2", "g3"])
        return ad.AnnData(X=X, obs=obs, var=var)

    snapshot_dir = tmp_path / "snap"
    regenerate_snapshots(
        builder_fn=fake_builder,
        fixture_path="ignored.tsv",
        snapshot_dir=snapshot_dir,
        builder_kwargs={"output_path": str(tmp_path / "out.h5ad")},
        input_path_kwarg="expr_path",
    )
    schema = json.loads((snapshot_dir / "schema.json").read_text())
    assert schema["n_obs"] == 2
    assert schema["n_vars"] == 3
    rows = json.loads((snapshot_dir / "sample_rows.json").read_text())
    assert "obs_head" in rows
    assert "var_head" in rows
    assert "X_corner" in rows
