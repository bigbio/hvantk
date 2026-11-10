"""AnnotationTable exports as deterministic gzipped TSV (pandas backend, no Java)."""

from __future__ import annotations

import gzip
from datetime import datetime, timezone

import click
import pandas as pd
import pytest

from hvantk.core import io as core_io
from hvantk.core.io._errors import ArtifactTypeError
from hvantk.core.io._formats import _tsv_sort_rows
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


def _table(df):
    return AnnotationTable.from_pandas(df, provenance=_prov())


def _text(path):
    with gzip.open(path, "rb") as fh:
        return fh.read().decode("utf-8")


def test_tsv_save_is_deterministic_and_joins_arrays(tmp_path):
    df = pd.DataFrame(
        {
            "gene": ["B", "A", "C", "D"],
            "tags": [["A", "B"], [], None, {"y", "x"}],
            "flag": [True, False, True, False],
            "n": pd.array([5, None, 7, 0], dtype="Int64"),
            "score": [1.5, None, 2.0, 0.25],
            "note": ["é", "plain", None, ""],
        }
    )
    ann = _table(df)
    first, second = tmp_path / "a.tsv.gz", tmp_path / "b.tsv.gz"
    core_io.save(ann, first)
    core_io.save(ann, second)

    raw = first.read_bytes()
    assert raw == second.read_bytes()
    # gzip header: MTIME (bytes 4-8) zero and no FNAME flag, so no clock or path
    assert raw[4:8] == b"\0\0\0\0" and not raw[3] & 0x08
    assert first.with_name("a.tsv.gz.provenance.json").exists()
    # sorted by every field left to right; arrays "|"-joined (sets sorted), empty
    # array and missing values empty, bools true/false, ints plain decimal.
    assert _text(first) == (
        "gene\ttags\tflag\tn\tscore\tnote\n"
        "A\t\tfalse\t\t\tplain\n"
        "B\tA|B\ttrue\t5\t1.5\té\n"
        "C\t\ttrue\t7\t2.0\t\n"
        "D\tx|y\tfalse\t0\t0.25\t\n"
    )


def test_tsv_float_column_holding_whole_numbers_for_an_int_schema_prints_plain(
    tmp_path,
):
    ann = _table(pd.DataFrame({"g": ["a", "b"], "n": [5, 6]}))
    # a missing value in a numpy int column makes pandas hold floats
    ann._table["n"] = pd.Series([5.0, float("nan")])
    out = tmp_path / "t.tsv.gz"
    core_io.save(ann, out)
    assert _text(out) == "g\tn\na\t5\nb\t\n"


@pytest.mark.parametrize("bad", ["a\tb", "a\nb", "a\rb"])
def test_tsv_save_refuses_a_tab_in_a_cell(tmp_path, bad):
    ann = _table(pd.DataFrame({"gene": ["G1", "G2"], "note": ["ok", bad]}))
    out = tmp_path / "t.tsv.gz"
    with pytest.raises(ValueError, match=r"note.*G2"):
        core_io.save(ann, out)
    assert not out.exists()


def test_tsv_cannot_be_loaded(tmp_path):
    out = tmp_path / "t.tsv.gz"
    core_io.save(_table(pd.DataFrame({"g": ["a"]})), out)
    with pytest.raises(ArtifactTypeError, match="export"):
        core_io.load(out)


def test_reprocess_accepts_tsv_gz_for_an_annotation_table():
    from types import SimpleNamespace

    from hvantk.tools.plugins.reprocess_cli import (
        _check_output_extension,
        _expected_extensions,
    )

    for backend, native in (("pandas", ".parquet"), ("hail", ".ht")):
        spec = SimpleNamespace(
            name="x:y",
            backend=backend,
            artifact_type=SimpleNamespace(__name__="AnnotationTable"),
        )
        assert _expected_extensions(spec) == (native, ".tsv.gz")
        _check_output_extension(spec, "out.tsv.gz")
        _check_output_extension(spec, "out" + native)
        with pytest.raises(click.UsageError):
            _check_output_extension(spec, "out.csv")
    vm = SimpleNamespace(
        name="x:y",
        backend="hail",
        artifact_type=SimpleNamespace(__name__="VariantMatrix"),
    )
    with pytest.raises(click.UsageError):
        _check_output_extension(vm, "out.tsv.gz")


def test_tsv_big_ints_in_an_all_numeric_table_keep_every_digit(tmp_path):
    big = 2**60 + 1
    ann = _table(pd.DataFrame({"n": [big + 2, big], "s": [0.5, 1.5]}))
    out = tmp_path / "t.tsv.gz"
    core_io.save(ann, out)
    assert _text(out) == f"n\ts\n{big}\t1.5\n{big + 2}\t0.5\n"


def test_tsv_nullable_int_column_prints_integers_and_an_empty_cell(tmp_path):
    ann = _table(
        pd.DataFrame({"g": ["a", "b", "c"], "n": pd.array([3, pd.NA, 12], "Int64")})
    )
    assert ann.schema["n"] == "int"
    out = tmp_path / "t.tsv.gz"
    core_io.save(ann, out)
    assert _text(out) == "g\tn\na\t3\nb\t\nc\t12\n"


def test_tsv_sorts_ties_by_the_next_field_with_missing_last(tmp_path):
    df = pd.DataFrame(
        {"a": ["x", "x", "x", "w"], "b": [None, "z", "y", "q"], "v": [1, 2, 3, 4]}
    )
    out = tmp_path / "t.tsv.gz"
    core_io.save(_table(df), out)
    assert _text(out) == "a\tb\tv\nw\tq\t4\nx\ty\t3\nx\tz\t2\nx\t\t1\n"


def test_tsv_sort_puts_key_fields_first_and_ignores_absent_key_names():
    cols = {"v": [1, 2, 3], "k": ["b", "a", "a"], "locus.pos": [9, 8, 7]}
    # key "locus" was flattened away: ignored; "k" leads, then v breaks the tie
    assert _tsv_sort_rows(cols, ["locus", "k"], {}) == [1, 2, 0]
    assert _tsv_sort_rows(cols, [], {}) == [0, 1, 2]


def test_tsv_float32_nan_is_a_missing_cell(tmp_path):
    import numpy as np

    df = pd.DataFrame({"g": ["a"], "x": np.array([np.nan], dtype="float32")})
    out = tmp_path / "t.tsv.gz"
    core_io.save(_table(df), out)
    assert _text(out) == "g\tx\na\t\n"


@pytest.mark.parametrize(
    "bad, why", [(["a|b"], r"contains '\|'"), ([["a"], "b"], "nested array")]
)
def test_tsv_refuses_ambiguous_array_cells(tmp_path, bad, why):
    ann = _table(pd.DataFrame({"gene": ["G1", "G2"], "tags": [["ok"], bad]}))
    out = tmp_path / "t.tsv.gz"
    with pytest.raises(ValueError, match=rf"tags.*G2.*{why}"):
        core_io.save(ann, out)
    assert not out.exists()


def test_tsv_empty_array_missing_and_blank_element_write_the_same_cell(tmp_path):
    df = pd.DataFrame({"g": ["a", "b", "c"], "t": [[], None, [""]]})
    out = tmp_path / "t.tsv.gz"
    core_io.save(_table(df), out)
    assert _text(out) == "g\tt\na\t\nb\t\nc\t\n"
