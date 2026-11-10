"""AnnotationTable exports as deterministic gzipped TSV (pandas backend, no Java)."""

from __future__ import annotations

import gzip
from datetime import datetime, timezone

import click
import pandas as pd
import pytest

from hvantk.core import io as core_io
from hvantk.core.io._errors import ArtifactTypeError
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
