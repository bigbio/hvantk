"""Regression tests for ``_resolve_ucsc_inputs``.

The builder accepts EITHER a mapping with ``expression_matrix`` / ``metadata`` keys
(direct Phase B callers) OR a path to the raw download directory -- the ``hvantk
reprocess`` contract for a dataset with no ``lifecycle.parse`` stage, where reprocess
passes the raw dir straight to the builder. Before the fix, ``ucsc-cellbrowser:default``
crashed end-to-end because the builder assumed a dict. These are self-contained (no
network, no snapshots) so they run without the plugin's snapshot-test harness.
"""
import pytest

from hvantk.skills.ucsc_cellbrowser.builder import _resolve_ucsc_inputs
from hvantk.skills.ucsc_cellbrowser.shared.constants import (
    EXPRESSION_MATRIX_FILE_NAME,
    METADATA_FILE_NAME,
)


def test_dict_input_passthrough():
    """Direct callers keep working unchanged."""
    out = _resolve_ucsc_inputs(
        {"expression_matrix": "x/exprMatrix.tsv.gz", "metadata": "x/meta.tsv"}
    )
    assert out == ("x/exprMatrix.tsv.gz", "x/meta.tsv")


def test_raw_dir_with_accession_subdir(tmp_path):
    """reprocess passes raw_dir; the downloader nests files under <accession>/."""
    acc = tmp_path / "some-accession"
    acc.mkdir()
    expr = acc / EXPRESSION_MATRIX_FILE_NAME
    meta = acc / METADATA_FILE_NAME
    expr.write_text("gene\tcell1\n")
    meta.write_text("cell\ttype\n")
    e, m = _resolve_ucsc_inputs(str(tmp_path))
    assert e == str(expr)
    assert m == str(meta)


def test_raw_dir_flat(tmp_path):
    """Files placed directly in the raw dir (e.g. test fixtures) also resolve."""
    (tmp_path / "expression_matrix.tsv").write_text("gene\tcell1\n")
    (tmp_path / "metadata.tsv").write_text("cell\ttype\n")
    e, m = _resolve_ucsc_inputs(str(tmp_path))
    assert e.endswith("expression_matrix.tsv")
    assert m.endswith("metadata.tsv")


def test_missing_files_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        _resolve_ucsc_inputs(str(tmp_path))
