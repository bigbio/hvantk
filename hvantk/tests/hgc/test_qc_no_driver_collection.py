"""Regression: QC output must not collect the variant table to the driver.

``Table.to_pandas()`` is ``table.aggregate(hl.struct(**{c: hl.agg.collect(c)}))`` --
every column becomes one array *on the driver*. That is O(n_variants) with no ceiling,
and it does not degrade gracefully: it kills the JVM.

Measured on real data (job 19925591): a 1005-sample chr20 dense MatrixTable has ~11.1 M
variants. Collecting its 28 variant-QC fields killed a 200 GB driver after ~38 minutes,
which took ``hvantk hgc compute-qc`` down and therefore made ``hvantk hgc qc-report``
unreachable -- the whole documented workflow was unusable at cohort scale.

A fixture-sized test cannot reproduce an OOM, so these tests assert the *mechanism*
instead: that the export path never touches ``to_pandas`` on the variant table, and that
the DataFrame path is bounded. Both are cheap and both fail loudly if the collect-based
implementation comes back.
"""

import pandas as pd
import pytest

from hvantk.algorithms.hgc.qc import (
    DEFAULT_VARIANT_DF_MAX_ROWS,
    QCMetrics,
    save_qc_metrics,
)


class _FakeTable:
    """Minimal stand-in for hl.Table, counting how it gets used."""

    def __init__(self, n_rows, exploding_to_pandas=False):
        self._n = n_rows
        self._exploding = exploding_to_pandas
        self.exported_to = None
        self.export_delimiter = None
        self.sampled_p = None
        self.to_pandas_calls = 0

    def count(self):
        return self._n

    def sample(self, p, seed=None):
        self.sampled_p = p
        out = _FakeTable(max(1, int(self._n * p)), self._exploding)
        out.to_pandas_calls = 0
        self._sampled_child = out
        return out

    def export(
        self, output, types_file=None, header=True, parallel=None, delimiter="\t"
    ):
        self.exported_to = output
        self.export_delimiter = delimiter

    def to_pandas(self):
        self.to_pandas_calls += 1
        if self._exploding:
            raise AssertionError(
                "to_pandas() was called on the variant table -- this is the collect "
                "that killed a 200 GB driver on one real chromosome. Use "
                "Table.export() for the file, or bound the row count for plotting."
            )
        return pd.DataFrame({"variant_qc.call_rate": [0.9] * self._n})


class _FakeMT:
    def __init__(self):
        self.written_to = None

    def write(self, path, overwrite=False):
        self.written_to = path


def test_save_qc_metrics_exports_variants_without_collecting(tmp_path):
    """The variant CSV must come from Table.export, never to_pandas."""
    # exploding_to_pandas: any collect on the variant table fails the test outright.
    variant = _FakeTable(11_100_000, exploding_to_pandas=True)
    sample = _FakeTable(1005)
    qc = QCMetrics(mt=_FakeMT(), sample_qc=sample, variant_qc=variant)

    saved = save_qc_metrics(qc, tmp_path, prefix="t")

    assert variant.exported_to == str(tmp_path / "t_variant_qc.csv")
    assert variant.export_delimiter == ",", "the file is named .csv, so export as CSV"
    assert variant.to_pandas_calls == 0
    assert saved["variant_qc"] == str(tmp_path / "t_variant_qc.csv")


def test_save_qc_metrics_still_collects_samples(tmp_path):
    """Sample QC is bounded by cohort size, so collecting it stays correct."""
    sample = _FakeTable(1005)
    qc = QCMetrics(mt=_FakeMT(), sample_qc=sample, variant_qc=None)

    save_qc_metrics(qc, tmp_path, prefix="t")

    assert sample.to_pandas_calls == 1
    assert (tmp_path / "t_sample_qc.csv").exists()


def test_no_save_mt_does_not_write_the_matrixtable(tmp_path):
    """--no-save-mt must skip the write, not write-then-delete.

    The MatrixTable here is ~10 GB and carries every entry field; writing it and
    calling rmtree afterwards still paid for the full materialisation.
    """
    mt = _FakeMT()
    qc = QCMetrics(mt=mt, sample_qc=_FakeTable(10), variant_qc=None)

    saved = save_qc_metrics(qc, tmp_path, prefix="t", save_mt=False)

    assert mt.written_to is None, "the MatrixTable must not be written at all"
    assert "matrix_table" not in saved


def test_save_mt_true_still_writes(tmp_path):
    mt = _FakeMT()
    qc = QCMetrics(mt=mt, sample_qc=_FakeTable(10), variant_qc=None)

    saved = save_qc_metrics(qc, tmp_path, prefix="t", save_mt=True)

    assert mt.written_to == str(tmp_path / "t_with_qc.mt")
    assert saved["matrix_table"] == str(tmp_path / "t_with_qc.mt")


def test_variant_dataframe_is_bounded_above_the_budget():
    """A table larger than the budget is sampled down, and says so."""
    variant = _FakeTable(11_100_000)
    qc = QCMetrics(mt=_FakeMT(), sample_qc=None, variant_qc=variant)

    df = qc.get_variant_metrics_df(max_rows=500_000)

    assert variant.sampled_p == pytest.approx(500_000 / 11_100_000)
    assert len(df) <= 500_000
    assert df.attrs["subsampled"] is True
    assert df.attrs["n_total_variants"] == 11_100_000


def test_variant_dataframe_untouched_below_the_budget():
    """Small tables must not be sampled -- no silent loss on ordinary data."""
    variant = _FakeTable(1000)
    qc = QCMetrics(mt=_FakeMT(), sample_qc=None, variant_qc=variant)

    df = qc.get_variant_metrics_df(max_rows=500_000)

    assert variant.sampled_p is None
    assert len(df) == 1000
    assert df.attrs["subsampled"] is False


def test_counts_do_not_build_dataframes():
    """Printing a count must not collect the table (the qc-report bug)."""
    variant = _FakeTable(11_100_000, exploding_to_pandas=True)
    sample = _FakeTable(1005)
    qc = QCMetrics(mt=_FakeMT(), sample_qc=sample, variant_qc=variant)

    assert qc.count_variants() == 11_100_000
    assert qc.count_samples() == 1005
    assert variant.to_pandas_calls == 0
    assert sample.to_pandas_calls == 0


def test_default_budget_is_sane():
    assert 0 < DEFAULT_VARIANT_DF_MAX_ROWS <= 5_000_000
