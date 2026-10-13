"""Tests for hvantk.core.models.anndata_utils."""

import os
from unittest.mock import MagicMock

import anndata as ad
import click
import numpy as np
import pandas as pd
import pytest


def _make_test_adata() -> ad.AnnData:
    """Create a small AnnData for testing (20 obs, 10 var)."""
    np.random.seed(42)
    X = np.random.rand(20, 10)
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(np.random.choice(["A", "B", "C"], size=20)),
            "n_counts": np.random.rand(20) * 1000,
        },
        index=[f"cell_{i}" for i in range(20)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(10)])
    return ad.AnnData(X=X, obs=obs, var=var)


class TestAnnotateColumnSummary:
    def test_adds_column_summary_to_uns(self):
        from hvantk.core.models.anndata_utils import annotate_column_summary_ad

        adata = _make_test_adata()
        annotate_column_summary_ad(adata)
        assert "column_summary" in adata.uns

    def test_numeric_column(self):
        from hvantk.core.models.anndata_utils import annotate_column_summary_ad

        adata = _make_test_adata()
        annotate_column_summary_ad(adata)
        summary = adata.uns["column_summary"]
        n_counts = summary["n_counts"]
        assert n_counts["dtype"] == "numeric"
        assert "min" in n_counts
        assert "max" in n_counts
        assert "mean" in n_counts
        assert "n_missing" in n_counts

    def test_categorical_column(self):
        from hvantk.core.models.anndata_utils import annotate_column_summary_ad

        adata = _make_test_adata()
        annotate_column_summary_ad(adata)
        summary = adata.uns["column_summary"]
        ct = summary["cell_type"]
        assert ct["dtype"] == "categorical"
        assert ct["n_unique"] == 3
        assert sorted(ct["levels"]) == ["A", "B", "C"]
        assert "n_missing" in ct

    def test_high_cardinality_categorical(self):
        from hvantk.core.models.anndata_utils import annotate_column_summary_ad

        adata = _make_test_adata()
        # Add a high-cardinality column
        adata.obs["cell_id"] = [f"id_{i}" for i in range(20)]
        annotate_column_summary_ad(adata, max_levels=5, top_n_levels=3)
        summary = adata.uns["column_summary"]
        cid = summary["cell_id"]
        assert cid["dtype"] == "categorical"
        assert cid["n_unique"] == 20
        assert "top_levels" in cid
        assert len(cid["top_levels"]) <= 3
        assert "levels" not in cid


class TestSaveLoadAnndata:
    def test_roundtrip(self, tmp_path):
        from hvantk.core.io.anndata_io import load_anndata, save_anndata

        adata = _make_test_adata()
        path = str(tmp_path / "test.h5ad")
        save_anndata(adata, path)
        loaded = load_anndata(path)
        assert loaded.shape == adata.shape
        np.testing.assert_array_almost_equal(loaded.X, adata.X)

    def test_overwrite_false_raises(self, tmp_path):
        from hvantk.core.io.anndata_io import save_anndata

        adata = _make_test_adata()
        path = str(tmp_path / "test.h5ad")
        save_anndata(adata, path)
        with pytest.raises(FileExistsError):
            save_anndata(adata, path, overwrite=False)

    def test_overwrite_true_succeeds(self, tmp_path):
        from hvantk.core.io.anndata_io import save_anndata

        adata = _make_test_adata()
        path = str(tmp_path / "test.h5ad")
        save_anndata(adata, path)
        save_anndata(adata, path, overwrite=True)
        assert os.path.exists(path)


# --- Regression: issue #383 -------------------------------------------------------
# A transcript-level Expression Atlas build indexes `var` by the transcript id
# (var.index.name == "transcript_id"); every generic AnnData consumer that reads
# var_names as the gene id must refuse one outright rather than silently emitting
# per-transcript results labelled as per-gene. One fixture shape per consumer --
# `hvantk annotate --matrix` consumes an already-summarized (groups x genes) AnnData,
# the other three consume raw (cells x genes) AnnData -- but all share the same
# transcript-keyed var index that require_gene_keyed is supposed to catch.


def _transcript_keyed_cells_adata() -> ad.AnnData:
    """Raw cells/samples x transcripts AnnData, as markers/summarize/ptm-constraint see it."""
    n_obs, n_vars = 4, 3
    X = (np.arange(n_obs * n_vars, dtype=np.float32) + 1.0).reshape(n_obs, n_vars)
    obs = pd.DataFrame(
        {"grp": ["a", "a", "b", "b"]}, index=[f"s{i}" for i in range(n_obs)]
    )
    var = pd.DataFrame(
        index=pd.Index([f"t{i}" for i in range(n_vars)], name="transcript_id")
    )
    return ad.AnnData(X=X, obs=obs, var=var)


def _transcript_keyed_summary_adata() -> ad.AnnData:
    """Groups x transcripts AnnData with a 'mean' layer, as reduce_matrix_to_gene expects."""
    var = pd.DataFrame(index=pd.Index(["t0", "t1"], name="transcript_id"))
    mean = np.array([[1.0, 2.0], [3.0, 4.0]])
    return ad.AnnData(
        X=None, obs=pd.DataFrame(index=["a", "b"]), var=var, layers={"mean": mean}
    )


def _run_markers(tmp_path, monkeypatch):
    import hvantk.algorithms.expression.matrix_utils as matrix_utils_mod
    import hvantk.tools.expression.summarize_expression_cli as cli_mod

    # markers_cmd does `from ...matrix_utils import require_scanpy` INSIDE the function
    # body, so the patch target is the source module, not this one -- the guard must
    # fire before scanpy is ever touched, and stubbing it out means the test does not
    # need the optional 'expression' extra installed to exercise that ordering.
    monkeypatch.setattr(matrix_utils_mod, "require_scanpy", lambda: MagicMock())
    h5ad_path = tmp_path / "transcript.h5ad"
    _transcript_keyed_cells_adata().write_h5ad(h5ad_path)
    cli_mod.markers_cmd.callback(
        matrix_path=str(h5ad_path),
        group_by="grp",
        filter_by=(),
        method="wilcoxon",
        top_n=10,
        output=str(tmp_path / "out.json"),
        results_tsv=None,
        overwrite=False,
    )


def _run_summarize(tmp_path, monkeypatch):
    import hvantk.algorithms.expression.matrix_utils as matrix_utils_mod
    import hvantk.tools.expression.summarize_expression_cli as cli_mod

    monkeypatch.setattr(matrix_utils_mod, "require_scanpy", lambda: MagicMock())
    h5ad_path = tmp_path / "transcript.h5ad"
    _transcript_keyed_cells_adata().write_h5ad(h5ad_path)
    cli_mod.summarize_expression_cmd.callback(
        matrix_path=str(h5ad_path),
        group_by=("grp",),
        filter_by=(),
        min_cells=1,
        output=str(tmp_path / "out.h5ad"),
        overwrite=False,
    )


def _run_ptm_constraint(tmp_path, monkeypatch):
    from hvantk.algorithms.ptm.constraint_expression import load_gene_by_group_matrix

    h5ad_path = tmp_path / "transcript.h5ad"
    _transcript_keyed_cells_adata().write_h5ad(h5ad_path)
    load_gene_by_group_matrix(
        "anndata", str(h5ad_path), grouping="grp", min_cells_per_group=1
    )


def _run_annotate_matrix(tmp_path, monkeypatch):
    from hvantk.algorithms.annotation.matrix import reduce_matrix_to_gene
    from hvantk.algorithms.annotation.spec import MatrixSpec

    reduce_matrix_to_gene(
        _transcript_keyed_summary_adata(), MatrixSpec(group_axis="grp", atlas="asp")
    )


@pytest.mark.parametrize(
    "consumer, run",
    [
        ("hvantk expression markers", _run_markers),
        ("hvantk expression summarize", _run_summarize),
        ("ptm constraint (anndata source)", _run_ptm_constraint),
        ("hvantk annotate --matrix", _run_annotate_matrix),
    ],
)
def test_require_gene_keyed_rejects_transcript_keyed_adata_everywhere(
    consumer, run, tmp_path, monkeypatch
):
    """Every gene-level consumer refuses a transcript-keyed AnnData, naming the fix.

    `click.ClickException` covers the CLI consumers (markers/summarize), which catch
    `require_gene_keyed`'s `ValueError` and re-raise it as one; its `str()` is the same
    message either way.
    """
    with pytest.raises((ValueError, click.ClickException), match="level=gene"):
        run(tmp_path, monkeypatch)
