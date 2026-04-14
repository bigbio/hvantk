"""Smoke tests for the PTM constraint pipeline (M1-M6).

Minimum checks that ``hvantk.utils.tissue_specificity``,
``hvantk.ptm.constraint_expression``, ``hvantk.ptm.constraint``, and the
``hvantk ptm constraint`` CLI subcommand import, agree with a reference
implementation, and respect their contracts. No Hail, no network.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner


def _toy_expr_matrix() -> pd.DataFrame:
    """5 genes x 6 tissues fixture with mixed broad / tissue-specific patterns."""
    rng = np.random.default_rng(0)
    tissues = [f"t{i}" for i in range(6)]
    genes = [f"g{i}" for i in range(5)]
    # Broad gene, tissue-specific gene, 3 noisy rows.
    data = np.vstack([
        np.full(6, 10.0),                      # housekeeping
        np.array([0.1, 0.1, 0.1, 0.1, 0.1, 50.0]),  # highly specific
        rng.uniform(1.0, 8.0, size=6),
        rng.uniform(0.5, 5.0, size=6),
        np.array([2.0, 0.2, 0.2, 0.2, 0.2, 0.2]),   # specific to t0
    ])
    return pd.DataFrame(data, index=genes, columns=tissues)


def test_tspex_wrapper_matches_yanai():
    from hvantk.utils.tissue_specificity import compute_specificity, tau_yanai_reference

    df = _toy_expr_matrix()
    tspex_tau = compute_specificity(df, method="tau")
    ref_tau = tau_yanai_reference(df)
    for gene in df.index:
        assert tspex_tau[gene] == pytest.approx(ref_tau[gene], abs=0.01)


def test_tabular_adapter_contract():
    from hvantk.ptm.constraint_expression import load_gene_by_group_matrix

    df = _toy_expr_matrix()
    with tempfile.TemporaryDirectory() as tmp:
        tsv_path = Path(tmp) / "expr.tsv"
        df.to_csv(tsv_path, sep="\t")
        out = load_gene_by_group_matrix("tabular", str(tsv_path), grouping="x")

    assert out.shape == df.shape
    assert out.index.name == "gene_id"
    assert not out.isna().any().any()
    assert (out.values >= 0).all()


def test_anndata_adapter_contract():
    import anndata as ad
    from hvantk.ptm.constraint_expression import load_gene_by_group_matrix

    n_cells = 6
    n_genes = 4
    X = np.abs(np.random.default_rng(42).normal(5.0, 1.0, size=(n_cells, n_genes)))
    obs = pd.DataFrame(
        {"tissue": ["liver", "liver", "brain", "brain", "kidney", "kidney"]},
        index=[f"c{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    with tempfile.TemporaryDirectory() as tmp:
        h5ad_path = Path(tmp) / "tiny.h5ad"
        adata.write_h5ad(h5ad_path)
        out = load_gene_by_group_matrix(
            "anndata",
            str(h5ad_path),
            grouping="tissue",
            min_cells_per_group=1,
        )

    assert out.shape[0] == n_genes
    assert set(out.columns) == {"liver", "brain", "kidney"}
    assert out.index.name == "gene_id"


def test_constraint_config_validate():
    from hvantk.ptm.constraint import PTMConstraintConfig

    cfg = PTMConstraintConfig(
        variants_ht_path="/does/not/exist.ht",
        expression_source="bogus-source",
        expression_path="/also/missing.tsv",
        grouping="tissue",
        output_dir="/tmp/out",
        label_filter="NOT_A_LABEL",
        min_cells_per_group=-5,
        min_variants_per_group=0,
    )
    errors = cfg.validate()
    assert errors, "Expected validation errors for bad config"
    joined = " ".join(errors)
    assert "expression-source" in joined
    assert "label-filter" in joined
    assert "variants-ht" in joined


def test_gene_features_compute():
    from hvantk.ptm.constraint import _compute_gene_features

    df = _toy_expr_matrix()
    features = _compute_gene_features(df, expressed_threshold=1.0, gene_id_map=None)

    assert list(features.columns) == ["tau", "primary_group", "max_expr", "expressed"]
    assert len(features) == len(df)
    # Row 0 (housekeeping) has all equal values, idxmax returns the first column.
    assert features.loc["g0", "primary_group"] == df.columns[0]
    # Row 1 peaks at t5.
    assert features.loc["g1", "primary_group"] == "t5"
    # Row 4 peaks at t0.
    assert features.loc["g4", "primary_group"] == "t0"
    assert bool(features.loc["g0", "expressed"]) is True


def test_cli_help_renders():
    from hvantk.commands.ptm_cli import ptm_group

    runner = CliRunner()
    result = runner.invoke(ptm_group, ["constraint", "--help"])
    assert result.exit_code == 0, result.output
    assert "--variants-ht" in result.output
