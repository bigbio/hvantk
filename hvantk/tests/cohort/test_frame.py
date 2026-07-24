"""The pandas path onto a cohort manifest -- the ``attach()`` sibling for ``rerank``.

``rerank`` is pandas + scikit-learn with zero Hail, so a ``CohortManifest`` reaches it
through :mod:`hvantk.algorithms.cohort.frame` rather than through
:func:`hvantk.algorithms.cohort.attach.attach`. Pure Python -- no Hail -- so these run
in the default fast suite.
"""
import gzip

import pytest

from hvantk.algorithms.cohort.frame import load_cohort_frame, load_prior_frame
from hvantk.algorithms.cohort.spec import CohortAxis, CohortManifest, CohortPrior


def _manifest(table, key="gene", key_column=None, prior_col="minp", axes=()):
    return CohortManifest(
        name="demo",
        key=key,
        key_column=key_column,
        table=str(table),
        prior=CohortPrior(column=prior_col, direction="lower_is_better"),
        cohort_axes=axes,
    )


def _write_tsv(path, header, rows):
    lines = ["\t".join(header)]
    lines.extend("\t".join(str(v) for v in row) for row in rows)
    path.write_text("\n".join(lines) + "\n")


def test_load_cohort_frame_reads_a_plain_tsv_and_renames_key_column(tmp_path):
    """key='symbol' names the identifier SPACE; key_column='symbol_col' names the
    actual column in the table -- the common real-world case where the two differ."""
    p = tmp_path / "cohort.tsv"
    _write_tsv(
        p,
        ["symbol_col", "minp", "n_case_var"],
        [("A", 0.01, 5), ("B", 0.20, 3)],
    )
    m = _manifest(
        p,
        key="symbol",
        key_column="symbol_col",
        axes=(CohortAxis(axis="burden", columns=("n_case_var",)),),
    )

    frame = load_cohort_frame(m)

    assert "gene" in frame.columns
    assert "symbol_col" not in frame.columns
    assert set(frame["gene"]) == {"A", "B"}


def test_load_cohort_frame_returns_exactly_declared_columns_plus_gene(tmp_path):
    p = tmp_path / "cohort.tsv"
    _write_tsv(
        p,
        ["gene", "minp", "n_case_var", "unused_col"],
        [("A", 0.01, 5, "x"), ("B", 0.20, 3, "y")],
    )
    m = _manifest(p, axes=(CohortAxis(axis="burden", columns=("n_case_var",)),))

    frame = load_cohort_frame(m)

    assert set(frame.columns) == {"gene", "minp", "n_case_var"}


def test_load_cohort_frame_raises_on_missing_declared_column(tmp_path):
    p = tmp_path / "cohort.tsv"
    _write_tsv(p, ["gene", "minp"], [("A", 0.01)])
    m = _manifest(p, axes=(CohortAxis(axis="burden", columns=("n_case_var",)),))

    with pytest.raises(ValueError, match="n_case_var"):
        load_cohort_frame(m)


def test_load_cohort_frame_raises_on_missing_key_column(tmp_path):
    p = tmp_path / "cohort.tsv"
    _write_tsv(p, ["symbol", "minp"], [("A", 0.01)])
    m = _manifest(p, key="gene", key_column="gene")

    with pytest.raises(ValueError, match="gene"):
        load_cohort_frame(m)


def test_load_cohort_frame_raises_on_duplicate_gene_and_names_it(tmp_path):
    p = tmp_path / "cohort.tsv"
    _write_tsv(p, ["gene", "minp"], [("A", 0.01), ("B", 0.20), ("A", 0.50)])
    m = _manifest(p)

    with pytest.raises(ValueError, match="A"):
        load_cohort_frame(m)


def test_load_cohort_frame_reads_a_gzipped_tsv(tmp_path):
    p = tmp_path / "cohort.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("gene\tminp\tn_case_var\nA\t0.01\t5\nB\t0.20\t3\n")
    m = _manifest(p, axes=(CohortAxis(axis="burden", columns=("n_case_var",)),))

    frame = load_cohort_frame(m)

    assert set(frame["gene"]) == {"A", "B"}
    assert set(frame.columns) == {"gene", "minp", "n_case_var"}


def test_load_prior_frame_matches_priorspec_load_column_names(tmp_path):
    """Must match hvantk.algorithms.rerank.config.PriorSpec.load()'s shape exactly --
    Task 2 swaps a CohortManifest in for a PriorSpec and nothing downstream may
    notice the difference."""
    p = tmp_path / "cohort.tsv"
    _write_tsv(p, ["gene", "minp"], [("A", 0.01), ("B", 0.20)])
    m = _manifest(p)

    frame = load_prior_frame(m)

    assert list(frame.columns) == ["unit", "prior_stat"]
    assert set(frame["unit"]) == {"A", "B"}
    assert set(frame["prior_stat"]) == {0.01, 0.20}
