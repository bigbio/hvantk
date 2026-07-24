"""The pandas path onto a cohort manifest -- the ``attach()`` sibling for ``rerank``.

``rerank`` is pandas + scikit-learn with zero Hail, so a ``CohortManifest`` reaches it
through :mod:`hvantk.algorithms.cohort.frame` rather than through
:func:`hvantk.algorithms.cohort.attach.attach`. Pure Python -- no Hail -- so these run
in the default fast suite.
"""
import gzip

import pandas as pd
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


def test_load_cohort_frame_excludes_the_prior_column_when_asked(tmp_path):
    # Findings 1+2 (whole-branch review): the audit merge/collision-check must be able
    # to see the cohort's axis columns without the prior column tagging along -- the
    # prior is already consumed elsewhere (as prior_stat) by the time that matters.
    p = tmp_path / "cohort.tsv"
    _write_tsv(p, ["gene", "minp", "n_case_var"], [("A", 0.01, 5), ("B", 0.20, 3)])
    m = _manifest(p, axes=(CohortAxis(axis="burden", columns=("n_case_var",)),))

    with_prior = load_cohort_frame(m)
    without_prior = load_cohort_frame(m, include_prior=False)

    assert set(with_prior.columns) == {"gene", "minp", "n_case_var"}
    assert set(without_prior.columns) == {"gene", "n_case_var"}


def test_load_cohort_frame_agrees_with_validate_on_a_comma_header_with_a_space(
    tmp_path,
):
    # Finding 6 (whole-branch review): check_declared_columns_exist/read_header
    # strips whitespace from header names ('validate' promises this); load must agree
    # or a table that passes 'hvantk cohort validate' can still blow up on load with a
    # bare KeyError naming neither the manifest nor the real column name.
    p = tmp_path / "cohort.csv"
    p.write_text("gene, minp\nA,0.01\nB,0.20\n")
    m = _manifest(p, prior_col="minp")

    from hvantk.algorithms.cohort.checks import (
        check_declared_columns_exist,
        check_key_column_exists,
        read_header,
    )

    header = read_header(str(p))
    check_key_column_exists(m, header)  # validate() succeeds ...
    check_declared_columns_exist(m, header)

    frame = load_cohort_frame(m)  # ... so load() must succeed too.
    assert set(frame.columns) == {"gene", "minp"}
    assert dict(zip(frame["gene"], frame["minp"])) == {"A": 0.01, "B": 0.20}


def test_load_cohort_frame_raises_on_empty_gene_keys_and_counts_them(tmp_path):
    # Finding 7 (whole-branch review): value_counts() drops NaN by default, so blank
    # gene cells used to sail through the duplicate check, enter the frame as NaN
    # genes, and vanish silently in a downstream left-merge. They must be rejected
    # instead, with a count of how many rows are affected.
    p = tmp_path / "cohort.tsv"
    p.write_text("gene\tminp\nA\t0.01\n\t0.20\nB\t0.3\n   \t0.4\n")
    m = _manifest(p)

    with pytest.raises(ValueError, match="2"):
        load_cohort_frame(m)


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


def test_load_prior_frame_matches_priorspec_load_exactly(tmp_path):
    """Task 1 review (Minor): the test above only checks column names and value
    *sets* and never exercises the real PriorSpec.load() -- two loaders that agree on
    names/sets could still disagree on row order or dtype and this would not catch
    it. Compare full frames instead, so the two implementations cannot drift apart
    silently."""
    from hvantk.algorithms.rerank.config import PriorSpec

    p = tmp_path / "cohort.tsv"
    _write_tsv(p, ["gene", "minp"], [("A", 0.01), ("B", 0.20), ("C", 0.55)])
    m = _manifest(p)

    manifest_frame = load_prior_frame(m)
    priorspec_frame = PriorSpec(path=str(p), unit_col="gene", stat_col="minp").load()

    pd.testing.assert_frame_equal(manifest_frame, priorspec_frame)
