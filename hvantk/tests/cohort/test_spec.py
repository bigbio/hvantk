"""Cohort manifest parsing and validation.

All pure Python -- no Hail -- so these run in the default fast suite. The contract's
validation logic is deliberately Hail-free (design G1).
"""
from pathlib import Path

import pytest

from hvantk.algorithms.cohort.spec import (
    CohortManifest,
    load_cohort,
)

MINIMAL = """\
name: demo
key: symbol
table: /data/demo_genes.tsv
prior:
  column: minp
  direction: lower_is_better
"""

FULL = """\
name: demo
key: gene_id
table: /data/demo_genes.tsv
min_mapping_rate: 0.95
prior:
  column: p_burden
  direction: lower_is_better
labels:
  gene_set: /data/panel.json
  set_name: cardiac_panel
  min_mapping_rate: 0.8
cohort_axes:
  - axis: burden
    columns: [n_case_var, conc]
  - axis: varagg
    columns: [revel_mean]
"""


def _write(tmp_path: Path, text: str) -> str:
    p = tmp_path / "cohort.yaml"
    p.write_text(text)
    return str(p)


def test_load_cohort_reads_the_minimal_manifest(tmp_path):
    m = load_cohort(_write(tmp_path, MINIMAL))

    assert isinstance(m, CohortManifest)
    assert m.name == "demo"
    assert m.key == "symbol"
    assert m.table == "/data/demo_genes.tsv"
    assert m.prior.column == "minp"
    assert m.prior.direction == "lower_is_better"
    # Optional blocks default to absent/empty -- a cohort is key + prior (D2).
    assert m.labels is None
    assert m.cohort_axes == ()
    assert m.min_mapping_rate == 0.9


def test_load_cohort_reads_every_optional_block(tmp_path):
    m = load_cohort(_write(tmp_path, FULL))

    assert m.min_mapping_rate == 0.95
    assert m.labels.gene_set == "/data/panel.json"
    assert m.labels.set_name == "cardiac_panel"
    assert m.labels.min_mapping_rate == 0.8
    assert [a.axis for a in m.cohort_axes] == ["burden", "varagg"]
    assert m.axis("burden").columns == ("n_case_var", "conc")


def test_declared_columns_is_prior_plus_every_axis_column(tmp_path):
    m = load_cohort(_write(tmp_path, FULL))
    assert m.declared_columns() == ("p_burden", "n_case_var", "conc", "revel_mean")


def test_axis_raises_keyerror_for_an_unknown_axis(tmp_path):
    m = load_cohort(_write(tmp_path, FULL))
    with pytest.raises(KeyError, match="no cohort axis 'nope'"):
        m.axis("nope")


def test_prior_is_required(tmp_path):
    text = "name: demo\nkey: symbol\ntable: /data/t.tsv\n"
    with pytest.raises(Exception):
        load_cohort(_write(tmp_path, text))


def test_prior_direction_is_required(tmp_path):
    text = "name: demo\nkey: symbol\ntable: /data/t.tsv\n" "prior:\n  column: minp\n"
    with pytest.raises(Exception):
        load_cohort(_write(tmp_path, text))


def test_prior_direction_must_be_one_of_the_two_values(tmp_path):
    text = (
        "name: demo\nkey: symbol\ntable: /data/t.tsv\n"
        "prior:\n  column: minp\n  direction: ascending\n"
    )
    with pytest.raises(Exception):
        load_cohort(_write(tmp_path, text))


def test_key_must_be_a_spine_mappable_identifier(tmp_path):
    text = (
        "name: demo\nkey: protein_id\ntable: /data/t.tsv\n"
        "prior:\n  column: minp\n  direction: lower_is_better\n"
    )
    with pytest.raises(Exception):
        load_cohort(_write(tmp_path, text))


def test_unknown_top_level_key_is_rejected(tmp_path):
    text = MINIMAL + "universe: all\n"
    with pytest.raises(Exception):
        load_cohort(_write(tmp_path, text))


def test_duplicate_axis_label_is_rejected(tmp_path):
    text = MINIMAL + (
        "cohort_axes:\n"
        "  - {axis: burden, columns: [a]}\n"
        "  - {axis: burden, columns: [b]}\n"
    )
    with pytest.raises(ValueError, match="duplicate cohort axis"):
        load_cohort(_write(tmp_path, text))


def test_same_column_declared_by_two_axes_is_rejected(tmp_path):
    text = MINIMAL + (
        "cohort_axes:\n"
        "  - {axis: burden, columns: [shared]}\n"
        "  - {axis: varagg, columns: [shared]}\n"
    )
    with pytest.raises(ValueError, match="duplicate declared column 'shared'"):
        load_cohort(_write(tmp_path, text))


def test_axis_column_colliding_with_the_prior_column_is_rejected(tmp_path):
    """The minp trap: the prior stat and a same-named model feature are different
    transforms of the same quantity (raw p vs -log10 p). Declaring both silently
    feeds the untransformed value in as a feature."""
    text = MINIMAL + ("cohort_axes:\n" "  - {axis: burden, columns: [minp]}\n")
    with pytest.raises(ValueError, match="duplicate declared column 'minp'"):
        load_cohort(_write(tmp_path, text))
