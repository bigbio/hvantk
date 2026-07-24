"""Attaching a cohort onto the Layer-1 matrix.

The collision check is pure Python (schema introspection only) and is deliberately
unmarked so it runs in the fast suite; the join tests need Hail.
"""
import pytest

from hvantk.algorithms.cohort.attach import (
    as_source_entry,
    attach,
    check_no_layer1_collisions,
)
from hvantk.algorithms.cohort.spec import (
    CohortAxis,
    CohortManifest,
    CohortPrior,
)


def _manifest(axes=(), prior_col="minp", key="gene_id"):
    return CohortManifest(
        name="demo",
        key=key,
        table="/data/demo.tsv",
        prior=CohortPrior(column=prior_col, direction="lower_is_better"),
        cohort_axes=axes,
    )


def _layer1():
    import hail as hl

    return hl.Table.parallelize(
        [
            {"gene_id": "ENSG1", "gene_name": "A", "mis_z": 1.0},
            {"gene_id": "ENSG2", "gene_name": "B", "mis_z": 2.0},
            {"gene_id": "ENSG3", "gene_name": "C", "mis_z": 3.0},
        ],
        hl.tstruct(gene_id=hl.tstr, gene_name=hl.tstr, mis_z=hl.tfloat64),
        key=["gene_id"],
    )


def _cohort():
    """ENSG3 is deliberately absent -- it was never tested by this cohort."""
    import hail as hl

    return hl.Table.parallelize(
        [
            {"gene_id": "ENSG1", "minp": 0.01, "n_case_var": 5},
            {"gene_id": "ENSG2", "minp": 0.20, "n_case_var": 3},
        ],
        hl.tstruct(gene_id=hl.tstr, minp=hl.tfloat64, n_case_var=hl.tint32),
        key=["gene_id"],
    )


def test_as_source_entry_carries_every_declared_column():
    m = _manifest(axes=(CohortAxis(axis="burden", columns=("n_case_var",)),))
    entry = as_source_entry(m)
    assert entry.key == "gene_id"
    assert entry.columns == ("minp", "n_case_var")
    assert entry.min_mapping_rate == m.min_mapping_rate


def test_collision_with_a_layer1_column_fails_loud():
    """Pure-Python schema check -- must run without Hail (no @pytest.mark.hail)."""

    class _FakeLayer1:
        row = {"gene_id": None, "gene_name": None, "mis_z": None}

    m = _manifest(axes=(CohortAxis(axis="constraint", columns=("mis_z",)),))
    with pytest.raises(ValueError, match="mis_z"):
        check_no_layer1_collisions(_FakeLayer1(), m)


def test_no_collision_passes():
    class _FakeLayer1:
        row = {"gene_id": None, "gene_name": None, "mis_z": None}

    check_no_layer1_collisions(_FakeLayer1(), _manifest())


@pytest.mark.hail
def test_attach_keeps_every_spine_gene_and_flags_tested(hail_session):
    m = _manifest(axes=(CohortAxis(axis="burden", columns=("n_case_var",)),))
    attached, report = attach(_layer1(), _cohort(), m)

    assert list(attached.key) == ["gene_id"]
    assert attached.count() == 3
    rows = {r.gene_id: r for r in attached.collect()}
    assert rows["ENSG1"].cohort_tested is True
    assert rows["ENSG2"].cohort_tested is True
    assert rows["ENSG3"].cohort_tested is False
    assert report["n_genes"] == 3
    assert report["n_tested"] == 2


@pytest.mark.hail
def test_attach_never_imputes(hail_session):
    """A gene absent from the cohort table stays missing -- not zero."""
    m = _manifest(axes=(CohortAxis(axis="burden", columns=("n_case_var",)),))
    attached, _ = attach(_layer1(), _cohort(), m)

    rows = {r.gene_id: r for r in attached.collect()}
    assert rows["ENSG3"].minp is None
    assert rows["ENSG3"].n_case_var is None
    assert rows["ENSG1"].minp == pytest.approx(0.01)


@pytest.mark.hail
def test_attach_preserves_layer1_columns(hail_session):
    attached, _ = attach(_layer1(), _cohort(), _manifest())
    assert set(attached.row) == {
        "gene_id",
        "gene_name",
        "mis_z",
        "minp",
        "cohort_tested",
    }


@pytest.mark.hail
def test_report_records_prior_direction_and_per_column_rates(hail_session):
    m = _manifest(axes=(CohortAxis(axis="burden", columns=("n_case_var",)),))
    _, report = attach(_layer1(), _cohort(), m)

    assert report["prior"]["column"] == "minp"
    assert report["prior"]["direction"] == "lower_is_better"
    assert report["axes"]["burden"]["n_case_var"]["non_null_rate"] == pytest.approx(
        2 / 3
    )


class _FakeHGNC:
    """Minimal HGNC stand-in: resolves symbols straight through to gene ids.

    Mirrors the fake-collaborator convention in hvantk/tests/annotation/test_mapping.py --
    tests never build a real HGNC table.
    """

    _SYMBOL_TO_HGNC = {"A": "HGNC:1", "B": "HGNC:2", "GHOST": "HGNC:9"}
    _HGNC_TO_ENSEMBL = {"HGNC:1": "ENSG1", "HGNC:2": "ENSG2", "HGNC:9": "ENSG99"}

    def resolve_to_canonical(self, symbol):
        return symbol

    def map_to_hgnc(self, ids, source_type):
        return {i: self._SYMBOL_TO_HGNC.get(i) for i in ids}

    def map_from_hgnc(self, hgnc_ids, target_type):
        return {h: self._HGNC_TO_ENSEMBL.get(h) for h in hgnc_ids}


def _labelled_manifest(tmp_path, genes, min_rate=0.5):
    import json

    from hvantk.algorithms.cohort.spec import CohortLabels

    payload = {
        "gene_sets": {
            "panel": {
                "name": "panel",
                "genes": sorted(genes),
                "source": "prepare-geneset",
                "n_genes": len(genes),
                "metadata": {},
            }
        },
        "background_genes": sorted(genes),
        "n_gene_sets": 1,
        "n_background": len(genes),
        "source_description": "test",
        "metadata": {},
    }
    p = tmp_path / "panel.json"
    p.write_text(json.dumps(payload))

    m = _manifest()
    return CohortManifest(
        name=m.name,
        key=m.key,
        table=m.table,
        prior=m.prior,
        labels=CohortLabels(gene_set=str(p), min_mapping_rate=min_rate),
    )


@pytest.mark.hail
def test_attach_emits_a_label_column_from_the_declared_gene_set(hail_session, tmp_path):
    manifest = _labelled_manifest(tmp_path, {"A"})
    attached, report = attach(_layer1(), _cohort(), manifest, hgnc=_FakeHGNC())

    rows = {r.gene_id: r for r in attached.collect()}
    assert rows["ENSG1"].label is True
    # False, not missing: a label set is a closed-world assertion.
    assert rows["ENSG2"].label is False
    assert rows["ENSG3"].label is False
    assert report["labels"]["n_mapped"] == 1


@pytest.mark.hail
def test_attach_fails_loud_when_too_few_label_symbols_resolve(hail_session, tmp_path):
    """GHOST maps to ENSG99, which is not on the spine -- so only 1 of 2 resolves."""
    from hvantk.algorithms.annotation.mapping import MappingRateError

    manifest = _labelled_manifest(tmp_path, {"A", "GHOST"}, min_rate=0.9)
    with pytest.raises(MappingRateError):
        attach(_layer1(), _cohort(), manifest, hgnc=_FakeHGNC())
