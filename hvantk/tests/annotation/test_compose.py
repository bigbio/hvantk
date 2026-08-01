"""compose: generic Stage-2 left-join of prepared Layer-1 axes onto the spine."""
from __future__ import annotations

import pytest

from hvantk.algorithms.annotation.spec import FeatureSpec, SourceEntry


def _spine():
    import hail as hl

    return hl.Table.parallelize(
        [
            {"gene_id": "ENSG1", "gene_name": "A"},
            {"gene_id": "ENSG2", "gene_name": "B"},
            {"gene_id": "ENSG3", "gene_name": "C"},
        ],
        hl.tstruct(gene_id=hl.tstr, gene_name=hl.tstr),
        key=["gene_id"],
    )


def _axis(rows, schema):
    import hail as hl

    return hl.Table.parallelize(rows, schema, key=["gene_id"])


CONSTRAINT_ENTRY = SourceEntry(
    axis="constraint",
    source="gnomad-metrics:metrics",
    key="gene_id",
    columns=("mis_z",),
)
GEVIR_ENTRY = SourceEntry(
    axis="gevir",
    source="gevir:table",
    key="gene_id",
    columns=("gevir_pct",),
)


def _two_axis_spec():
    return FeatureSpec(name="test-spec", layer1=(CONSTRAINT_ENTRY, GEVIR_ENTRY))


def _constraint_axis_ht():
    import hail as hl

    # ENSG1/ENSG2 only -- ENSG3 absent.
    return _axis(
        [{"gene_id": "ENSG1", "mis_z": 2.0}, {"gene_id": "ENSG2", "mis_z": -1.0}],
        hl.tstruct(gene_id=hl.tstr, mis_z=hl.tfloat64),
    )


def _gevir_axis_ht():
    import hail as hl

    # ENSG2/ENSG3 only -- ENSG1 absent.
    return _axis(
        [
            {"gene_id": "ENSG2", "gevir_pct": 0.5},
            {"gene_id": "ENSG3", "gevir_pct": 0.9},
        ],
        hl.tstruct(gene_id=hl.tstr, gevir_pct=hl.tfloat64),
    )


@pytest.mark.hail
def test_compose_left_joins_all_axes_onto_the_spine(hail_session):
    from hvantk.algorithms.annotation.compose import compose

    prepared_by_axis = {
        "constraint": _constraint_axis_ht(),
        "gevir": _gevir_axis_ht(),
    }
    composed, manifest = compose(_spine(), prepared_by_axis, _two_axis_spec())

    assert list(composed.key) == ["gene_id"]
    assert composed.count() == 3
    assert set(composed.row) == {
        "gene_id",
        "gene_name",
        "mis_z",
        "gevir_pct",
        "constraint_present",
        "gevir_present",
    }


@pytest.mark.hail
def test_compose_never_imputes_missing_stays_missing(hail_session):
    from hvantk.algorithms.annotation.compose import compose

    prepared_by_axis = {
        "constraint": _constraint_axis_ht(),
        "gevir": _gevir_axis_ht(),
    }
    composed, _manifest = compose(_spine(), prepared_by_axis, _two_axis_spec())
    rows = {r.gene_id: r for r in composed.collect()}

    # ENSG3 absent from constraint -> mis_z missing (not 0), presence flag False.
    assert rows["ENSG3"].mis_z is None
    assert rows["ENSG3"].constraint_present is False

    # ENSG1 absent from gevir -> gevir_pct missing, presence flag False.
    assert rows["ENSG1"].gevir_pct is None
    assert rows["ENSG1"].gevir_present is False


@pytest.mark.hail
def test_compose_presence_flag_true_when_axis_has_the_gene(hail_session):
    from hvantk.algorithms.annotation.compose import compose

    prepared_by_axis = {
        "constraint": _constraint_axis_ht(),
        "gevir": _gevir_axis_ht(),
    }
    composed, _manifest = compose(_spine(), prepared_by_axis, _two_axis_spec())
    rows = {r.gene_id: r for r in composed.collect()}

    # ENSG2 is present in BOTH axes.
    row = rows["ENSG2"]
    assert row.constraint_present is True
    assert row.gevir_present is True
    assert row.mis_z == pytest.approx(-1.0)
    assert row.gevir_pct == pytest.approx(0.5)


@pytest.mark.hail
def test_compose_manifest_reports_nonnull_and_positive_rate(hail_session):
    from hvantk.algorithms.annotation.compose import compose

    # Single-axis spec so the manifest arithmetic is easy to hand-check:
    # mis_z: ENSG1=2.0, ENSG2=-1.0, ENSG3 missing -> non-null 2/3, positive 1/2.
    spec = FeatureSpec(name="test-spec", layer1=(CONSTRAINT_ENTRY,))
    composed, manifest = compose(_spine(), {"constraint": _constraint_axis_ht()}, spec)

    stats = manifest["axes"]["constraint"]["mis_z"]
    assert stats["non_null_rate"] == pytest.approx(2 / 3)
    assert stats["positive_rate"] == pytest.approx(1 / 2)
    assert manifest["n_genes"] == 3


def test_compose_raises_on_duplicate_output_column_across_axes():
    # Pure-Python collision check -- must run without Hail (no @pytest.mark.hail).
    from hvantk.algorithms.annotation.compose import compose

    dup_entry = SourceEntry(
        axis="gevir",
        source="gevir:table",
        key="gene_id",
        columns=("mis_z",),  # collides with CONSTRAINT_ENTRY's "mis_z"
    )
    spec = FeatureSpec(name="test-spec", layer1=(CONSTRAINT_ENTRY, dup_entry))

    with pytest.raises(ValueError, match="mis_z"):
        compose(spine=None, prepared_by_axis={}, spec=spec)


@pytest.mark.hail
def test_compose_carries_matrix_vector_columns_not_just_declared(hail_session):
    """A matrix axis contributes every column its prepared table holds.

    Regression: compose annotated only ``entry.columns``, so a matrix source's per-group
    specificity vector -- whose names come from the atlas's group labels at runtime and
    so cannot be enumerated in the spec -- was silently dropped. That made
    ``emit="vector"``, the documented default, a no-op end to end.
    """
    import hail as hl

    from hvantk.algorithms.annotation.compose import compose

    entry = SourceEntry(
        axis="expr",
        source="ucsc-cellbrowser:asp_2019",
        key="symbol",
        columns=("asp_cm_spec",),  # the spec names ONLY the roll-up
    )
    # ...while the prepared table also carries the per-group vector.
    prepared = _axis(
        [
            {"gene_id": "ENSG1", "asp_cm_spec": 1.0, "asp_cm": 1.0, "asp_other": 0.0},
            {"gene_id": "ENSG2", "asp_cm_spec": 0.5, "asp_cm": 0.5, "asp_other": 0.5},
        ],
        hl.tstruct(
            gene_id=hl.tstr,
            asp_cm_spec=hl.tfloat64,
            asp_cm=hl.tfloat64,
            asp_other=hl.tfloat64,
        ),
    )
    spec = FeatureSpec(name="test-spec", layer1=(entry,))
    composed, manifest = compose(_spine(), {"expr": prepared}, spec)

    assert {"asp_cm_spec", "asp_cm", "asp_other"} <= set(composed.row)
    assert set(manifest["axes"]["expr"]) == {"asp_cm_spec", "asp_cm", "asp_other"}
    rows = {r.gene_id: r for r in composed.collect()}
    assert rows["ENSG1"].asp_other == pytest.approx(0.0)
    assert rows["ENSG2"].asp_other == pytest.approx(0.5)
    assert rows["ENSG3"].asp_cm is None  # absent from the axis, left missing
