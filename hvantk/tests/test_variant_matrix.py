"""Tests for VariantMatrix: multi-sample variant cohort backed by a Hail MatrixTable.

All tests here require a live Hail session and are explicitly marked
``@pytest.mark.hail``.
"""
from __future__ import annotations

from datetime import datetime, timezone

import pytest

from hvantk.core.models._expr import col
from hvantk.core.models.variant_matrix import VariantMatrix
from hvantk.core.models.provenance import Provenance


def _prov():
    return Provenance(
        plugin="t",
        dataset="t:mt",
        plugin_version="0.0",
        source_fingerprint="sha256:abc",
        schema_id="t-mt-v1",
        build_timestamp=datetime(2026, 5, 21, tzinfo=timezone.utc),
        builder_commit=None,
    )


@pytest.fixture
def small_mt(hail_session):
    """A 3×2 Hail MatrixTable (3 variants / rows × 2 samples / cols).

    Entry value layout (row_idx × col_idx):
        col 0  col 1
    row 0:  0.0    1.0
    row 1: 10.0   11.0
    row 2: 20.0   21.0

    After VariantMatrix wrapping (n_samples=2, n_variants=3).
    """
    import hail as hl
    rows = []
    for r in range(3):   # rows = variants
        for c in range(2):   # cols = samples
            rows.append({"row_idx": r, "col_idx": c, "value": float(r * 10 + c)})
    ht = hl.Table.parallelize(
        rows,
        hl.tstruct(row_idx=hl.tint32, col_idx=hl.tint32, value=hl.tfloat64),
    )
    return ht.to_matrix_table(
        row_key=["row_idx"],
        col_key=["col_idx"],
    )


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

@pytest.mark.hail
def test_from_hail_mt_attributes(small_mt):
    vm = VariantMatrix.from_hail_mt(small_mt, provenance=_prov())
    assert vm.n_samples == 2
    assert vm.n_variants == 3
    assert vm.provenance == _prov()


@pytest.mark.hail
def test_from_hail_mt_stores_matrix(small_mt):
    import hail as hl
    vm = VariantMatrix.from_hail_mt(small_mt, provenance=_prov())
    assert isinstance(vm._mt, hl.MatrixTable)


# ---------------------------------------------------------------------------
# samples / variants accessors
# ---------------------------------------------------------------------------

@pytest.mark.hail
def test_samples_is_annotation_table_hail(small_mt):
    from hvantk.core.models.annotation_table import AnnotationTable
    vm = VariantMatrix.from_hail_mt(small_mt, provenance=_prov())
    samples = vm.samples
    assert isinstance(samples, AnnotationTable)
    assert samples.backend == "hail"


@pytest.mark.hail
def test_variants_is_annotation_table_hail(small_mt):
    from hvantk.core.models.annotation_table import AnnotationTable
    vm = VariantMatrix.from_hail_mt(small_mt, provenance=_prov())
    variants = vm.variants
    assert isinstance(variants, AnnotationTable)
    assert variants.backend == "hail"


@pytest.mark.hail
def test_samples_count_matches_n_samples(small_mt):
    vm = VariantMatrix.from_hail_mt(small_mt, provenance=_prov())
    assert vm.samples.count() == vm.n_samples


@pytest.mark.hail
def test_variants_count_matches_n_variants(small_mt):
    vm = VariantMatrix.from_hail_mt(small_mt, provenance=_prov())
    assert vm.variants.count() == vm.n_variants


# ---------------------------------------------------------------------------
# to_hail_mt (escape hatch)
# ---------------------------------------------------------------------------

@pytest.mark.hail
def test_to_hail_mt_returns_matrix_table(small_mt):
    """to_hail_mt() returns the underlying Hail MatrixTable."""
    import hail as hl
    vm = VariantMatrix.from_hail_mt(small_mt, provenance=_prov())
    mt = vm.to_hail_mt()
    assert isinstance(mt, hl.MatrixTable)
    n_rows, n_cols = mt.count()
    assert n_rows == 3
    assert n_cols == 2


# ---------------------------------------------------------------------------
# subset_samples / subset_variants
# ---------------------------------------------------------------------------

@pytest.mark.hail
def test_subset_samples_hail_mt(small_mt):
    """Filter cols (samples) by col_idx == 0 should leave 1 sample."""
    vm = VariantMatrix.from_hail_mt(small_mt, provenance=_prov())
    sub = vm.subset_samples(col("col_idx") == 0)
    assert sub.n_samples == 1
    assert sub.n_variants == 3


@pytest.mark.hail
def test_subset_variants_hail_mt(small_mt):
    """Filter rows (variants) by row_idx == 0 should leave 1 variant."""
    vm = VariantMatrix.from_hail_mt(small_mt, provenance=_prov())
    sub = vm.subset_variants(col("row_idx") == 0)
    assert sub.n_samples == 2
    assert sub.n_variants == 1


@pytest.mark.hail
def test_subset_samples_gt_predicate(small_mt):
    """col_idx > 0 selects samples with col_idx=1 only."""
    vm = VariantMatrix.from_hail_mt(small_mt, provenance=_prov())
    sub = vm.subset_samples(col("col_idx") > 0)
    assert sub.n_samples == 1


@pytest.mark.hail
def test_subset_variants_gt_predicate(small_mt):
    """row_idx > 0 selects variants with row_idx=1 and row_idx=2."""
    vm = VariantMatrix.from_hail_mt(small_mt, provenance=_prov())
    sub = vm.subset_variants(col("row_idx") > 0)
    assert sub.n_variants == 2


# ---------------------------------------------------------------------------
# core/io round-trip: save → load
# ---------------------------------------------------------------------------

@pytest.mark.hail
def test_save_load_mt_round_trip(tmp_path, small_mt):
    """Save a VariantMatrix as .mt/ via core/io; load it back."""
    from hvantk.core import io as core_io

    vm = VariantMatrix.from_hail_mt(small_mt, provenance=_prov())
    out = tmp_path / "variants.mt"
    core_io.save(vm, out)

    # The .mt directory must exist and a sidecar provenance file alongside it
    assert out.exists()
    assert (out.with_name("variants.mt.provenance.json")).exists()

    loaded = core_io.load(out)
    assert isinstance(loaded, VariantMatrix)
    assert loaded.provenance == vm.provenance
    assert loaded.n_samples == 2
    assert loaded.n_variants == 3


@pytest.mark.hail
def test_save_load_mt_entry_count_parity(tmp_path, small_mt):
    """Entry counts survive a save/load round-trip."""
    from hvantk.core import io as core_io

    vm = VariantMatrix.from_hail_mt(small_mt, provenance=_prov())
    out = tmp_path / "variants2.mt"
    core_io.save(vm, out)

    loaded = core_io.load(out)
    assert loaded.n_variants == vm.n_variants
    assert loaded.n_samples == vm.n_samples
    assert loaded.provenance.schema_id == vm.provenance.schema_id


# ---------------------------------------------------------------------------
# Lazy count behavior
# ---------------------------------------------------------------------------

@pytest.mark.hail
def test_from_hail_mt_does_not_eagerly_count(monkeypatch, hail_session):
    """from_hail_mt should not trigger count_cols / count_rows at construction."""
    import hail as hl

    rows = []
    for r in range(3):
        for c in range(2):
            rows.append({"row_idx": r, "col_idx": c, "value": float(r * 10 + c)})
    ht = hl.Table.parallelize(
        rows,
        hl.tstruct(row_idx=hl.tint32, col_idx=hl.tint32, value=hl.tfloat64),
    )
    mt = ht.to_matrix_table(row_key=["row_idx"], col_key=["col_idx"])

    # Track count calls via monkeypatching the type
    counts = {"cols": 0, "rows": 0}
    orig_count_cols = type(mt).count_cols
    orig_count_rows = type(mt).count_rows

    def spy_count_cols(self):
        counts["cols"] += 1
        return orig_count_cols(self)

    def spy_count_rows(self):
        counts["rows"] += 1
        return orig_count_rows(self)

    monkeypatch.setattr(type(mt), "count_cols", spy_count_cols)
    monkeypatch.setattr(type(mt), "count_rows", spy_count_rows)

    vm = VariantMatrix.from_hail_mt(mt, provenance=_prov())

    # Construction must NOT have triggered any counts
    assert counts == {"cols": 0, "rows": 0}, (
        f"from_hail_mt triggered eager counts: {counts}"
    )

    # Accessing n_samples triggers exactly one count_cols call, then caches
    _ = vm.n_samples
    _ = vm.n_samples  # second access — should NOT re-count
    assert counts["cols"] == 1, "n_samples should trigger count_cols exactly once"
    assert counts["rows"] == 0, "n_samples must not trigger count_rows"

    # Accessing n_variants triggers exactly one count_rows call, then caches
    _ = vm.n_variants
    _ = vm.n_variants  # second access — should NOT re-count
    assert counts["rows"] == 1, "n_variants should trigger count_rows exactly once"
    assert counts["cols"] == 1, "n_variants must not re-trigger count_cols"
