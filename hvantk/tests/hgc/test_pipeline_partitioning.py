"""The HGC pipeline must be able to tune gVCF-combiner partitioning.

`hvantk hgc pipeline` is the documented, recommended entry point, but it previously
passed ``kwargs={}`` to ``combine_gvcfs`` unconditionally -- so a pipeline user was stuck
with Hail's 1.2 Mb genome default and had no way to raise the partition count above their
core count. These tests pin the wiring so that cannot silently regress.
"""

import pytest

from hvantk.algorithms.hgc.pipeline import PipelineConfig


def _cfg(tmp_path, **overrides):
    return PipelineConfig(
        input_dir=str(tmp_path), output_dir=str(tmp_path / "out"), **overrides
    )


def test_combiner_kwargs_empty_by_default(tmp_path):
    """No tuning options -> empty kwargs -> Hail's genome default (legacy behaviour)."""
    assert _cfg(tmp_path).combiner_kwargs() == {}


def test_combiner_kwargs_carries_import_interval_size(tmp_path):
    cfg = _cfg(tmp_path, import_interval_size=600_000)
    assert cfg.combiner_kwargs() == {"import_interval_size": 600_000}


def test_combiner_kwargs_carries_tree_merge_options(tmp_path):
    cfg = _cfg(tmp_path, gvcf_batch_size=100, branch_factor=50)
    assert cfg.combiner_kwargs() == {"gvcf_batch_size": 100, "branch_factor": 50}


def test_combiner_kwargs_carries_exome_default(tmp_path):
    cfg = _cfg(tmp_path, use_exome_default_intervals=True)
    assert cfg.combiner_kwargs() == {"use_exome_default_intervals": True}


@pytest.mark.parametrize(
    "overrides, expected",
    [
        ({"import_interval_size": 0}, "at least 1 bp"),
        ({"import_interval_size": -5}, "at least 1 bp"),
        ({"gvcf_batch_size": 0}, "at least 1"),
        ({"branch_factor": 1}, "at least 2"),
        (
            {"import_interval_size": 600_000, "use_exome_default_intervals": True},
            "mutually exclusive",
        ),
    ],
)
def test_invalid_combiner_tuning_is_rejected_by_validate(tmp_path, overrides, expected):
    """Bad values must be caught by config validation, not deep inside Hail.

    Hail does reject these, but only after init + gVCF validation, and the error is then
    wrapped in generic "check your Spark version / GVCF files" guidance that points the
    user at the wrong cause entirely.
    """
    errors = _cfg(tmp_path, **overrides).validate()
    assert any(expected in e for e in errors), errors


def test_valid_combiner_tuning_passes_validate(tmp_path):
    errors = _cfg(
        tmp_path, import_interval_size=600_000, gvcf_batch_size=50, branch_factor=100
    ).validate()
    combiner_errors = [
        e for e in errors if "interval" in e or "batch" in e or "branch" in e
    ]
    assert combiner_errors == []
