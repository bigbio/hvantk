"""Unit tests for `combine_gvcfs` partitioning handling.

These exercise the algorithm directly (not through the CLI), because the CLI does not
expose ``intervals`` and therefore cannot reach some of these branches at all.
Hail is stubbed out: every assertion here is about what `combine_gvcfs` decides *before*
it hands off to ``hl.vds.new_combiner``.
"""

from unittest.mock import MagicMock, patch

import pytest

from hvantk.algorithms.hgc.combiners import combine_gvcfs


def _run(kwargs, gvcfs=("/gvcfs/a.g.vcf.gz",)):
    """Call combine_gvcfs with Hail + path validation stubbed; return the new_combiner mock."""
    with patch(
        "hvantk.algorithms.hgc.combiners.validate_vcfs_paths", return_value=list(gvcfs)
    ):
        with patch("hvantk.algorithms.hgc.combiners.hl") as mock_hl:
            mock_hl.vds.new_combiner.return_value = MagicMock()
            combine_gvcfs(
                gvcf_dir="/gvcfs",
                vds_output_path="/out.vds",
                tmp_path="/tmp",
                save_path=None,
                vdses=[],
                kwargs=kwargs,
            )
            return mock_hl.vds.new_combiner


def test_no_partitioning_option_defaults_to_genome_intervals():
    """With no options, Hail's genome default is applied (unchanged legacy behaviour)."""
    nc = _run({})
    call = nc.call_args[1]
    assert call["use_genome_default_intervals"] is True
    assert call["import_interval_size"] is None
    assert call["intervals"] is None


def test_import_interval_size_is_forwarded_and_suppresses_genome_default():
    """import_interval_size must NOT also trigger use_genome_default_intervals.

    Hail only *warns* about colliding partition arguments and silently picks one, so a
    regression here would quietly ignore the user's requested interval size.
    """
    nc = _run({"import_interval_size": 600_000})
    call = nc.call_args[1]
    assert call["import_interval_size"] == 600_000
    assert call["use_genome_default_intervals"] is False
    assert call["use_exome_default_intervals"] is False


def test_explicit_intervals_suppress_genome_default():
    nc = _run({"intervals": ["iv1", "iv2"]})
    call = nc.call_args[1]
    assert call["intervals"] == ["iv1", "iv2"]
    assert call["use_genome_default_intervals"] is False


def test_intervals_may_be_any_iterable():
    """A generator must not blow up (len() is used for logging)."""
    nc = _run({"intervals": (f"iv{i}" for i in range(3))})
    assert nc.call_args[1]["intervals"] == ["iv0", "iv1", "iv2"]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"import_interval_size": 600_000, "use_genome_default_intervals": True},
        {"import_interval_size": 600_000, "use_exome_default_intervals": True},
        {"intervals": ["iv"], "import_interval_size": 600_000},
        {"use_genome_default_intervals": True, "use_exome_default_intervals": True},
    ],
)
def test_colliding_partitioning_options_raise(kwargs):
    """Colliding options must fail loudly rather than let Hail silently pick one."""
    with pytest.raises(ValueError, match="Only one of"):
        _run(kwargs)


@pytest.mark.parametrize("bad", [0, -1])
def test_non_positive_import_interval_size_raises(bad):
    """0 would reach Hail as ceil(contig_length / 0) -> ZeroDivisionError."""
    with pytest.raises(ValueError, match="at least 1"):
        _run({"import_interval_size": bad})


def test_caller_kwargs_dict_is_not_mutated():
    """combine_gvcfs must not consume the caller's dict.

    Regression guard: the interval keys are read with .pop(). If they were popped from
    the caller's own dict, a natural pattern like

        opts = {"import_interval_size": 600_000}
        for contig in contigs:
            combine_gvcfs(..., kwargs=opts)

    would silently lose the setting after the first iteration and fall back to the
    genome default -- exactly the parallelism ceiling these options exist to avoid.
    """
    opts = {"import_interval_size": 600_000, "gvcf_batch_size": 100}
    nc = _run(opts)

    assert opts == {"import_interval_size": 600_000, "gvcf_batch_size": 100}
    # ...and the setting really did reach Hail on this call
    assert nc.call_args[1]["import_interval_size"] == 600_000


def test_tree_merge_options_are_forwarded_unchanged():
    nc = _run({"gvcf_batch_size": 100, "branch_factor": 50})
    call = nc.call_args[1]
    assert call["gvcf_batch_size"] == 100
    assert call["branch_factor"] == 50
