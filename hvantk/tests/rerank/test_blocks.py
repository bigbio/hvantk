"""The blocks Pass B's folds are built from, and the one case where blocking must refuse.

HGNC `gene_group` is multi-membership and pipe-separated. The obvious blocker -- connected
components over shared membership -- is wrong here, and wrong in a way that produces
plausible output: membership chains transitively, so the closure can collapse a large
fraction of a gene universe into one component.
"""
from __future__ import annotations

import numpy as np
import pytest

from hvantk.algorithms.rerank.blocks import (
    DEFAULT_MAX_BLOCK_FRAC,
    BlockPolicy,
    DominantBlockError,
    gene_blocks,
    load_gene_groups,
)
from hvantk.tests.rerank._synth import paralogue_groups


def test_genes_in_one_primary_family_share_a_block():
    genes, mapping = paralogue_groups(n=60, family_size=6)
    rep = gene_blocks(genes, mapping)
    by_family = {}
    for g, b in zip(genes, rep.blocks):
        by_family.setdefault(mapping[g].split("|")[0], set()).add(int(b))
    for fam, blocks in by_family.items():
        if fam:
            assert len(blocks) == 1, (fam, blocks)


def test_an_ungrouped_gene_is_its_own_block():
    """A gene with no known family has no paralogue to leak through, so it must be
    unconstrained rather than pooled with every other ungrouped gene."""
    # max_block_frac=1.0: this fixture is about grouping correctness, not the ceiling, and
    # at n=3 a single gene is already 33% of the universe -- the ceiling is exercised on its
    # own purpose-built (n=100) fixtures below, not on these tiny illustrative ones.
    rep = gene_blocks(["a", "b", "c"], {"a": "", "b": "", "c": "Fam"}, max_block_frac=1.0)
    assert len({int(rep.blocks[0]), int(rep.blocks[1])}) == 2


def test_only_the_primary_group_is_used_not_connected_components():
    """Every gene here shares a SECONDARY group, so a components blocker would return one
    block for the whole matrix. Primary-only keeps them apart."""
    mapping = {"a": "Fam A|Shared", "b": "Fam B|Shared", "c": "Fam A|Other"}
    # max_block_frac=1.0: see test_an_ungrouped_gene_is_its_own_block -- two of three genes
    # sharing a primary group is inherently a 67% block at this scale, which is the point
    # being tested, not a ceiling violation.
    rep = gene_blocks(["a", "b", "c"], mapping, max_block_frac=1.0)
    assert rep.blocks[0] == rep.blocks[2]
    assert rep.blocks[1] != rep.blocks[0]
    assert rep.n_blocks == 2


def test_whitespace_and_missing_values_are_treated_as_ungrouped():
    rep = gene_blocks(["a", "b", "c"], {"a": "  ", "b": None, "c": "Fam"}, max_block_frac=1.0)
    assert rep.n_blocks == 3


def test_a_gene_absent_from_the_mapping_is_a_singleton():
    rep = gene_blocks(["a", "b"], {"a": "Fam"}, max_block_frac=1.0)
    assert rep.n_blocks == 2


def test_the_report_describes_the_blocking():
    genes, mapping = paralogue_groups(n=60, family_size=6, seed=1)
    rep = gene_blocks(genes, mapping)
    assert rep.largest == max(np.bincount(rep.blocks))
    assert rep.largest_frac == pytest.approx(rep.largest / len(genes))
    assert 0 < rep.n_in_multi <= len(genes)
    assert rep.n_blocks == len(set(rep.blocks.tolist()))
    # D1: WHICH blocking, not only whether one was used -- two different groupings of the
    # same genes must not collide, and re-blocking the same genes the same way must
    # reproduce exactly.
    assert gene_blocks(genes, mapping).digest == rep.digest
    _, other_mapping = paralogue_groups(n=60, family_size=6, seed=2)
    assert gene_blocks(genes, other_mapping).digest != rep.digest


# --- the hard abort -------------------------------------------------------------------------


def test_a_dominant_block_aborts_above_the_ceiling():
    """Hard failure, not a warning. StratifiedGroupKFold must put a whole block in one
    fold, so an oversized block forces one fold to hold a large fraction of the data and the
    pooled OOF AUC stops estimating the same quantity as the ungrouped run. It shows as a
    'harder' grouped model scoring HIGHER than random folds, which blocking cannot do."""
    genes, mapping = paralogue_groups(n=100, family_size=4, dominant=40)
    with pytest.raises(DominantBlockError) as info:
        gene_blocks(genes, mapping)
    msg = str(info.value)
    assert "40" in msg and "100" in msg
    assert "10%" in msg or "0.1" in msg
    assert "max_block_frac" in msg


def test_the_ceiling_is_configurable_upwards():
    genes, mapping = paralogue_groups(n=100, family_size=4, dominant=40)
    rep = gene_blocks(genes, mapping, max_block_frac=0.5)
    assert rep.largest == 40


def test_just_under_the_ceiling_is_allowed():
    genes, mapping = paralogue_groups(n=100, family_size=4, dominant=10)
    rep = gene_blocks(genes, mapping, max_block_frac=0.10)
    assert rep.largest == 10  # == the ceiling exactly; the abort is on strictly greater


@pytest.mark.parametrize("bad", [0.0, -0.1, 1.5, float("nan"), True])
def test_an_impossible_ceiling_is_rejected(bad):
    with pytest.raises(ValueError, match="max_block_frac"):
        gene_blocks(["a"], {"a": ""}, max_block_frac=bad)


def test_the_default_ceiling_is_ten_percent():
    assert DEFAULT_MAX_BLOCK_FRAC == 0.10


# --- a table that cannot possibly be blocking this universe ----------------------------------


def test_a_table_that_matches_no_gene_is_an_error():
    with pytest.raises(ValueError, match="none of"):
        gene_blocks(["ENSG1", "ENSG2"], {"A1BG": "Fam"})


# --- the HGNC reader --------------------------------------------------------------------------


def test_load_gene_groups_reads_an_hgnc_shaped_table(tmp_path):
    p = tmp_path / "hgnc.txt"
    p.write_text(
        "symbol\tgene_group\tstatus\n"
        "AAA\tMyosin light chains\tApproved\n"
        "BBB\t\tApproved\n"
        "CCC\tKinases|Other\tApproved\n"
        "DDD\tKinases\tEntry Withdrawn\n"
    )
    m = load_gene_groups(p)
    assert m["AAA"] == "Myosin light chains"
    assert m["CCC"] == "Kinases|Other"
    assert m["BBB"] == ""
    assert "DDD" not in m, "withdrawn entries must not define a family"


def test_load_gene_groups_reports_a_missing_column(tmp_path):
    p = tmp_path / "bad.txt"
    p.write_text("symbol\tstatus\nAAA\tApproved\n")
    with pytest.raises(ValueError, match="gene_group"):
        load_gene_groups(p)


# --- config wiring ----------------------------------------------------------------------------


def test_config_defaults_to_no_blocks():
    from hvantk.algorithms.rerank.config import Config

    assert Config.__dataclass_fields__["blocks"].default is None


def test_config_rejects_a_bare_path_instead_of_a_policy():
    from hvantk.algorithms.rerank.config import Config

    with pytest.raises(TypeError, match="BlockPolicy"):
        Config(name="x", features=[], labels=None, blocks="/some/hgnc.txt").__post_init__()


def test_config_rejects_an_impossible_ceiling_at_construction_time():
    with pytest.raises(ValueError, match="max_block_frac"):
        BlockPolicy(table="hgnc.txt", max_block_frac=1.5)
    with pytest.raises(ValueError, match="table"):
        BlockPolicy(table="")
