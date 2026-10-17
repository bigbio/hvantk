"""Regression test for the coloc stage's gene selection in ``qtlcascade run``.

``CascadePipeline._stage_coloc`` aggregated ``ht.gene_id`` over ``ht.filter(...)``,
an expression from a different Hail source, so every ``qtlcascade run`` given
all-pairs tables failed with ``'Table.aggregate': source mismatch``. The CLI tests
mock the pipeline, so the defect only showed on real tables.
"""

import pytest

pytestmark = pytest.mark.hail


def test_coloc_gene_ids_reads_the_filtered_table():
    import hail as hl

    from hvantk.algorithms.qtlcascade.pipeline import coloc_gene_ids

    ht = hl.Table.parallelize(
        [
            {"gene_id": "ENSG3", "cascade_class": "eqtl_mediated"},
            {"gene_id": "ENSG1", "cascade_class": "discordant"},
            {"gene_id": "ENSG3", "cascade_class": "discordant"},
            {"gene_id": "ENSG2", "cascade_class": "eqtl_only"},
            {"gene_id": "ENSG4", "cascade_class": "pqtl_only"},
        ],
        hl.tstruct(gene_id=hl.tstr, cascade_class=hl.tstr),
    )
    assert coloc_gene_ids(ht) == ["ENSG1", "ENSG3"]
