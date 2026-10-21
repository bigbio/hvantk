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


def test_run_coloc_per_gene_window_is_contig_aware(tmp_path):
    """``prepare_coloc_data`` computed a ``contig`` field and then left it out of
    the columns it exports, so ``run_coloc_per_gene``'s per-gene +/-window filter
    compared ``pos`` only. A gene_id with rows on two chromosomes at the same
    numeric position (e.g. a pseudoautosomal gene) pulled the other chromosome's
    variant into the window.
    """
    import hail as hl

    from hvantk.algorithms.qtlcascade.coloc import run_coloc_per_gene

    def _row(contig, pos, p):
        return {
            "locus": hl.Locus(contig, pos, "GRCh38"),
            "alleles": ["A", "T"],
            "gene_id": "ENSG_PAR",
            "beta": 0.5,
            "se": 0.1,
            "p_value": p,
        }

    eqtl_rows = [
        _row("chr1", 1_000, 1e-10),  # lead eQTL variant (smallest p), chr1
        _row("chr1", 1_001, 0.5),  # second real variant, same contig as the lead
        _row("chr2", 1_000, 0.9),  # different chromosome, same numeric position
    ]
    pqtl_rows = [{k: v for k, v in r.items() if k != "p_value"} for r in eqtl_rows]

    eqtl_ht = hl.Table.parallelize(
        eqtl_rows,
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            gene_id=hl.tstr,
            beta=hl.tfloat64,
            se=hl.tfloat64,
            p_value=hl.tfloat64,
        ),
        key=["locus", "alleles", "gene_id"],
    )
    pqtl_ht = hl.Table.parallelize(
        pqtl_rows,
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            gene_id=hl.tstr,
            beta=hl.tfloat64,
            se=hl.tfloat64,
        ),
        key=["locus", "alleles", "gene_id"],
    )

    eqtl_path = str(tmp_path / "eqtl_allpairs.ht")
    pqtl_path = str(tmp_path / "pqtl_allpairs.ht")
    eqtl_ht.write(eqtl_path)
    pqtl_ht.write(pqtl_path)

    df = run_coloc_per_gene(
        eqtl_allpairs_ht_path=eqtl_path,
        pqtl_allpairs_ht_path=pqtl_path,
        cascade_genes=["ENSG_PAR"],
    )

    assert len(df) == 1
    # Only the two chr1 variants belong in ENSG_PAR's window; the chr2 variant
    # at the same numeric position must not be pulled in.
    assert int(df.iloc[0]["n_variants"]) == 2


def test_stage_coloc_writes_results_even_when_no_cascade_genes_qualify(tmp_path):
    """``_stage_coloc`` returned early, without writing ``coloc_results.tsv``,
    when no cascade gene had both an eQTL and a pQTL signal. An earlier run's
    file at the same per-tissue path was then silently left in place and read
    back as this run's result.
    """
    import hail as hl
    import pandas as pd

    from hvantk.algorithms.qtlcascade.pipeline import (
        CascadeConfig,
        CascadePipeline,
        CascadeResult,
    )

    ht = hl.Table.parallelize(
        [
            {"gene_id": "ENSG1", "cascade_class": "eqtl_only"},
            {"gene_id": "ENSG2", "cascade_class": "pqtl_only"},
        ],
        hl.tstruct(gene_id=hl.tstr, cascade_class=hl.tstr),
    )
    cascade_path = str(tmp_path / "cascade.ht")
    ht.write(cascade_path)

    config = CascadeConfig(
        eqtl_ht="unused_eqtl.ht",
        pqtl_ht="unused_pqtl.ht",
        eqtl_allpairs_ht="unused_eqtl_allpairs.ht",
        pqtl_allpairs_ht="unused_pqtl_allpairs.ht",
        output_dir=str(tmp_path / "out"),
        generate_plots=False,
        generate_report=False,
    )
    pipeline = CascadePipeline(config)

    # A previous run's coloc results, sitting at the exact path this run would
    # write to.
    stale_path = pipeline._per_tissue_dir / "Liver" / "coloc_results.tsv"
    stale_path.parent.mkdir(parents=True, exist_ok=True)
    stale_path.write_text(
        "gene_id\ttissue\tH0\tH1\tH2\tH3\tH4\tn_variants\n"
        "ENSG_OLD\tLiver\t0.0\t0.0\t0.0\t0.0\t1.0\t50\n"
    )

    result = CascadeResult(tissue="Liver", cascade_ht_path=cascade_path)
    result = pipeline._stage_coloc(result, tissue="Liver")

    written = pd.read_csv(stale_path, sep="\t")
    assert "ENSG_OLD" not in written["gene_id"].tolist()
    assert result.coloc_df is not None and result.coloc_df.empty
