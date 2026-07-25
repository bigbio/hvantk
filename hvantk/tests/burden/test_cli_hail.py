import pytest
import pandas as pd
import hail as hl
from click.testing import CliRunner

from hvantk.tools.cohort.cohort_cli import cohort_group
from hvantk.tests.burden.test_aggregate_hail import _toy_mt

pytestmark = pytest.mark.hail


def test_cohort_burden_cli_writes_gene_table(tmp_path):
    mt_path = str(tmp_path / "cohort.mt")
    _toy_mt().write(mt_path, overwrite=True)
    out_path = str(tmp_path / "genes.tsv")
    res = CliRunner().invoke(
        cohort_group,
        [
            "burden",
            "--mt",
            mt_path,
            "--gene-col",
            "SYMBOL",
            "--route-col",
            "csq_group",
            "--arm-col",
            "is_case",
            "--key",
            "symbol",
            "--mtc",
            "bh",
            "--output",
            out_path,
        ],
    )
    assert res.exit_code == 0, res.output
    df = pd.read_csv(out_path, sep="\t")
    assert {"gene", "route", "minp", "n_case_var", "conc", "driver_af", "p_adj"} <= set(
        df.columns
    )
    assert "GENEA" in set(df["gene"])


def test_cohort_burden_cli_emit_manifest(tmp_path):
    """--emit-manifest writes a stub whose declared columns are present in the output."""
    from hvantk.algorithms.cohort.spec import load_cohort

    mt_path = str(tmp_path / "cohort.mt")
    _toy_mt().write(mt_path, overwrite=True)
    out_path = str(tmp_path / "burden.tsv")
    man_path = str(tmp_path / "cohort.yaml")
    res = CliRunner().invoke(
        cohort_group,
        [
            "burden",
            "--mt",
            mt_path,
            "--gene-col",
            "SYMBOL",
            "--route-col",
            "csq_group",
            "--arm-col",
            "is_case",
            "--key",
            "symbol",
            "--output",
            out_path,
            "--emit-manifest",
            man_path,
        ],
    )
    assert res.exit_code == 0, res.output
    manifest = load_cohort(man_path)
    assert manifest.prior.column == "minp"
    assert manifest.key == "symbol" and manifest.key_column == "gene"
    assert manifest.table == out_path
    # the emitted manifest's declared columns really are in the burden output (the
    # hand-off guarantee: `cohort burden --emit-manifest` -> a valid, consumable cohort)
    cols = set(pd.read_csv(out_path, sep="\t").columns)
    assert set(manifest.declared_columns()) <= cols
