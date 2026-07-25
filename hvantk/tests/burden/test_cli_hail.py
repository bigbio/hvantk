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
