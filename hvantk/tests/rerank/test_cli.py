# hvantk/tests/rerank/test_cli.py
import logging

import numpy as np, pandas as pd, yaml
from click.testing import CliRunner
from hvantk.tools.rerank import rerank_cmd

_RERANK_CLI_LOGGER = "hvantk.tools.rerank.rerank_cli"


def _toy_fixtures(tmp_path, n=150, extra_cohort_cols=None):
    """The synthetic gene/label/feature fixtures shared by the CLI tests below.

    ``extra_cohort_cols`` is a dict of {column_name: array-like} merged into the
    cohort table alongside the prior column, so callers can add case/control
    architecture columns without duplicating the rest of the setup.
    """
    rng = np.random.default_rng(0)
    genes = [f"g{i}" for i in range(n)]
    y = (rng.random(n) < 0.3).astype(int)
    cohort_cols = {"gene": genes, "p": rng.random(n)}
    if extra_cohort_cols:
        cohort_cols.update(extra_cohort_cols)
    pd.DataFrame(cohort_cols).to_csv(tmp_path / "cohort.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": genes, "x": y + rng.normal(0, 0.5, n)}).to_parquet(
        tmp_path / "feat.parquet"
    )
    (tmp_path / "labels.txt").write_text(
        "\n".join(g for g, yy in zip(genes, y) if yy == 1)
    )
    return genes, y


def _write_config(tmp_path, cohort_spec):
    (tmp_path / "cohort.yaml").write_text(yaml.safe_dump(cohort_spec))
    cfg = {
        "name": "toy",
        "cohort": str(tmp_path / "cohort.yaml"),
        "features": [{"name": "constraint", "path": str(tmp_path / "feat.parquet")}],
        "labels": {"path": str(tmp_path / "labels.txt")},
        "min_label_coverage": 0.0,
    }
    (tmp_path / "c.yaml").write_text(yaml.safe_dump(cfg))
    return tmp_path / "c.yaml"


def test_rerank_cli_end_to_end(tmp_path, caplog):
    n = 150
    genes, y = _toy_fixtures(tmp_path, n=n)
    cohort_spec = {
        "name": "toy",
        "key": "symbol",
        "key_column": "gene",
        "table": str(tmp_path / "cohort.tsv"),
        "prior": {"column": "p", "direction": "lower_is_better"},
    }
    config_path = _write_config(tmp_path, cohort_spec)
    out = tmp_path / "out.tsv"
    with caplog.at_level(logging.WARNING, logger=_RERANK_CLI_LOGGER):
        r = CliRunner().invoke(rerank_cmd, ["-c", str(config_path), "-o", str(out)])
    assert r.exit_code == 0, r.output
    t = pd.read_csv(out, sep="\t")
    assert len(t) == n and {"gene", "score", "verdict"} <= set(t.columns)
    # No cohort_axes entry named "architecture" is declared at all -- the near-miss
    # warning must stay silent here, otherwise it becomes noise on every plain,
    # prior-only cohort (the common case).
    assert [r for r in caplog.records if r.name == _RERANK_CLI_LOGGER] == []


def test_rerank_cli_warns_on_a_near_miss_architecture_axis(tmp_path, caplog):
    # A plausible rename typo: the axis is named "architecture" and the column really
    # exists in the table (so manifest validation passes), but it's spelled
    # "driver_afx" rather than "driver_af" -- CaseControlArchitectureAudit never gets
    # wired up, and before this fix nothing said why.
    n = 150
    genes, y = _toy_fixtures(
        tmp_path,
        n=n,
        extra_cohort_cols={
            "n_case_var": np.full(n, 10),
            "conc": np.full(n, 0.1),
            "driver_afx": np.zeros(n),
        },
    )
    cohort_spec = {
        "name": "toy",
        "key": "symbol",
        "key_column": "gene",
        "table": str(tmp_path / "cohort.tsv"),
        "prior": {"column": "p", "direction": "lower_is_better"},
        "cohort_axes": [
            {"axis": "architecture", "columns": ["n_case_var", "conc", "driver_afx"]}
        ],
    }
    config_path = _write_config(tmp_path, cohort_spec)
    out = tmp_path / "out.tsv"
    with caplog.at_level(logging.WARNING, logger=_RERANK_CLI_LOGGER):
        r = CliRunner().invoke(rerank_cmd, ["-c", str(config_path), "-o", str(out)])
    assert r.exit_code == 0, r.output

    warnings = [rec for rec in caplog.records if rec.name == _RERANK_CLI_LOGGER]
    assert len(warnings) == 1
    message = warnings[0].getMessage()
    assert "toy" in message  # names the cohort
    assert "architecture" in message  # names the axis
    assert "driver_af" in message  # names the missing required column

    # The near-miss falls back to NoAudit, exactly like a cohort with no architecture
    # axis at all -- nothing is flagged.
    t = pd.read_csv(out, sep="\t")
    assert not t["flag"].any()


def test_rerank_cli_wires_up_the_architecture_audit_when_declared(tmp_path):
    # A cohort_axes entry that genuinely covers all three columns
    # CaseControlArchitectureAudit needs must make the CLI wire up the real audit, not
    # NoAudit -- crafted so gene "g0" trips the insufficient_data rule (n_case_var <=
    # 2, per hvantk/algorithms/rerank/audit.py) and nothing else does.
    n = 150
    n_case_var = np.full(n, 10)
    n_case_var[0] = 2
    genes, y = _toy_fixtures(
        tmp_path,
        n=n,
        extra_cohort_cols={
            "n_case_var": n_case_var,
            "conc": np.full(n, 0.1),
            "driver_af": np.zeros(n),
        },
    )
    cohort_spec = {
        "name": "toy",
        "key": "symbol",
        "key_column": "gene",
        "table": str(tmp_path / "cohort.tsv"),
        "prior": {"column": "p", "direction": "lower_is_better"},
        "cohort_axes": [
            {"axis": "architecture", "columns": ["n_case_var", "conc", "driver_af"]}
        ],
    }
    config_path = _write_config(tmp_path, cohort_spec)
    out = tmp_path / "out.tsv"
    r = CliRunner().invoke(rerank_cmd, ["-c", str(config_path), "-o", str(out)])
    assert r.exit_code == 0, r.output

    t = pd.read_csv(out, sep="\t")
    # NoAudit would flag nothing, ever -- a flagged gene here is only possible if
    # CaseControlArchitectureAudit genuinely ran.
    assert t["flag"].any()
    flagged = t[t["flag"]]
    assert "g0" in set(flagged["gene"])
    assert set(flagged["flag_reason"]) <= {
        "insufficient_data",
        "recurrent_variant",
        "common_driver",
    }


def test_rerank_cli_rejects_a_stale_prior_block(tmp_path):
    # Finding 4 (whole-branch review): a config migrated to declare 'cohort:' that
    # still carries a leftover pre-migration 'prior:' block must fail loud, not be
    # silently ignored -- the prior it names is never read, and prior_stat would come
    # from the cohort manifest's own table instead, with no warning at all.
    n = 150
    genes, y = _toy_fixtures(tmp_path, n=n)
    cohort_spec = {
        "name": "toy",
        "key": "symbol",
        "key_column": "gene",
        "table": str(tmp_path / "cohort.tsv"),
        "prior": {"column": "p", "direction": "lower_is_better"},
    }
    (tmp_path / "cohort.yaml").write_text(yaml.safe_dump(cohort_spec))
    cfg = {
        "name": "toy",
        "cohort": str(tmp_path / "cohort.yaml"),
        "features": [{"name": "constraint", "path": str(tmp_path / "feat.parquet")}],
        "labels": {"path": str(tmp_path / "labels.txt")},
        "min_label_coverage": 0.0,
        # Stale: pre-migration key, points at a file this run never touches.
        "prior": {
            "path": str(tmp_path / "my_curated_burden.tsv"),
            "unit_col": "gene",
            "stat_col": "minp",
        },
    }
    config_path = tmp_path / "c.yaml"
    config_path.write_text(yaml.safe_dump(cfg))
    out = tmp_path / "out.tsv"
    r = CliRunner().invoke(rerank_cmd, ["-c", str(config_path), "-o", str(out)])
    assert r.exit_code != 0
    assert "prior" in r.output
    assert not out.exists()


def test_rerank_cli_rejects_any_unknown_top_level_key(tmp_path):
    n = 150
    genes, y = _toy_fixtures(tmp_path, n=n)
    cohort_spec = {
        "name": "toy",
        "key": "symbol",
        "key_column": "gene",
        "table": str(tmp_path / "cohort.tsv"),
        "prior": {"column": "p", "direction": "lower_is_better"},
    }
    config_path = _write_config(tmp_path, cohort_spec)
    cfg = yaml.safe_load(config_path.read_text())
    cfg["not_a_real_key"] = True
    config_path.write_text(yaml.safe_dump(cfg))
    out = tmp_path / "out.tsv"
    r = CliRunner().invoke(rerank_cmd, ["-c", str(config_path), "-o", str(out)])
    assert r.exit_code != 0
    assert "not_a_real_key" in r.output


def test_rerank_cli_no_near_miss_warning_when_architecture_columns_are_split(
    tmp_path, caplog
):
    # Finding 5 (whole-branch review): a manifest splitting the three architecture
    # columns across an "architecture" axis (n_case_var, conc) and a differently-named
    # axis (driver_af here, under "qc") is a WORKING configuration --
    # CaseControlArchitectureAudit genuinely wires up (has_architecture_columns()
    # looks at declared_columns() as a whole, not any one axis). The near-miss warning
    # must stay silent here; it used to fire unconditionally before has_architecture
    # was computed.
    n = 150
    n_case_var = np.full(n, 10)
    n_case_var[0] = 2
    genes, y = _toy_fixtures(
        tmp_path,
        n=n,
        extra_cohort_cols={
            "n_case_var": n_case_var,
            "conc": np.full(n, 0.1),
            "driver_af": np.zeros(n),
        },
    )
    cohort_spec = {
        "name": "toy",
        "key": "symbol",
        "key_column": "gene",
        "table": str(tmp_path / "cohort.tsv"),
        "prior": {"column": "p", "direction": "lower_is_better"},
        "cohort_axes": [
            {"axis": "architecture", "columns": ["n_case_var", "conc"]},
            {"axis": "qc", "columns": ["driver_af"]},
        ],
    }
    config_path = _write_config(tmp_path, cohort_spec)
    out = tmp_path / "out.tsv"
    with caplog.at_level(logging.WARNING, logger=_RERANK_CLI_LOGGER):
        r = CliRunner().invoke(rerank_cmd, ["-c", str(config_path), "-o", str(out)])
    assert r.exit_code == 0, r.output

    assert [rec for rec in caplog.records if rec.name == _RERANK_CLI_LOGGER] == []

    # The audit genuinely ran (not NoAudit): gene "g0" trips insufficient_data.
    t = pd.read_csv(out, sep="\t")
    assert t["flag"].any()
    assert "g0" in set(t.loc[t["flag"], "gene"])


def test_rerank_cli_falls_back_cleanly_when_the_prior_column_is_the_missing_third(
    tmp_path, caplog
):
    # Finding 3 (re-review): has_architecture_columns() used to be fed
    # declared_columns() (prior column + axis columns), but engine.rerank() merges
    # the cohort frame with include_prior=False -- the two disagree by exactly the
    # prior column. Reproduction: prior.column='n_case_var', and a cohort_axes
    # "architecture" entry declaring only the other two required columns ('conc',
    # 'driver_af'). declared_columns() looks like a superset of the three required
    # columns (it counts the prior column too), so the audit used to get wired up --
    # then engine.rerank() merged only the axis columns, CaseControlArchitectureAudit
    # never saw 'n_case_var', and it raised "requires case/control architecture
    # columns ['n_case_var']" for a manifest the user believed declared all three.
    # Post-fix, eligibility is judged on axis_columns() alone, which agrees with what
    # the engine actually merges: this manifest is correctly NOT eligible, falls back
    # to NoAudit, and the near-miss warning explains why -- no raise-after-fitting.
    n = 150
    genes, y = _toy_fixtures(
        tmp_path,
        n=n,
        extra_cohort_cols={
            "n_case_var": np.full(n, 10),
            "conc": np.full(n, 0.1),
            "driver_af": np.zeros(n),
        },
    )
    cohort_spec = {
        "name": "toy",
        "key": "symbol",
        "key_column": "gene",
        "table": str(tmp_path / "cohort.tsv"),
        "prior": {"column": "n_case_var", "direction": "lower_is_better"},
        "cohort_axes": [{"axis": "architecture", "columns": ["conc", "driver_af"]}],
    }
    config_path = _write_config(tmp_path, cohort_spec)
    out = tmp_path / "out.tsv"
    with caplog.at_level(logging.WARNING, logger=_RERANK_CLI_LOGGER):
        r = CliRunner().invoke(rerank_cmd, ["-c", str(config_path), "-o", str(out)])
    assert r.exit_code == 0, r.output
    assert "requires case/control architecture columns" not in r.output

    warnings = [rec for rec in caplog.records if rec.name == _RERANK_CLI_LOGGER]
    assert len(warnings) == 1
    message = warnings[0].getMessage()
    assert "architecture" in message
    assert "n_case_var" in message

    # NoAudit ran, not CaseControlArchitectureAudit -- nothing is flagged.
    t = pd.read_csv(out, sep="\t")
    assert not t["flag"].any()
    assert len(t) == n


def test_rerank_cli_reports_a_missing_cohort_manifest_clearly(tmp_path):
    cfg = {
        "name": "toy",
        "cohort": str(tmp_path / "does-not-exist.yaml"),
        "features": [],
        "labels": {"path": str(tmp_path / "labels.txt")},
    }
    (tmp_path / "c.yaml").write_text(yaml.safe_dump(cfg))
    r = CliRunner().invoke(
        rerank_cmd, ["-c", str(tmp_path / "c.yaml"), "-o", str(tmp_path / "out.tsv")]
    )
    assert r.exit_code != 0
    assert "does-not-exist.yaml" in r.output
    assert "Traceback" not in r.output
