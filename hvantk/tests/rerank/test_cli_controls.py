"""#247 reaching the command line.

A control that exists only as a Config field is a control nobody running `hvantk rerank`
can use, which is the state leakage.py (#244) was left in.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import yaml
from click.testing import CliRunner

from hvantk.tools.rerank import rerank_cmd
from hvantk.tests.rerank.test_cli import _toy_fixtures, _write_config


def _cohort(tmp_path):
    return {
        "name": "toy",
        "key": "symbol",
        "key_column": "gene",
        "table": str(tmp_path / "cohort.tsv"),
        "prior": {"column": "p", "direction": "lower_is_better"},
    }


def _run(tmp_path, *extra):
    # `tmp_path / "a"`, `tmp_path / "b"` etc. (used by tests that need two isolated
    # invocations) are not created by the `tmp_path` fixture itself.
    tmp_path.mkdir(parents=True, exist_ok=True)
    _toy_fixtures(tmp_path, n=150)
    cfg = _write_config(tmp_path, _cohort(tmp_path))
    out = tmp_path / "out.tsv"
    return CliRunner().invoke(rerank_cmd, ["-c", str(cfg), "-o", str(out), *extra]), out


def _write_two_axis_config(tmp_path, genes, y):
    """A config with a SECOND feature axis, for the null-specific tests below.

    `_toy_fixtures`/`_write_config` (shared with test_cli.py, and reused by every other
    test in this file) declare exactly one feature axis. That axis is then also the
    baseline, so there is no candidate axis left for a permutation null to search over
    and `--n-perm` is refused (see the single-axis test below). Written locally, rather
    than by editing the shared fixture, which every other rerank CLI test relies on
    staying single-axis.
    """
    rng = np.random.default_rng(1)
    pd.DataFrame({"gene": genes, "z": y + rng.normal(0, 0.5, len(genes))}).to_parquet(
        tmp_path / "feat2.parquet"
    )
    (tmp_path / "cohort.yaml").write_text(yaml.safe_dump(_cohort(tmp_path)))
    spec = {
        "name": "toy",
        "cohort": str(tmp_path / "cohort.yaml"),
        "features": [
            {"name": "constraint", "path": str(tmp_path / "feat.parquet")},
            {"name": "expression", "path": str(tmp_path / "feat2.parquet")},
        ],
        "labels": {"path": str(tmp_path / "labels.txt")},
        "min_label_coverage": 0.0,
    }
    (tmp_path / "c.yaml").write_text(yaml.safe_dump(spec))
    return tmp_path / "c.yaml"


def test_the_defaults_are_unchanged(tmp_path):
    r, out = _run(tmp_path)
    assert r.exit_code == 0, r.output
    assert len(pd.read_csv(out, sep="\t")) == 150


def test_seed_is_a_flag(tmp_path):
    a, out_a = _run(tmp_path / "a", "--seed", "1")
    b, out_b = _run(tmp_path / "b", "--seed", "2")
    assert a.exit_code == 0 and b.exit_code == 0, (a.output, b.output)
    sa = pd.read_csv(out_a, sep="\t").score.round(9).tolist()
    sb = pd.read_csv(out_b, sep="\t").score.round(9).tolist()
    assert sa != sb, "--seed did not reach the partition"


def test_seed_sweep_adds_the_envelope_columns_only_when_sweeping(tmp_path):
    plain, _ = _run(tmp_path / "p")
    swept, _ = _run(tmp_path / "s", "--seed-sweep", "3")
    assert swept.exit_code == 0, swept.output
    assert "d_lo_env" in swept.output
    assert "d_lo_env" not in plain.output


def test_n_perm_prints_a_selection_corrected_p_and_null_out_writes_the_summary(
    tmp_path,
):
    """One CLI run covers both --n-perm's console output and --null-out's TSV: the two were
    previously separate tests that each paid for their own permutation run over the same
    two-axis config."""
    genes, y = _toy_fixtures(tmp_path, n=150)
    cfg = _write_two_axis_config(tmp_path, genes, y)
    nul = tmp_path / "null.tsv"
    r = CliRunner().invoke(
        rerank_cmd,
        [
            "-c",
            str(cfg),
            "-o",
            str(tmp_path / "out.tsv"),
            "--n-perm",
            "5",
            "--null-out",
            str(nul),
        ],
    )
    assert r.exit_code == 0, r.output
    assert "p_selected_max" in r.output
    summary = pd.read_csv(nul, sep="\t")
    assert {"axis", "n_perm", "n_candidates", "selmax_median"} <= set(summary.columns)
    # No --blocks: the unblocked null is anti-conservative wherever labels cluster by
    # family, and the console must say so ONCE, on stderr (never inside a stdout table).
    assert r.stderr.count("anti-conservative") == 1, r.stderr
    assert "--blocks" in r.stderr


def test_n_perm_with_blocks_runs_a_blocked_null_and_does_not_warn(tmp_path):
    """The unblocked-null warning is about the MISSING blocking; with --blocks the null
    permutes by block and the warning must stay silent -- otherwise it becomes noise on
    exactly the run that did the right thing."""
    genes, y = _toy_fixtures(tmp_path, n=150)
    cfg = _write_two_axis_config(tmp_path, genes, y)
    hgnc = tmp_path / "hgnc.txt"
    rows = ["symbol\tgene_group\tstatus"]
    for i, g in enumerate(genes):
        rows.append(f"{g}\tFamily {i // 5}\tApproved")
    hgnc.write_text("\n".join(rows) + "\n")
    r = CliRunner().invoke(
        rerank_cmd,
        [
            "-c",
            str(cfg),
            "-o",
            str(tmp_path / "out.tsv"),
            "--n-perm",
            "2",
            "--blocks",
            str(hgnc),
        ],
    )
    assert r.exit_code == 0, r.output
    assert "p_selected_max" in r.output
    assert "block_digest='" in r.output  # the Setting line records the blocking
    assert "anti-conservative" not in r.stderr, r.stderr


def test_n_perm_on_a_single_axis_config_fails_cleanly(tmp_path):
    """The Minor 6 guard: `_write_config`/`_toy_fixtures` (shared with `_run` above)
    declare exactly one feature axis, which is then also the baseline, so there is no
    candidate axis left for `--n-perm` to search over. This must be a clean CLI failure
    before any output is written, not a silent no-op (the engine used to return
    `res.nulls is None` and the CLI just skipped the Multiplicity block). It is knowable
    from the config alone, so the CLI rejects it before any work starts."""
    r, out = _run(tmp_path, "--n-perm", "3")
    assert r.exit_code != 0
    assert "Traceback" not in r.output
    assert "--n-perm" in r.output
    assert not out.exists(), "a failed run must not leave a scored table behind"


def test_n_perm_with_an_axis_that_contributes_no_column_fails_cleanly(tmp_path):
    """The engine's own refusal must reach the user as a clean CLI error, not a traceback
    and not a scored table with the correction silently missing. This config passes the
    one-axis preflight (it names two axes), but the second table carries only a 'gene'
    column, so after scoring the engine is left with the baseline alone -- the case the
    CLI cannot see from the config and `engine._run_nulls` raises on."""
    genes, y = _toy_fixtures(tmp_path, n=150)
    pd.DataFrame({"gene": genes}).to_parquet(tmp_path / "empty.parquet")
    (tmp_path / "cohort.yaml").write_text(yaml.safe_dump(_cohort(tmp_path)))
    spec = {
        "name": "toy",
        "cohort": str(tmp_path / "cohort.yaml"),
        "features": [
            {"name": "constraint", "path": str(tmp_path / "feat.parquet")},
            {"name": "expression", "path": str(tmp_path / "empty.parquet")},
        ],
        "labels": {"path": str(tmp_path / "labels.txt")},
        "min_label_coverage": 0.0,
    }
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump(spec))
    out = tmp_path / "out.tsv"
    r = CliRunner().invoke(
        rerank_cmd, ["-c", str(cfg), "-o", str(out), "--n-perm", "3"]
    )
    assert r.exit_code != 0
    assert "Traceback" not in r.output
    assert "'expression'" in r.output and "no feature column" in r.output, r.output
    assert not out.exists(), "a failed run must not leave a scored table behind"


def test_blocks_flag_wires_the_block_builder(tmp_path):
    genes, _ = _toy_fixtures(tmp_path, n=150)
    hgnc = tmp_path / "hgnc.txt"
    rows = ["symbol\tgene_group\tstatus"]
    for i, g in enumerate(genes):
        rows.append(f"{g}\tFamily {i // 5}\tApproved")
    hgnc.write_text("\n".join(rows) + "\n")
    cfg = _write_config(tmp_path, _cohort(tmp_path))
    r = CliRunner().invoke(
        rerank_cmd,
        ["-c", str(cfg), "-o", str(tmp_path / "out.tsv"), "--blocks", str(hgnc)],
    )
    assert r.exit_code == 0, r.output
    assert "blocks" in r.output.lower()


def test_a_dominant_block_fails_the_command_rather_than_scoring(tmp_path):
    """The abort must reach the exit code. A run that continued past it would produce a
    plausible table whose pooled OOF AUC estimates something else."""
    genes, _ = _toy_fixtures(tmp_path, n=150)
    hgnc = tmp_path / "hgnc.txt"
    rows = ["symbol\tgene_group\tstatus"]
    for i, g in enumerate(genes):
        rows.append(f"{g}\t{'Huge' if i < 100 else f'Family {i}'}\tApproved")
    hgnc.write_text("\n".join(rows) + "\n")
    cfg = _write_config(tmp_path, _cohort(tmp_path))
    out = tmp_path / "out.tsv"
    r = CliRunner().invoke(
        rerank_cmd, ["-c", str(cfg), "-o", str(out), "--blocks", str(hgnc)]
    )
    assert r.exit_code != 0
    assert "ceiling" in r.output or "max_block_frac" in r.output
    assert "Traceback" not in r.output
    assert not out.exists(), "a failed run must not leave a scored table behind"


def test_max_block_frac_can_accept_the_dominant_block_deliberately(tmp_path):
    genes, _ = _toy_fixtures(tmp_path, n=150)
    hgnc = tmp_path / "hgnc.txt"
    rows = ["symbol\tgene_group\tstatus"]
    for i, g in enumerate(genes):
        rows.append(f"{g}\t{'Huge' if i < 20 else f'Family {i // 5}'}\tApproved")
    hgnc.write_text("\n".join(rows) + "\n")
    cfg = _write_config(tmp_path, _cohort(tmp_path))
    r = CliRunner().invoke(
        rerank_cmd,
        [
            "-c",
            str(cfg),
            "-o",
            str(tmp_path / "out.tsv"),
            "--blocks",
            str(hgnc),
            "--max-block-frac",
            "0.5",
        ],
    )
    assert r.exit_code == 0, r.output


@pytest.mark.parametrize(
    "flag",
    [
        "--seed",
        "--seed-sweep",
        "--blocks",
        "--max-block-frac",
        "--n-perm",
        "--null-out",
    ],
)
def test_every_new_flag_is_documented_in_help(flag):
    r = CliRunner().invoke(rerank_cmd, ["--help"])
    assert flag in r.output
    if flag == "--max-block-frac":
        # The help text hardcodes the default as a literal (avoids importing the rerank
        # package, and scikit-learn with it, just to build a help string); this assertion
        # is what stops that literal drifting from the real runtime default.
        from hvantk.algorithms.rerank.blocks import DEFAULT_MAX_BLOCK_FRAC

        assert f"{DEFAULT_MAX_BLOCK_FRAC:g}" in r.output


def test_a_gene_group_table_matching_no_gene_fails_cleanly(tmp_path):
    """A gene-group table with no symbols in the universe must fail cleanly with a
    user-facing error, not a traceback. The run must not write a scored table."""
    genes, _ = _toy_fixtures(tmp_path, n=150)
    hgnc = tmp_path / "hgnc.txt"
    # Write a table with symbols that don't exist in the universe
    rows = ["symbol\tgene_group\tstatus"]
    rows.append("NOTAGENE1\tFamily1\tApproved")
    rows.append("NOTAGENE2\tFamily2\tApproved")
    rows.append("NOTAGENE3\tFamily3\tApproved")
    hgnc.write_text("\n".join(rows) + "\n")
    cfg = _write_config(tmp_path, _cohort(tmp_path))
    out = tmp_path / "out.tsv"
    r = CliRunner().invoke(
        rerank_cmd, ["-c", str(cfg), "-o", str(out), "--blocks", str(hgnc)]
    )
    assert r.exit_code != 0
    assert "Traceback" not in r.output
    assert "none of" in r.output
    assert not out.exists(), "a failed run must not leave a scored table behind"


def test_the_docs_page_documents_every_new_flag():
    """docs_site is where a user looks first, and a flag that only exists in --help is a
    flag only someone who already knew about it will find."""
    from pathlib import Path

    import hvantk

    doc = (
        Path(hvantk.__file__).resolve().parents[1] / "docs_site" / "tools" / "rerank.md"
    )
    text = doc.read_text()
    for flag in (
        "--seed",
        "--seed-sweep",
        "--blocks",
        "--max-block-frac",
        "--n-perm",
        "--null-out",
    ):
        assert flag in text, flag
    assert "selected-maximum" in text
