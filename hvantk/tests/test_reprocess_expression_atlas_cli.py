"""Real end-to-end `hvantk reprocess expression-atlas:dataset` CliRunner tests.

Lives under hvantk/tests/, not hvantk/skills/expression_atlas/tests/, because it
imports hvantk.tools.plugins.reprocess_cli: skills/ must not import from tools/
(test_dependency_directions.py::test_skills_does_not_import_algorithms_or_tools
is a blanket, no-exceptions AST scan of every .py file under skills/, tests
included, and failed CI when this module first lived there).

Unlike hvantk/tests/test_reprocess_cli.py (mocked specs via an autouse,
registry-resetting fixture), this drives the REAL plugin registry and REAL
builder against the committed expression_atlas fixture -- only network is
avoided: --skip-download, since the fixture's files are already on disk, and a
stubbed drift probe.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import anndata as ad
import pandas as pd
from click.testing import CliRunner

from hvantk.core.plugin import loader as plugin_loader
from hvantk.tools.plugins.reprocess_cli import reprocess_cmd

FIXTURE_DIR = Path("hvantk/skills/expression_atlas/tests/testdata/raw/expression-atlas")
EXPRESSION = str(FIXTURE_DIR / "E-MTAB-6798-transcripts-tpms.tsv")
SDRF = str(FIXTURE_DIR / "E-MTAB-6798.condensed-sdrf.tsv")


# --- Regression: #383 review, "the documented reprocess command cannot build" ---------
# With no lifecycle.parse declared, `hvantk reprocess expression-atlas:dataset
# --skip-download` passed the bare --raw-dir STRING to build_expression_atlas, which
# indexed it like a dict (`parsed_input["expression_matrix"]`) and crashed with
# "string indices must be integers": the documented `--plugin-arg level=gene` route
# was unreachable through the real CLI.


def test_reprocess_cli_builds_from_the_raw_fixture_with_and_without_level_gene(
    tmp_path,
):
    """`hvantk reprocess expression-atlas:dataset --skip-download` must actually build."""
    # `--no-check-drift` only skips reprocess's post-build comparison;
    # run_builder_for_spec still calls spec.drift_probe() to stamp provenance, and
    # the real probe fetches the Expression Atlas experiment index. reprocess_cmd
    # fetches this same cached DatasetSpec, so stubbing it keeps every build below
    # offline.
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("expression-atlas:dataset")
    object.__setattr__(
        spec, "drift_probe", lambda: {"fingerprint": "sha256:test-expression-atlas"}
    )

    runner = CliRunner()

    out_transcript = tmp_path / "transcript.h5ad"
    result = runner.invoke(
        reprocess_cmd,
        [
            "expression-atlas:dataset",
            "--raw-dir",
            str(FIXTURE_DIR),
            "--output",
            str(out_transcript),
            "--skip-download",
            "--no-check-drift",
        ],
    )
    assert result.exit_code == 0, result.output
    built = ad.read_h5ad(out_transcript)
    assert built.var.index.name == "transcript_id"
    assert built.n_vars == 21, "21 transcripts in the fixture, default level=transcript"

    out_gene = tmp_path / "gene.h5ad"
    result = runner.invoke(
        reprocess_cmd,
        [
            "expression-atlas:dataset",
            "--raw-dir",
            str(FIXTURE_DIR),
            "--output",
            str(out_gene),
            "--skip-download",
            "--no-check-drift",
            "--plugin-arg",
            "level=gene",
        ],
    )
    assert result.exit_code == 0, result.output
    built_gene = ad.read_h5ad(out_gene)
    assert built_gene.var.index.name == "gene_id"
    assert built_gene.n_vars == 20, "21 transcripts collapse onto 20 distinct genes"

    # --- a per-run gene-level file is preferred over the transcript one -----------
    # A plain `--plugin-arg accession=...` download fetches every file in the
    # experiment's FTP directory (hvantk/skills/expression_atlas/cli.py, ftp.nlst()),
    # so both "*-tpms.tsv" files can be present. When both have the SDRF's runs as
    # columns, the gene-level one is used: it needs no aggregation.
    raw_dir = tmp_path / "both_files"
    raw_dir.mkdir()
    shutil.copy(EXPRESSION, raw_dir / "E-MTAB-6798-transcripts-tpms.tsv")
    shutil.copy(SDRF, raw_dir / "E-MTAB-6798.condensed-sdrf.tsv")
    # Fabricate a per-run gene-level file: one row per distinct gene (dropping the
    # transcript id column and the duplicate Cdc45 row).
    transcripts_df = pd.read_csv(EXPRESSION, sep="\t")
    gene_level_df = transcripts_df.drop(columns=["GeneID"]).drop_duplicates(
        subset=["Gene ID"], keep="first"
    )
    gene_level_df.to_csv(raw_dir / "E-MTAB-6798-tpms.tsv", sep="\t", index=False)

    out_preferred = tmp_path / "preferred.h5ad"
    result = runner.invoke(
        reprocess_cmd,
        [
            "expression-atlas:dataset",
            "--raw-dir",
            str(raw_dir),
            "--output",
            str(out_preferred),
            "--skip-download",
            "--no-check-drift",
        ],
    )
    assert result.exit_code == 0, result.output
    built_preferred = ad.read_h5ad(out_preferred)
    assert built_preferred.var.index.name == "gene_id"
    assert built_preferred.n_vars == 20, (
        "picked the gene-level file, not the transcript one"
    )

    # --- a gene-level file of assay-group summaries is passed over -----------------
    # E-MTAB-6798's real gene-level file has one column per assay group (g1 ... g94,
    # each cell a string of quartiles), none of them a run in the SDRF, so the build
    # must use the transcript-level file, whose columns are the runs.
    assay_group_dir = tmp_path / "assay_group_gene_file"
    assay_group_dir.mkdir()
    shutil.copy(EXPRESSION, assay_group_dir / "E-MTAB-6798-transcripts-tpms.tsv")
    shutil.copy(SDRF, assay_group_dir / "E-MTAB-6798.condensed-sdrf.tsv")
    (assay_group_dir / "E-MTAB-6798-tpms.tsv").write_text(
        "GeneID\tGene Name\tg1\tg2\n"
        "ENSMUSG00000000001\tGnai3\t268,272,276,314,351\t12,15,16,19,24\n"
        "ENSMUSG00000000003\tPbsn\t0,0,0,0,0\t0,0,0,0,1\n"
    )
    out_assay_group = tmp_path / "assay_group.h5ad"
    result = runner.invoke(
        reprocess_cmd,
        [
            "expression-atlas:dataset",
            "--raw-dir",
            str(assay_group_dir),
            "--output",
            str(out_assay_group),
            "--skip-download",
            "--no-check-drift",
        ],
    )
    assert result.exit_code == 0, result.output
    built_assay_group = ad.read_h5ad(out_assay_group)
    assert built_assay_group.var.index.name == "transcript_id"
    assert built_assay_group.n_vars == 21, (
        "picked the transcript file, not the assay-group gene file"
    )
    assert built_assay_group.n_obs == 4

    # --- mismatched accession must raise, not silently combine (PR #455 review) -----
    # SDRF and TPM candidates used to be counted independently ("exactly one of
    # each"), so a --raw-dir holding one accession's SDRF and a DIFFERENT accession's
    # TPM file (a reused directory, two downloads merged by hand, ...) silently
    # passed and would have combined mismatched expression + metadata.
    mismatched_dir = tmp_path / "mismatched_accession"
    mismatched_dir.mkdir()
    shutil.copy(SDRF, mismatched_dir / "E-MTAB-6798.condensed-sdrf.tsv")
    shutil.copy(
        EXPRESSION, mismatched_dir / "E-MTAB-9999-transcripts-tpms.tsv"
    )  # a different accession
    result = runner.invoke(
        reprocess_cmd,
        [
            "expression-atlas:dataset",
            "--raw-dir",
            str(mismatched_dir),
            "--output",
            str(tmp_path / "mismatched.h5ad"),
            "--skip-download",
            "--no-check-drift",
        ],
    )
    assert result.exit_code != 0
    assert isinstance(result.exception, FileNotFoundError)
    assert "E-MTAB-6798" in str(result.exception)
    assert "same accession" in str(result.exception)

    # --- several gene-level candidates must raise, not fall back to transcript ------
    # A second gene-level-shaped file for the same accession makes the gene-level
    # pick ambiguous; falling back to the lone transcript-level file in that case
    # would silently prefer it over an ambiguous gene-level set.
    ambiguous_dir = tmp_path / "ambiguous_gene_level"
    ambiguous_dir.mkdir()
    shutil.copy(SDRF, ambiguous_dir / "E-MTAB-6798.condensed-sdrf.tsv")
    shutil.copy(EXPRESSION, ambiguous_dir / "E-MTAB-6798-transcripts-tpms.tsv")
    gene_level_df.to_csv(ambiguous_dir / "E-MTAB-6798-tpms.tsv", sep="\t", index=False)
    gene_level_df.to_csv(
        ambiguous_dir / "E-MTAB-6798-old-tpms.tsv", sep="\t", index=False
    )
    result = runner.invoke(
        reprocess_cmd,
        [
            "expression-atlas:dataset",
            "--raw-dir",
            str(ambiguous_dir),
            "--output",
            str(tmp_path / "ambiguous.h5ad"),
            "--skip-download",
            "--no-check-drift",
        ],
    )
    assert result.exit_code != 0
    assert isinstance(result.exception, FileNotFoundError)
    assert "could not find exactly one gene-level" in str(result.exception)
