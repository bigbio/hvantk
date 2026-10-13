"""Real end-to-end `hvantk reprocess expression-atlas:dataset` CliRunner tests.

Lives under hvantk/tests/, not hvantk/skills/expression_atlas/tests/, because it
imports hvantk.tools.plugins.reprocess_cli: skills/ must not import from tools/
(test_dependency_directions.py::test_skills_does_not_import_algorithms_or_tools
is a blanket, no-exceptions AST scan of every .py file under skills/, tests
included, and failed CI when this module first lived there).

Unlike hvantk/tests/test_reprocess_cli.py (mocked specs via an autouse,
registry-resetting fixture), this drives the REAL plugin registry and REAL
builder against the committed expression_atlas fixture -- only network is
avoided, via --skip-download, since the fixture's files are already on disk.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import anndata as ad
import pandas as pd
from click.testing import CliRunner

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

    # --- gene-level file preferred when both are present (#383 review follow-up) ---
    # A plain `--plugin-arg accession=...` download fetches every file in the
    # experiment's FTP directory (hvantk/skills/expression_atlas/cli.py, ftp.nlst()).
    # Confirmed live against E-MTAB-6798 itself: its directory ships BOTH
    # E-MTAB-6798-tpms.tsv (gene-level) and E-MTAB-6798-transcripts-tpms.tsv
    # (transcript-level), so the documented, unfiltered `reprocess ... --plugin-arg
    # accession=E-MTAB-6798` must build from the gene-level one without needing
    # level=gene at all -- not raise on finding two "*-tpms.tsv" candidates.
    raw_dir = tmp_path / "both_files"
    raw_dir.mkdir()
    shutil.copy(EXPRESSION, raw_dir / "E-MTAB-6798-transcripts-tpms.tsv")
    shutil.copy(SDRF, raw_dir / "E-MTAB-6798.condensed-sdrf.tsv")
    # Fabricate the gene-level file the real accession also ships: one row per
    # distinct gene (dropping the transcript id column and the duplicate Cdc45 row).
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
