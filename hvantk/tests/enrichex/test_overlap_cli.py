"""
Tests for the EnrichEx `overlap` CLI command.

Overlap enrichment uses scipy, not Hail, so these exercise the Click
command end-to-end (CliRunner) with no Hail session needed.
"""

from pathlib import Path

from click.testing import CliRunner

from hvantk.core.utils.gene_sets import GeneSet, GeneSetCollection
from hvantk.tools.enrichex.overlap_cli import overlap_test


def _write_query_genes(path: Path) -> Path:
    path.write_text("A\nB\nC\n")
    return path


def _write_gene_sets(path: Path) -> Path:
    gs = GeneSet(name="set1", genes={"A", "B", "C"})
    collection = GeneSetCollection(
        gene_sets={gs.name: gs},
        background_genes=set("ABCDEFGHIJ"),
        source_description="Test gene sets",
    )
    collection.save(path)
    return path


def test_overlap_cli_creates_missing_output_directory(tmp_path):
    """-o under a not-yet-created directory tree must not crash.

    No --generate-report involved: the plain TSV write used to raise
    FileNotFoundError because nothing created -o's parent directory.
    """
    gene_list = _write_query_genes(tmp_path / "query.txt")
    gene_sets = _write_gene_sets(tmp_path / "gene_sets.json")
    output = tmp_path / "nested" / "sub" / "results.tsv"
    assert not output.parent.exists()

    runner = CliRunner()
    result = runner.invoke(
        overlap_test,
        ["-g", str(gene_list), "-s", str(gene_sets), "-o", str(output)],
    )

    assert result.exit_code == 0, result.output
    assert output.exists()


def test_overlap_cli_json_output_with_report_does_not_crash(tmp_path):
    """-o results.json --generate-report must not crash.

    The HTML report reader (_load_results) always parses its input with
    ``pd.read_csv(path, sep="\\t")``; feeding it the literal JSON output
    path misread JSON as TSV (or, depending on --output-format, pointed at
    a file that --output-format never wrote at that literal path).
    """
    gene_list = _write_query_genes(tmp_path / "query.txt")
    gene_sets = _write_gene_sets(tmp_path / "gene_sets.json")
    output = tmp_path / "results.json"

    runner = CliRunner()
    result = runner.invoke(
        overlap_test,
        [
            "-g",
            str(gene_list),
            "-s",
            str(gene_sets),
            "-o",
            str(output),
            "--output-format",
            "json",
            "--generate-report",
        ],
    )

    assert result.exit_code == 0, result.output
    report_path = tmp_path / "enrichex_overlap_report.html"
    assert report_path.exists()
