"""
Tests for the EnrichEx `overlap` CLI command.

Overlap enrichment uses scipy, not Hail, so these exercise the Click
command end-to-end (CliRunner) with no Hail session needed.
"""

from pathlib import Path

from click.testing import CliRunner

from hvantk.core.utils.gene_sets import GeneSet, GeneSetCollection
from hvantk.tools.enrichex.overlap_cli import overlap_test


def _invoke_overlap(gene_list, gene_sets, output, extra_args=()):
    runner = CliRunner()
    return runner.invoke(
        overlap_test,
        [
            "-g",
            str(gene_list),
            "-s",
            str(gene_sets),
            "-o",
            str(output),
            "--generate-report",
            *extra_args,
        ],
    )


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


def test_overlap_cli_rerun_report_reflects_second_run_not_a_stale_one(tmp_path):
    """Rerunning the same -o with a different gene list must report the
    SECOND run's results, not the first.

    The primary-output block treats a suffix-less -o as a file prefix
    (writes results.tsv); the report block used to re-derive the output
    location with different rules (treating "results" as a directory) and
    staged its own TSV copy there only `if not results_path.exists()` -- so
    a rerun into the same -o kept building the HTML report from the first
    run's results.
    """
    gene_sets = _write_gene_sets(tmp_path / "gene_sets.json")
    output = tmp_path / "results"  # no suffix -- a file prefix, not a directory
    report_path = tmp_path / "enrichex_overlap_report.html"

    gene_list_1 = _write_query_genes(tmp_path / "query1.txt")  # A,B,C -> n_overlap=3
    result1 = _invoke_overlap(gene_list_1, gene_sets, output)
    assert result1.exit_code == 0, result1.output
    assert "<td>3</td>" in report_path.read_text()

    gene_list_2 = tmp_path / "query2.txt"
    gene_list_2.write_text("A\n")  # A only -> n_overlap=1
    result2 = _invoke_overlap(gene_list_2, gene_sets, output)
    assert result2.exit_code == 0, result2.output

    html = report_path.read_text()
    assert "<td>1</td>" in html
    assert "<td>3</td>" not in html
