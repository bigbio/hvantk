"""Per-gene structural summary parsed from an Ensembl GTF."""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from hvantk.skills.ensembl_gene.structure.parse import parse_gtf_structure

FIXTURE = str(
    Path("hvantk/skills/ensembl_gene/structure/tests/testdata/raw")
    / "ensembl-structure"
    / "mini.gtf"
)


@pytest.fixture(scope="module")
def df():
    return parse_gtf_structure(FIXTURE)


def test_one_row_per_gene_sorted_by_gene_id(df):
    assert list(df.gene_id) == [
        "ENSG00000000001",
        "ENSG00000000002",
        "ENSG00000000003",
    ]


def test_mane_wins_over_a_longer_non_mane_transcript(df):
    """ENST...001 is MANE with 300 bp CDS; ENST...002 is non-MANE with 1000 bp."""
    row = df[df.gene_id == "ENSG00000000001"].iloc[0]
    assert row.mane_select == "ENST00000000001"
    assert row.cds_transcript == "ENST00000000001"
    assert row.cds_length == 300  # (1099-1000+1) + (2199-2000+1)
    assert row.n_coding_exons == 2
    assert row.n_transcripts == 2


def test_longest_cds_wins_when_there_is_no_mane(df):
    """ENST...003 is 50 bp, ENST...004 is 300 bp, neither is MANE."""
    row = df[df.gene_id == "ENSG00000000002"].iloc[0]
    assert row.mane_select == ""
    assert row.cds_transcript == "ENST00000000004"
    assert row.cds_length == 300  # 1299-1000+1
    assert row.n_coding_exons == 1


def test_non_coding_gene_has_zero_cds(df):
    row = df[df.gene_id == "ENSG00000000003"].iloc[0]
    assert row.gene_biotype == "lncRNA"
    assert row.cds_transcript == ""
    assert row.cds_length == 0
    assert row.n_coding_exons == 0
    assert row.n_transcripts == 1


def test_version_suffixes_are_stripped(tmp_path):
    gtf = tmp_path / "v.gtf"
    gtf.write_text(
        '1\te\ttranscript\t1\t9\t.\t+\t.\tgene_id "ENSG00000000009.7"; '
        'transcript_id "ENST00000000009.3"; gene_biotype "protein_coding";\n'
        '1\te\tCDS\t1\t9\t.\t+\t0\tgene_id "ENSG00000000009.7"; '
        'transcript_id "ENST00000000009.3";\n'
    )
    out = parse_gtf_structure(str(gtf))
    assert out.iloc[0].gene_id == "ENSG00000000009"
    assert out.iloc[0].cds_transcript == "ENST00000000009"


def test_tie_break_deterministic_on_equal_cds_length(tmp_path):
    """When two non-MANE transcripts have equal CDS length, output must be deterministic.

    Without transcript ID in the tie-break key, max() would resolve ties by set iteration
    order (Python's per-process hash seed), causing non-deterministic winners across runs.
    This test verifies that output is identical across at least 8 subprocess runs with
    different PYTHONHASHSEED values, proving the tuple-key (length, transcript_id) is
    load-bearing for determinism.

    The winner should be the transcript with the larger ID (lexicographically), since
    (length, transcript_id) tuple comparison sorts on the second element to break ties.
    """
    gtf = tmp_path / "tie.gtf"
    gtf.write_text(
        '1\te\ttranscript\t1\t9\t.\t+\t.\tgene_id "ENSG00000000010"; '
        'transcript_id "ENST00000000010"; gene_biotype "protein_coding";\n'
        '1\te\tCDS\t1\t9\t.\t+\t0\tgene_id "ENSG00000000010"; '
        'transcript_id "ENST00000000010";\n'
        '1\te\ttranscript\t1\t9\t.\t+\t.\tgene_id "ENSG00000000010"; '
        'transcript_id "ENST00000000011"; gene_biotype "protein_coding";\n'
        '1\te\tCDS\t1\t9\t.\t+\t0\tgene_id "ENSG00000000010"; '
        'transcript_id "ENST00000000011";\n'
    )

    # Snippet to run in subprocess: parse the GTF and print the winner cds_transcript.
    snippet = f"""
import sys
sys.path.insert(0, "/fs/dss/home/heto4575/projects/pyvatk")
from hvantk.skills.ensembl_gene.structure.parse import parse_gtf_structure
df = parse_gtf_structure({str(gtf)!r})
print(df.iloc[0]["cds_transcript"])
"""

    # Run the parser 8 times, each with a different PYTHONHASHSEED.
    winners = []
    for seed in range(8):
        env = {**os.environ, "PYTHONHASHSEED": str(seed)}
        result = subprocess.run(
            [sys.executable, "-c", snippet],
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        winner = result.stdout.strip()
        winners.append(winner)

    # All runs must produce the identical winner; if they differ, the determinism is broken.
    assert len(set(winners)) == 1, (
        f"Determinism broken: got different winners across seeds: {winners}. "
        "The tuple-key (cds_bp[t], t) is not being used."
    )

    # Verify the winner is the larger transcript ID (the expected tie-break result).
    expected_winner = "ENST00000000011"
    assert winners[0] == expected_winner
