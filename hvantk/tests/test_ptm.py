"""Tests for hvantk.ptm module.

Covers:
1. Coordinate mapper correctness (TP53 known positions from Phase 0/0.5)
2. Transcript resolution cascade (3-strategy)
3. map_ptm_sites pipeline round-trip (TSV → mapped TSV → verify output)
4. CLI wiring (ptm group registers and renders help)
"""

import csv
import os

import pytest
from click.testing import CliRunner

from hvantk.algorithms.ptm.mapper import (
    CodonMapping,
    GTFData,
    map_residue_to_genomic,
    map_protein_sites,
    resolve_transcript,
)
from hvantk.algorithms.ptm.pipeline import PTMBuildConfig, map_ptm_sites
from hvantk.core.ptm_constants import PTM_OUTPUT_COLUMNS


# ---------- Fixtures ----------

# Real TP53 CDS exons from Ensembl GTF release 113 (ENST00000269305, chr17, minus strand)
_TP53_CDS = [
    ("17", 7669612, 7669690, "-", 1),
    ("17", 7670609, 7670715, "-", 0),
    ("17", 7673535, 7673608, "-", 2),
    ("17", 7673701, 7673837, "-", 1),
    ("17", 7674181, 7674290, "-", 0),
    ("17", 7674859, 7674971, "-", 2),
    ("17", 7675053, 7675236, "-", 0),
    ("17", 7675994, 7676272, "-", 0),
    ("17", 7676382, 7676403, "-", 1),
    ("17", 7676521, 7676594, "-", 0),
]


@pytest.fixture
def tp53_cds_lookup():
    return {"ENST00000269305": _TP53_CDS}


@pytest.fixture
def gtf_data(tp53_cds_lookup):
    return GTFData(
        cds_by_transcript=tp53_cds_lookup,
        mane_transcripts={"ENST00000269305"},
        gene_to_mane={"TP53": "ENST00000269305"},
        transcript_to_gene={"ENST00000269305": "TP53"},
    )


# ---------- Test 1: Mapper correctness ----------

# Ground truth from Phase 0 pilot and Phase 0.5 validation (Ensembl REST API)
_TP53_GROUND_TRUTH = [
    (6, "17", 7676577, 7676579),
    (15, "17", 7676550, 7676552),
    (46, "17", 7676231, 7676233),
    (150, "17", 7675162, 7675164),
    (315, "17", 7673583, 7673585),  # Ser315 — key PTM site
    (382, "17", 7669645, 7669647),
]


@pytest.mark.parametrize("residue_pos,exp_chrom,exp_start,exp_end", _TP53_GROUND_TRUTH)
def test_mapper_tp53_positions(
    tp53_cds_lookup, residue_pos, exp_chrom, exp_start, exp_end
):
    """Mapper reproduces Phase 0/0.5 ground truth for TP53 positions."""
    result = map_residue_to_genomic("ENST00000269305", residue_pos, tp53_cds_lookup)
    assert result is not None
    assert result.chrom == exp_chrom
    assert result.codon_start == exp_start
    assert result.codon_end == exp_end
    assert result.strand == "-"


def test_mapper_out_of_bounds(tp53_cds_lookup):
    """Positions beyond CDS length return None (isoform mismatch case)."""
    # TP53 has 393 codons, position 999 is out of bounds
    assert map_residue_to_genomic("ENST00000269305", 999, tp53_cds_lookup) is None


def test_mapper_unknown_transcript(tp53_cds_lookup):
    """Unknown transcript returns None."""
    assert map_residue_to_genomic("ENST_FAKE", 1, tp53_cds_lookup) is None


# ---------- Test 2: Transcript resolution cascade ----------


def test_resolve_transcript_cascade(gtf_data):
    """3-strategy resolution cascade works as documented."""
    # Strategy 1: MANE Select xref (versioned ID — strip suffix)
    enst, method = resolve_transcript([{"id": "ENST00000269305.8"}], "TP53", gtf_data)
    assert enst == "ENST00000269305"
    assert method == "xref_mane"

    # Strategy 2: any xref (when no MANE transcripts in GTF)
    gtf_no_mane = GTFData(
        cds_by_transcript=gtf_data.cds_by_transcript,
        mane_transcripts=set(),
        gene_to_mane={},
        transcript_to_gene=gtf_data.transcript_to_gene,
    )
    enst2, method2 = resolve_transcript(
        [{"id": "ENST00000269305.8"}], "TP53", gtf_no_mane
    )
    assert enst2 == "ENST00000269305"
    assert method2 == "xref_any"

    # Strategy 3: gene name → MANE fallback (no xrefs available)
    enst3, method3 = resolve_transcript([], "TP53", gtf_data)
    assert enst3 == "ENST00000269305"
    assert method3 == "gene_mane"

    # Unresolved
    enst4, method4 = resolve_transcript([], "NONEXISTENT", gtf_data)
    assert enst4 is None
    assert method4 == "unresolved"


# ---------- Test 3: map_ptm_sites round-trip ----------


def test_map_ptm_sites_roundtrip(tmp_path, gtf_data):
    """map_ptm_sites reads a UniProt TSV and writes a mapped TSV with correct output."""
    # Create a minimal input TSV (UniProt format)
    input_tsv = tmp_path / "input.tsv"
    with open(input_tsv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "accession",
                "gene_symbol",
                "position",
                "description",
                "amino_acid",
                "ensembl_xrefs",
                "sequence_length",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "accession": "P04637",
                "gene_symbol": "TP53",
                "position": "315",
                "description": "Phosphoserine",
                "amino_acid": "S",
                "ensembl_xrefs": "ENST00000269305.8",
                "sequence_length": "393",
            }
        )
        writer.writerow(
            {
                "accession": "P04637",
                "gene_symbol": "TP53",
                "position": "6",
                "description": "Phosphoserine",
                "amino_acid": "S",
                "ensembl_xrefs": "ENST00000269305.8",
                "sequence_length": "393",
            }
        )

    output_bgz = tmp_path / "mapped.tsv.bgz"
    result = map_ptm_sites(str(input_tsv), gtf_data, str(output_bgz))

    assert result.n_total == 2
    assert result.n_mapped == 2
    assert result.n_failed == 0
    assert result.resolution_counts["xref_mane"] == 1

    # Read back and verify (BGZF is gzip-compatible)
    import gzip
    with gzip.open(output_bgz, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 2
    # Ser315 row
    ser315 = next(r for r in rows if r["residue_pos"] == "315")
    assert ser315["chrom"] == "17"
    assert ser315["codon_start"] == "7673583"
    assert ser315["codon_end"] == "7673585"
    assert ser315["ptm_category"] == "phosphorylation"
    assert ser315["gene_symbol"] == "TP53"
    assert ser315["source_db"] == "UniProt"
    assert ser315["evidence_type"] == "curated"
    assert ser315["n_observations"] == "0"

    # Check all expected columns are present
    assert set(reader.fieldnames) == set(PTM_OUTPUT_COLUMNS)


# ---------- Test 4: CLI registration ----------


def test_ptm_cli_help():
    """PTM command group is registered and renders help."""
    from hvantk.tools.ptm.ptm_cli import ptm_group

    runner = CliRunner()
    result = runner.invoke(ptm_group, ["--help"])
    assert result.exit_code == 0
    assert "build" in result.output
    assert "annotate" in result.output
    assert "landscape" in result.output

    # Build subcommand help
    result2 = runner.invoke(ptm_group, ["build", "--help"])
    assert result2.exit_code == 0
    assert "--output-dir" in result2.output
    assert "--output-ht" in result2.output
    assert "--gtf-path" in result2.output


def test_ptm_build_config_validation():
    """PTMBuildConfig.validate catches missing required fields."""
    config = PTMBuildConfig()
    errors = config.validate()
    assert len(errors) == 2
    assert any("output_dir" in e for e in errors)
    assert any("output_ht" in e for e in errors)

    config2 = PTMBuildConfig(output_dir="/tmp/test", output_ht="/tmp/test/out.ht")
    assert config2.validate() == []
