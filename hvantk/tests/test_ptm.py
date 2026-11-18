"""Tests for hvantk.ptm module.

Covers:
1. Coordinate mapper correctness (TP53 known positions from Ensembl REST lookups)
2. Transcript resolution cascade (3-strategy)
3. map_ptm_sites pipeline round-trip (TSV → mapped TSV → verify output), and the
   isoform guard: a UniProt site maps only through a transcript encoding its isoform
4. CLI wiring (ptm group registers and renders help)
5. Build core: records the sources it mapped and needs no output_ht
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
from hvantk.algorithms.ptm.constants import PTM_OUTPUT_COLUMNS
from hvantk.skills.uniprot_ptm.shared.datasets import (
    _TSV_COLUMNS as UNIPROT_TSV_COLUMNS,
    _extract_ptm_records,
)


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

# Ground truth from Ensembl REST API lookups
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
    """Mapper reproduces Ensembl REST ground truth for TP53 positions."""
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
                "isoform",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for position in ("315", "6", "25"):
            writer.writerow(
                {
                    "accession": "P04637",
                    "gene_symbol": "TP53",
                    "position": position,
                    "description": "Phosphoserine",
                    "amino_acid": "S",
                    "ensembl_xrefs": "ENST00000269305.8",
                    "sequence_length": "393",
                    "isoform": "",
                }
            )

    output_bgz = tmp_path / "mapped.tsv.bgz"
    result = map_ptm_sites(str(input_tsv), gtf_data, str(output_bgz))

    assert result.n_total == 3
    assert result.n_mapped == 3
    assert result.n_failed == 0
    assert result.resolution_counts["xref_mane"] == 1

    # Read back and verify (BGZF is gzip-compatible)
    import gzip

    with gzip.open(output_bgz, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    assert len(rows) == 3
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
    assert [ser315[f"codon_pos{i}"] for i in (1, 2, 3)] == [
        "7673585",
        "7673584",
        "7673583",
    ]

    # Residue 25's codon is split by an intron: two bases end one exon, the third
    # starts the next. codon_start/codon_end span the intron; the bases do not.
    split = next(r for r in rows if r["residue_pos"] == "25")
    assert (split["codon_start"], split["codon_end"]) == ("7676403", "7676522")
    assert [split[f"codon_pos{i}"] for i in (1, 2, 3)] == [
        "7676522",
        "7676521",
        "7676403",
    ]

    # Check all expected columns are present
    assert set(reader.fieldnames) == set(PTM_OUTPUT_COLUMNS)


def test_uniprot_sites_map_only_through_their_isoforms_transcripts(tmp_path):
    """A UniProt position is a residue of one isoform, so it must be placed on a
    transcript encoding that isoform (#486). The MANE Select transcript here encodes
    another isoform, whose residue 20 sits on a different codon: a canonical site placed
    there would land on the wrong codon without error. A site on an isoform no
    transcript encodes is dropped, never moved to the gene's MANE transcript.
    """
    # Canonical X-1: one exon, residue 20 = bases 1057-1059. X-3 (MANE): an exon of
    # 10 codons, a 30-base intron, then the rest, so its residue 20 = bases 1087-1089.
    gtf = GTFData(
        cds_by_transcript={
            "ENST_CAN": [("1", 1000, 1299, "+", 0)],
            "ENST_MANE": [("1", 1000, 1029, "+", 0), ("1", 1060, 1299, "+", 0)],
        },
        mane_transcripts={"ENST_MANE"},
        gene_to_mane={"GENEX": "ENST_MANE"},
        transcript_to_gene={"ENST_CAN": "GENEX", "ENST_MANE": "GENEX"},
    )

    def site(position, sequence=None):
        location = {"start": {"value": position}, "end": {"value": position}}
        if sequence:
            location["sequence"] = sequence
        return {
            "type": "Modified residue",
            "location": location,
            "description": "Phosphoserine",
        }

    entry = {
        "primaryAccession": "X",
        "genes": [{"geneName": {"value": "GENEX"}}],
        "sequence": {"value": "S" * 100, "length": 100},
        "comments": [
            {
                "commentType": "ALTERNATIVE PRODUCTS",
                "isoforms": [
                    {"isoformIds": ["X-1"], "isoformSequenceStatus": "Displayed"},
                    {"isoformIds": ["X-2"], "isoformSequenceStatus": "Described"},
                    {"isoformIds": ["X-3"], "isoformSequenceStatus": "Described"},
                ],
            }
        ],
        "uniProtKBCrossReferences": [
            {"database": "Ensembl", "id": "ENST_MANE.1", "isoformId": "X-3"},
            {"database": "Ensembl", "id": "ENST_CAN.1", "isoformId": "X-1"},
        ],
        "features": [site(20), site(5, "X-2"), site(20, "X-3")],
    }

    input_tsv = tmp_path / "uniprot.tsv"
    with open(input_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=UNIPROT_TSV_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(_extract_ptm_records(entry))

    result = map_ptm_sites(str(input_tsv), gtf, str(tmp_path / "mapped.tsv.bgz"))

    import gzip

    with gzip.open(result.mapped_tsv_path, "rt") as f:
        rows = {r["uniprot_id"]: r for r in csv.DictReader(f, delimiter="\t")}

    assert set(rows) == {"X", "X-3"}
    assert (rows["X"]["codon_start"], rows["X"]["amino_acid"]) == ("1057", "S")
    # Numbered on X-3: placed on X-3's transcript; the canonical residue is not its own
    assert (rows["X-3"]["codon_start"], rows["X-3"]["amino_acid"]) == ("1087", "")
    assert result.failure_counts == {"unresolved": 1}


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
    assert config.validate() == ["output_dir is required"]

    # output_ht is not required here: only hvantk.tools.ptm.pipeline writes the table.
    config2 = PTMBuildConfig(output_dir="/tmp/test")
    assert config2.validate() == []


# ---------- Test 5: build core ----------


def test_build_core_records_the_sources_it_mapped(tmp_path, gtf_data, monkeypatch):
    """`source_counts` gives the sites each source contributed, UniProt first.

    Runs without Hail, network or output_ht: the GTF parse is stubbed with the TP53
    fixture.
    """
    from hvantk.algorithms.ptm import pipeline as ptm_pipeline

    monkeypatch.setattr(ptm_pipeline, "parse_ensembl_gtf", lambda path: gtf_data)
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text("")
    row = {
        "accession": "P04637",
        "gene_symbol": "TP53",
        "position": "315",
        "description": "Phosphoserine",
        "amino_acid": "S",
        "ensembl_xrefs": "ENST00000269305.8",
        "isoform": "",
    }
    uniprot = tmp_path / "uniprot.tsv"
    peptideatlas = tmp_path / "peptideatlas.tsv"
    for path in (uniprot, peptideatlas):
        path.write_text("\t".join(row) + "\n" + "\t".join(row.values()) + "\n")

    result = ptm_pipeline.ptm_build_pipeline_core(
        PTMBuildConfig(
            output_dir=str(tmp_path / "out"),
            gtf_path=str(gtf),
            ptm_tsv=str(uniprot),
            peptideatlas_tsv=str(peptideatlas),
        )
    )

    assert result.source_counts == {"UniProt": 1, "PeptideAtlas": 1}
    assert result.n_mapped == 2
    assert result.mapped_tsv_path.endswith("ptm_sites_combined.tsv.bgz")

    # UniProt alone is the default `ptm build` run; it skips the concat step.
    alone = ptm_pipeline.ptm_build_pipeline_core(
        PTMBuildConfig(
            output_dir=str(tmp_path / "alone"),
            gtf_path=str(gtf),
            ptm_tsv=str(uniprot),
        )
    )
    assert alone.source_counts == {"UniProt": 1}
    assert alone.mapped_tsv_path.endswith("ptm_sites_mapped.tsv.bgz")


# ---------- Test 6: `ptm test`'s reader accepts compressed input ----------


def test_read_variants_table_accepts_bgz_compressed_tsv(tmp_path):
    """``ptm test``'s reader (``read_variants_tabular_file``, shared with
    ``ptm constraint``) relied on pandas' compression='infer' alone, which
    does not recognise the '.bgz' extension: pandas fell back to reading it
    as plain text and choked on the first gzip byte with UnicodeDecodeError.
    """
    import gzip

    from hvantk.algorithms.ptm.constraint import read_variants_tabular_file

    path = tmp_path / "variants.tsv.bgz"
    with gzip.open(path, "wt") as f:
        f.write("gene\taf\tis_ptm\tstratum\n")
        f.write("G1\t0.1\tTrue\tA\n")
        f.write("G2\t0.2\tFalse\tA\n")

    df = read_variants_tabular_file(str(path))
    assert list(df.columns) == ["gene", "af", "is_ptm", "stratum"]
    assert len(df) == 2
