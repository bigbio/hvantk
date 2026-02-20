"""
Hail integration tests for build_1k_genome_mt.

These tests create minimal synthetic bgzipped VCF files and run the full
MatrixTable build pipeline.  They require Hail/Spark and bgzip/tabix.
"""

import gzip
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

from hvantk.tables.genome_builders import build_1k_genome_mt

pytestmark = [pytest.mark.hail, pytest.mark.slow]

TEST_DIR = Path(__file__).parent / "testdata"
TMP_DIR = Path(__file__).parent / "tmp" / "genome_builders"

# Minimal VCF header + variants for chr1 (GRCh38 coordinates)
_VCF_TEMPLATE = textwrap.dedent("""\
    ##fileformat=VCFv4.1
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##contig=<ID=chr1,length=248956422>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\tNA12879
    chr1\t10177\t.\tA\tAC\t.\tPASS\t.\tGT\t0/1\t0/0
    chr1\t10352\t.\tT\tTA\t.\tPASS\t.\tGT\t1/1\t0/1
    chr1\t11008\t.\tC\tG\t.\tPASS\t.\tGT\t0/0\t0/1
""")

_PHENO_TSV = textwrap.dedent("""\
    sample_id\tpopulation\tsuper_population
    NA12878\tCEU\tEUR
    NA12879\tCEU\tEUR
""")


def _make_vcf(vcf_dir: Path, chrom: str, content: str) -> Path:
    """Write VCF content to a bgzipped file and index it with tabix."""
    vcf_path = vcf_dir / f"1kg_{chrom}.vcf.gz"
    # Write plain text to a tmp file then bgzip it
    tmp_vcf = vcf_dir / f"_{chrom}.vcf"
    tmp_vcf.write_text(content)
    subprocess.run(
        ["bgzip", "-f", str(tmp_vcf)],
        check=True,
        capture_output=True,
    )
    bgz_path = vcf_dir / f"_{chrom}.vcf.gz"
    bgz_path.rename(vcf_path)
    subprocess.run(
        ["tabix", "-p", "vcf", str(vcf_path)],
        check=True,
        capture_output=True,
    )
    return vcf_path


@pytest.fixture(autouse=True)
def setup_teardown():
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    yield
    if TMP_DIR.exists():
        shutil.rmtree(TMP_DIR)


@pytest.fixture()
def vcf_dir():
    d = TMP_DIR / "vcfs"
    d.mkdir(parents=True, exist_ok=True)
    _make_vcf(d, "chr1", _VCF_TEMPLATE)
    return d


def test_build_basic_mt(vcf_dir):
    """End-to-end: import a single-chromosome VCF and verify the MatrixTable."""
    output_mt = str(TMP_DIR / "out.mt")
    mt = build_1k_genome_mt(
        input_vcfs=str(vcf_dir),
        output_mt=output_mt,
        overwrite=True,
    )
    assert mt.count_rows() == 3
    assert mt.count_cols() == 2
    # Column key should be 's' (sample ID)
    assert "s" in mt.col_key.dtype


def test_build_with_phenotype(vcf_dir, tmp_path):
    """Phenotype TSV is joined to column annotations."""
    pheno_file = tmp_path / "pheno.tsv"
    pheno_file.write_text(_PHENO_TSV)
    output_mt = str(TMP_DIR / "out_pheno.mt")

    mt = build_1k_genome_mt(
        input_vcfs=str(vcf_dir),
        output_mt=output_mt,
        phenotype=str(pheno_file),
        overwrite=True,
    )
    assert "phenotype" in mt.col.dtype
    col_fields = list(mt.col.dtype["phenotype"])
    assert "population" in col_fields
    assert "super_population" in col_fields


def test_build_chromosome_filter(tmp_path):
    """Only the requested chromosomes are imported."""
    vcf_dir = TMP_DIR / "vcfs_two"
    vcf_dir.mkdir(parents=True, exist_ok=True)

    # chr1 has 3 variants; chr2 has 1 variant
    _make_vcf(vcf_dir, "chr1", _VCF_TEMPLATE)
    chr2_content = textwrap.dedent("""\
        ##fileformat=VCFv4.1
        ##FILTER=<ID=PASS,Description="All filters passed">
        ##contig=<ID=chr2,length=242193529>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\tNA12879
        chr2\t10000\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t0/0
    """)
    _make_vcf(vcf_dir, "chr2", chr2_content)

    output_mt = str(TMP_DIR / "out_filtered.mt")
    mt = build_1k_genome_mt(
        input_vcfs=str(vcf_dir),
        output_mt=output_mt,
        chromosomes=["chr1"],
        overwrite=True,
    )
    # Only chr1 variants should be present
    assert mt.count_rows() == 3