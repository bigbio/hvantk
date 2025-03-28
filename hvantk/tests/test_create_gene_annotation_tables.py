from pathlib import Path

from hvantk.utils.make_tables import (
    create_gnomad_constraint_gene_metrics_tb,
    create_interactome_tb,
    create_clinvar_tb,
    create_gevir_tb,
    create_ensembl_gene_tb,
)

# get the root directory of the project
PROJECT_DIR = Path(__file__).parent.parent


def test_create_gnomad_constraint_gene_metrics_tb():
    input_path = (
        PROJECT_DIR
        / "testdata/raw/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.chr20.tsv.bgz"
    )
    output_path = (
        PROJECT_DIR
        / "testdata/hail_tables/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.chr20.ht"
    )
    fields = ["oe_syn_upper", "oe_mis_upper", "oe_lof_upper"]  # fields to select
    overwrite = True
    export_tsv = True

    tb = create_gnomad_constraint_gene_metrics_tb(
        input_path=input_path,
        output_path=output_path,
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )

    tb.show()

    # check the number of rows
    assert tb.count() == 542
    # check that the .SUCCESS file exists
    success_file = Path(output_path) / "_SUCCESS"
    assert success_file.exists()


def test_create_interactome_tb():
    input_path = (
        PROJECT_DIR
        / "testdata/raw/interactome/Interactome_INSIDER_hg38_stripped.chr20.bed.bgz"
    )
    output_path = (
        PROJECT_DIR
        / "testdata/hail_tables/interactome/Interactome_INSIDER_hg38_stripped.chr20.ht"
    )
    overwrite = True
    export_tsv = True
    reference_genome = "GRCh38"

    tb = create_interactome_tb(
        input_path=input_path,
        output_path=output_path,
        overwrite=overwrite,
        export_tsv=export_tsv,
        reference_genome=reference_genome,
    )

    tb.show()

    # check the number of rows
    assert tb.count() == 32781
    # check that the .SUCCESS file exists
    success_file = Path(output_path) / "_SUCCESS"
    assert success_file.exists()


def test_create_clinvar_tb():
    input_path = PROJECT_DIR / "testdata/raw/clinvar/clinvar_20220403_chr20.vcf.bgz"
    output_path = PROJECT_DIR / "testdata/hail_tables/clinvar/clinvar_20220403_chr20.ht"
    overwrite = True
    export_tsv = True

    tb = create_clinvar_tb(
        input_path=input_path,
        output_path=output_path,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )

    tb.show()

    # check the number of rows
    assert tb.count() == 24554
    # check that the .SUCCESS file exists
    success_file = Path(output_path) / "_SUCCESS"
    assert success_file.exists()


def test_create_gevir_tb():
    input_path = PROJECT_DIR / "testdata/raw/gevir/gevir_metrics_pmid31873297.tsv.bgz"
    output_path = (
        PROJECT_DIR / "testdata/hail_tables/gevir/gevir_metrics_pmid31873297.ht"
    )
    overwrite = True
    export_tsv = True

    tb = create_gevir_tb(
        input_path=input_path,
        output_path=output_path,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )

    tb.show()

    # check the number of rows
    assert tb.count() == 19361
    # check that the .SUCCESS file exists
    success_file = Path(output_path) / "_SUCCESS"
    assert success_file.exists()


def test_create_ensembl_gene_tb():
    input_path = PROJECT_DIR / "testdata/raw/ensembl/ensembl_gene_biomart.tsv.bgz"
    output_path = PROJECT_DIR / "testdata/hail_tables/ensembl/ensembl_gene_biomart.ht"
    fields = None  # select all fields
    canonical = True
    overwrite = True

    tb = create_ensembl_gene_tb(
        input_path=input_path,
        output_path=output_path,
        fields=fields,
        canonical=canonical,
        overwrite=overwrite,
    )

    tb.show()

    # check the number of rows
    assert tb.count() == 86402
    # check that the .SUCCESS file exists
    success_file = Path(output_path) / "_SUCCESS"
    assert success_file.exists()
