"""Per-plugin conformance tests for the plugin builder contract.

Each test:
  1. Looks up the dataset's DatasetSpec from the live registry.
  2. Prepares minimal parsed input (using the plugin's own test fixtures).
  3. Calls run_builder_for_spec() with the spec.
  4. Verifies the saved artifact loads back via core/io with the expected type
     and provenance fields.

Each plugin has its own test.
"""

from __future__ import annotations

import csv
import gzip
import shutil
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from hvantk.core import io as core_io
from hvantk.core.models import AnnotationTable, ExpressionMatrix, VariantMatrix
from hvantk.core.plugin import loader as plugin_loader
from hvantk.core.plugin.run_builder import run_builder_for_spec
from hvantk.tools.plugins.reprocess_cli import _expected_extensions, reprocess_cmd


# ---------- peptideatlas:phospho ----------


@pytest.fixture
def peptideatlas_phospho_parsed_tsv(tmp_path):
    """Write a minimal intermediate phospho TSV with the schema produced by parse_raw_dir."""
    cols = [
        "gene_symbol",
        "protein_accession",
        "site_position",
        "residue",
        "n_observations",
        "peptide_sequences",
        "source_db",
        "evidence_type",
    ]
    rows = [
        {
            "gene_symbol": "TP53",
            "protein_accession": "P04637",
            "site_position": "315",
            "residue": "S",
            "n_observations": "2",
            "peptide_sequences": "ABCDE;FGHIJ",
            "source_db": "peptideatlas",
            "evidence_type": "phospho",
        },
        {
            "gene_symbol": "BRCA1",
            "protein_accession": "P38398",
            "site_position": "988",
            "residue": "S",
            "n_observations": "1",
            "peptide_sequences": "KLMNO",
            "source_db": "peptideatlas",
            "evidence_type": "phospho",
        },
    ]
    tsv_path = tmp_path / "peptideatlas-phospho.tsv"
    with open(tsv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return tsv_path


def test_peptideatlas_phospho_round_trip(tmp_path, peptideatlas_phospho_parsed_tsv):
    # Force a fresh registry so any module-level cache from other tests is bypassed.
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("peptideatlas:phospho")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "peptideatlas-phospho-v1"
    assert spec.plugin_version

    out = tmp_path / "phospho.parquet"

    # Avoid network calls in the drift_probe by monkey-patching the spec's probe
    # to return a fixed fingerprint. (Test isolation.)
    original_probe = spec.drift_probe
    object.__setattr__(spec, "drift_probe", lambda: {"fingerprint": "sha256:test"})
    try:
        prov = run_builder_for_spec(
            spec,
            parsed_input=peptideatlas_phospho_parsed_tsv,
            output_path=out,
            plugin_version=spec.plugin_version,
        )
    finally:
        object.__setattr__(spec, "drift_probe", original_probe)

    assert prov.plugin == "peptideatlas"
    assert prov.dataset == "peptideatlas:phospho"
    assert prov.schema_id == "peptideatlas-phospho-v1"
    assert prov.source_fingerprint  # nonempty

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.provenance == prov
    assert loaded.count() == 2


# ---------- clinvar:variants ----------


@pytest.mark.hail
def test_clinvar_variants_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("clinvar:variants")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "clinvar-variants-v1"
    assert spec.plugin_version

    # Use the bundled fixture VCF (chr20 subset)
    fixture = Path(
        "hvantk/skills/clinvar/tests/testdata/raw/clinvar/clinvar_20220403_chr20.vcf.bgz"
    )
    assert fixture.exists(), f"missing fixture: {fixture}"

    # Avoid network calls in the drift_probe (clinvar probe hits NCBI)
    object.__setattr__(
        spec, "drift_probe", lambda: {"fingerprint": "sha256:test-clinvar"}
    )

    out = tmp_path / "variants.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "clinvar"
    assert prov.dataset == "clinvar:variants"
    assert prov.schema_id == "clinvar-variants-v1"
    assert prov.source_fingerprint == "sha256:test-clinvar"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.backend == "hail"
    assert loaded.provenance == prov
    # Sanity: fixture has > 0 variants
    assert loaded.count() > 0


# ---------- hgnc:lookup ----------


@pytest.mark.hail
def test_hgnc_lookup_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("hgnc:lookup")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "hgnc-lookup-v1"
    assert spec.plugin_version

    fixture = Path("hvantk/skills/hgnc/tests/testdata/raw/hgnc/hgnc_test_sample.tsv")
    assert fixture.exists(), f"missing fixture: {fixture}"

    object.__setattr__(spec, "drift_probe", lambda: {"fingerprint": "sha256:test-hgnc"})

    out = tmp_path / "lookup.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "hgnc"
    assert prov.dataset == "hgnc:lookup"
    assert prov.schema_id == "hgnc-lookup-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.backend == "hail"
    assert loaded.count() > 0


# ---------- gencc:submissions ----------


@pytest.mark.hail
def test_gencc_submissions_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("gencc:submissions")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "gencc-submissions-v1"
    assert spec.plugin_version

    fixture = Path("hvantk/skills/gencc/tests/testdata/raw/gencc/gencc_test_sample.tsv")
    assert fixture.exists()

    object.__setattr__(
        spec, "drift_probe", lambda: {"fingerprint": "sha256:test-gencc"}
    )

    out = tmp_path / "submissions.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "gencc"
    assert prov.dataset == "gencc:submissions"
    assert prov.schema_id == "gencc-submissions-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.backend == "hail"
    assert loaded.count() > 0


# ---------- clingen:gene-disease ----------


@pytest.mark.hail
def test_clingen_gene_disease_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("clingen:gene-disease")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "clingen-gene-disease-v1"

    fixture = Path(
        "hvantk/skills/clingen/tests/testdata/raw/clingen/clingen_test_sample.csv"
    )
    assert fixture.exists()

    object.__setattr__(
        spec, "drift_probe", lambda: {"fingerprint": "sha256:test-clingen"}
    )

    out = tmp_path / "gene-disease.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "clingen"
    assert prov.schema_id == "clingen-gene-disease-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.backend == "hail"
    assert loaded.count() > 0


# ---------- gwas-catalog:associations ----------


@pytest.mark.hail
def test_gwas_catalog_associations_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("gwas-catalog:associations")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "gwas-catalog-associations-v1"

    fixture = Path(
        "hvantk/skills/gwas_catalog/tests/testdata/raw/gwas-catalog/gwas-catalog-sample.tsv"
    )
    assert fixture.exists()

    object.__setattr__(spec, "drift_probe", lambda: {"fingerprint": "sha256:test-gwas"})

    out = tmp_path / "associations.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "gwas-catalog"
    assert prov.schema_id == "gwas-catalog-associations-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.backend == "hail"
    assert loaded.count() > 0


# ---------- gtex-eqtl:eqtls ----------


@pytest.mark.hail
def test_gtex_eqtl_eqtls_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("gtex-eqtl:eqtls")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "gtex-eqtl-eqtls-v1"

    fixture = Path(
        "hvantk/skills/gtex_eqtl/tests/testdata/raw/gtex-eqtl/Liver.v11.eQTLs.signif_pairs.parquet"
    )
    assert fixture.exists()

    object.__setattr__(spec, "drift_probe", lambda: {"fingerprint": "sha256:test-gtex"})

    out = tmp_path / "eqtls.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=str(fixture.parent),  # parquet importer expects DIRECTORY
        output_path=out,
        plugin_version=spec.plugin_version,
        source="gtex_v11",
        p_threshold=0,  # keep all rows (fixture has few rows)
    )

    assert prov.plugin == "gtex-eqtl"
    assert prov.schema_id == "gtex-eqtl-eqtls-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.backend == "hail"


# ---------- insider:variants ----------


@pytest.mark.hail
def test_insider_variants_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("insider:variants")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "insider-variants-v1"

    fixture = Path(
        "hvantk/skills/insider/variants/tests/testdata/raw/insider/"
        "insider_sample-synthetic.bed"
    )
    assert fixture.exists()

    object.__setattr__(
        spec, "drift_probe", lambda: {"fingerprint": "sha256:test-insider"}
    )

    out = tmp_path / "variants.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "insider"
    assert prov.schema_id == "insider-variants-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.backend == "hail"
    assert loaded.count() > 0


# ---------- msigdb:genesets ----------


@pytest.mark.hail
def test_msigdb_genesets_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("msigdb:genesets")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "msigdb-genesets-v1"

    fixture = Path("hvantk/skills/msigdb/tests/testdata/raw/msigdb/c2.cp-sample.gmt")
    assert fixture.exists()

    object.__setattr__(
        spec, "drift_probe", lambda: {"fingerprint": "sha256:test-msigdb"}
    )

    out = tmp_path / "genesets.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "msigdb"
    assert prov.schema_id == "msigdb-genesets-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.backend == "hail"
    assert loaded.count() > 0


# ---------- uniprot-ptm:sites ----------


@pytest.fixture
def uniprot_ptm_sites_tsv(tmp_path):
    """Write a minimal mapped PTM coordinates TSV with the columns the builder expects."""
    cols = [
        "chrom",
        "codon_start",
        "codon_end",
        "strand",
        "uniprot_id",
        "gene_symbol",
        "residue_pos",
        "amino_acid",
        "ptm_type",
        "ptm_category",
        "source_db",
        "evidence_type",
        "n_observations",
        "tissue_type",
    ]
    rows = [
        # Use chr17 (TP53) and chr13 (BRCA2) for plausibility
        [
            "17",
            "7676272",
            "7676274",
            "-",
            "P04637",
            "TP53",
            "315",
            "S",
            "phosphoserine",
            "phosphorylation",
            "uniprot",
            "experimental",
            "10",
            "",
        ],
        [
            "13",
            "32316461",
            "32316463",
            "+",
            "P51587",
            "BRCA2",
            "988",
            "S",
            "phosphoserine",
            "phosphorylation",
            "uniprot",
            "experimental",
            "5",
            "",
        ],
    ]
    tsv = tmp_path / "ptm_sites.tsv"
    with open(tsv, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    return tsv


@pytest.mark.hail
def test_uniprot_ptm_sites_round_trip(tmp_path, uniprot_ptm_sites_tsv):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("uniprot-ptm:sites")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "uniprot-ptm-sites-v1"

    object.__setattr__(spec, "drift_probe", lambda: {"fingerprint": "sha256:test-ptm"})

    out = tmp_path / "sites.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=uniprot_ptm_sites_tsv,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "uniprot-ptm"
    assert prov.schema_id == "uniprot-ptm-sites-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.backend == "hail"
    assert loaded.count() == 2


# ---------- expression-atlas:dataset ----------


@pytest.fixture
def expression_atlas_inputs(tmp_path):
    """Generate minimal expression matrix TSV + SDRF metadata files."""
    # Expression matrix: 2 genes x 3 samples
    expr_lines = [
        "Gene ID\tGene Name\tSample1\tSample2\tSample3",
        "ENSG00000141510\tTP53\t10.5\t20.3\t15.7",
        "ENSG00000139618\tBRCA2\t5.2\t8.1\t3.4",
    ]
    expr_path = tmp_path / "expression.tsv"
    expr_path.write_text("\n".join(expr_lines) + "\n")

    # SDRF: long-format with no header.
    # Columns: accession, unused, sample_id, column_type, column_name, column_value
    # _import_sdrf reads exactly 6 columns in this order.
    sdrf_lines = [
        "E-MTAB-0001\t\tSample1\tcharacteristic\torganism\tHomo sapiens",
        "E-MTAB-0001\t\tSample1\tfactor\tdisease\thealthy",
        "E-MTAB-0001\t\tSample2\tcharacteristic\torganism\tHomo sapiens",
        "E-MTAB-0001\t\tSample2\tfactor\tdisease\tcancer",
        "E-MTAB-0001\t\tSample3\tcharacteristic\torganism\tHomo sapiens",
        "E-MTAB-0001\t\tSample3\tfactor\tdisease\thealthy",
    ]
    sdrf_path = tmp_path / "metadata.sdrf.tsv"
    sdrf_path.write_text("\n".join(sdrf_lines) + "\n")

    return {"expression_matrix": expr_path, "sdrf": sdrf_path}


def test_expression_atlas_dataset_round_trip(tmp_path, expression_atlas_inputs):
    """First ExpressionMatrix conformance test — pandas-driven (no Hail dep)."""
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("expression-atlas:dataset")

    assert spec.artifact_type is ExpressionMatrix
    assert spec.schema_id == "expression-atlas-dataset-v1"

    object.__setattr__(
        spec, "drift_probe", lambda: {"fingerprint": "sha256:test-expr-atlas"}
    )

    out = tmp_path / "dataset.h5ad"
    prov = run_builder_for_spec(
        spec,
        parsed_input=expression_atlas_inputs,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "expression-atlas"
    assert prov.schema_id == "expression-atlas-dataset-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, ExpressionMatrix)
    assert loaded.n_vars == 2  # 2 genes
    assert loaded.n_obs == 3  # 3 samples


# ---------- cptac:expression ----------


@pytest.fixture
def cptac_expression_inputs(tmp_path):
    """Long-format CPTAC expression: rows are (sample, gene) pairs with values."""
    expr_lines = [
        "SampleID\tGeneID\tGene Name\tExpression",
        "S1\tENSG00000141510\tTP53\t10.5",
        "S1\tENSG00000139618\tBRCA2\t5.2",
        "S2\tENSG00000141510\tTP53\t20.3",
        "S2\tENSG00000139618\tBRCA2\t8.1",
    ]
    expr_path = tmp_path / "expression.tsv"
    expr_path.write_text("\n".join(expr_lines) + "\n")

    meta_lines = [
        "SampleID\ttissue\tdisease",
        "S1\tliver\tnormal",
        "S2\tliver\ttumor",
    ]
    meta_path = tmp_path / "metadata.tsv"
    meta_path.write_text("\n".join(meta_lines) + "\n")

    return {"expression": expr_path, "metadata": meta_path}


def test_cptac_expression_round_trip(tmp_path, cptac_expression_inputs):
    from hvantk.core.models import ExpressionMatrix

    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("cptac:expression")

    assert spec.artifact_type is ExpressionMatrix
    assert spec.schema_id == "cptac-expression-v1"

    object.__setattr__(
        spec, "drift_probe", lambda: {"fingerprint": "sha256:test-cptac-expr"}
    )

    out = tmp_path / "expression.h5ad"
    prov = run_builder_for_spec(
        spec,
        parsed_input=cptac_expression_inputs,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "cptac"
    assert prov.dataset == "cptac:expression"
    assert prov.schema_id == "cptac-expression-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, ExpressionMatrix)
    assert loaded.n_obs == 2  # 2 samples
    assert loaded.n_vars == 2  # 2 genes


# ---------- cptac:phospho ----------


@pytest.fixture
def cptac_phospho_inputs(tmp_path):
    """Wide-format CPTAC phospho: rows are sites, columns are samples + metadata."""
    # Site column header matches what CPTACPhosphoDataset.write_matrix_csv writes
    # (df.index.name = "Site"); the builder's site_id_col default is "Site" (#198).
    expr_lines = [
        "Site\tS1\tS2",
        "TP53_S315\t1.2\t3.4",
        "BRCA2_S988\t0.5\t2.1",
    ]
    expr_path = tmp_path / "phospho.tsv"
    expr_path.write_text("\n".join(expr_lines) + "\n")

    meta_lines = [
        "SampleID\ttissue\tdisease",
        "S1\tliver\tnormal",
        "S2\tliver\ttumor",
    ]
    meta_path = tmp_path / "metadata.tsv"
    meta_path.write_text("\n".join(meta_lines) + "\n")

    return {"expression": expr_path, "metadata": meta_path}


def test_cptac_phospho_round_trip(tmp_path, cptac_phospho_inputs):
    from hvantk.core.models import ExpressionMatrix

    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("cptac:phospho")

    assert spec.artifact_type is ExpressionMatrix
    assert spec.schema_id == "cptac-phospho-v1"

    object.__setattr__(
        spec, "drift_probe", lambda: {"fingerprint": "sha256:test-cptac-phospho"}
    )

    out = tmp_path / "phospho.h5ad"
    prov = run_builder_for_spec(
        spec,
        parsed_input=cptac_phospho_inputs,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "cptac"
    assert prov.dataset == "cptac:phospho"
    assert prov.schema_id == "cptac-phospho-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, ExpressionMatrix)
    assert loaded.n_obs == 2
    assert loaded.n_vars == 2


# ---------- ucsc-cellbrowser (3 datasets sharing one builder) ----------


@pytest.fixture
def ucsc_cellbrowser_inputs(tmp_path):
    """Minimal UCSC expression matrix (genes x cells) + metadata."""
    expr_lines = [
        "gene\tcell1\tcell2",
        "ENSG00000141510|TP53\t10.5\t20.3",
        "ENSG00000139618|BRCA2\t5.2\t8.1",
    ]
    expr_path = tmp_path / "expression.tsv"
    expr_path.write_text("\n".join(expr_lines) + "\n")

    meta_lines = [
        "cellId\ttissue\tcell_type",
        "cell1\tliver\thepatocyte",
        "cell2\tliver\tkupffer",
    ]
    meta_path = tmp_path / "metadata.tsv"
    meta_path.write_text("\n".join(meta_lines) + "\n")

    return {"expression_matrix": expr_path, "metadata": meta_path}


@pytest.mark.parametrize(
    "dataset_name,expected_schema",
    [
        ("default", "ucsc-cellbrowser-default-v1"),
        ("adult-ctx", "ucsc-cellbrowser-adult-ctx-v1"),
        ("dev-ctx", "ucsc-cellbrowser-dev-ctx-v1"),
    ],
)
def test_ucsc_cellbrowser_round_trip(
    tmp_path, ucsc_cellbrowser_inputs, dataset_name, expected_schema
):
    from hvantk.core.models import ExpressionMatrix

    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset(f"ucsc-cellbrowser:{dataset_name}")

    assert spec.artifact_type is ExpressionMatrix
    assert spec.schema_id == expected_schema

    object.__setattr__(
        spec, "drift_probe", lambda: {"fingerprint": f"sha256:test-ucsc-{dataset_name}"}
    )

    out = tmp_path / f"{dataset_name}.h5ad"
    prov = run_builder_for_spec(
        spec,
        parsed_input=ucsc_cellbrowser_inputs,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "ucsc-cellbrowser"
    assert prov.dataset == f"ucsc-cellbrowser:{dataset_name}"
    assert prov.schema_id == expected_schema

    loaded = core_io.load(out)
    assert isinstance(loaded, ExpressionMatrix)
    assert loaded.n_vars == 2  # 2 genes
    assert loaded.n_obs == 2  # 2 cells


# ---------- gevir:metrics, gnomad-metrics:metrics, dbnsfp:variants ----------


@pytest.mark.hail
def test_gevir_metrics_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("gevir:metrics")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "gevir-metrics-v1"

    fixture = Path(
        "hvantk/tests/testdata/raw/gevir/gevir_metrics_pmid31873297-synthetic.tsv.bgz"
    )
    assert fixture.exists()

    object.__setattr__(spec, "drift_probe", lambda: {"source_version": "test"})

    out = tmp_path / "metrics.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )
    assert prov.plugin == "gevir"
    assert prov.schema_id == "gevir-metrics-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.count() > 0


@pytest.mark.hail
def test_gnomad_metrics_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("gnomad-metrics:metrics")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "gnomad-metrics-v1"

    fixture = Path(
        "hvantk/tests/testdata/raw/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.chr20.tsv.bgz"
    )
    assert fixture.exists()

    object.__setattr__(spec, "drift_probe", lambda: {"source_version": "test"})

    out = tmp_path / "metrics.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )
    assert prov.plugin == "gnomad-metrics"
    assert prov.schema_id == "gnomad-metrics-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.count() > 0


@pytest.mark.hail
def test_dbnsfp_variants_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("dbnsfp:variants")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "dbnsfp-v1"

    fixture = Path("hvantk/tests/testdata/raw/dbnsfp/dbNSFP4_v49a_example_variants.bgz")
    assert fixture.exists()

    object.__setattr__(spec, "drift_probe", lambda: {"source_version": "test"})

    out = tmp_path / "variants.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )
    assert prov.plugin == "dbnsfp"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.count() > 0


# ---------- onek-genomes:variants ----------


@pytest.mark.hail
def test_onek_genomes_variants_round_trip(tmp_path):
    """onek-genomes:variants builds a VariantMatrix from per-chromosome 1KG VCFs."""
    from hvantk.core.models import VariantMatrix

    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("onek-genomes:variants")

    assert spec.artifact_type is VariantMatrix
    assert spec.schema_id == "onek-genomes-variants-v1"

    fixture = Path("hvantk/skills/onek_genomes/tests/testdata/raw/onek_genomes")
    assert fixture.exists()

    object.__setattr__(
        spec, "drift_probe", lambda: {"fingerprint": "sha256:test-onek-variants"}
    )

    out = tmp_path / "variants.mt"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "onek-genomes"
    assert prov.dataset == "onek-genomes:variants"
    assert prov.schema_id == "onek-genomes-variants-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, VariantMatrix)
    assert loaded.n_variants > 0
    assert loaded.n_samples > 0


# ---------- onek-genomes:samples ----------


@pytest.mark.hail
def test_onek_genomes_samples_round_trip(tmp_path):
    """onek-genomes:samples builds an AnnotationTable from the IGSR sample TSV."""
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("onek-genomes:samples")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "onek-genomes-samples-v1"

    fixture = Path("hvantk/skills/onek_genomes/tests/testdata/raw/onek_genomes")
    assert fixture.exists()

    object.__setattr__(
        spec, "drift_probe", lambda: {"fingerprint": "sha256:test-onek-samples"}
    )

    out = tmp_path / "samples.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "onek-genomes"
    assert prov.dataset == "onek-genomes:samples"
    assert prov.schema_id == "onek-genomes-samples-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.count() > 0


# ---------- reprocess CLI: no-parse datasets must resolve their raw directory ----------
#
# Regression test for the bug fixed by `resolve_raw_file` (hvantk/core/utils/raw_inputs.py):
# a dataset that declares no `lifecycle.parse` has nothing to hand its builder but the raw
# directory -- `hvantk reprocess` forwards `raw_dir` verbatim as `parsed_input`
# (reprocess_cli.py: "No parse stage declared: the builder consumes the raw dir
# directly"), and a bring-your-own dataset gets the identical directory from the user's
# own first `--raw-dir`. hgnc, clinvar, clingen, gencc, msigdb, dbnsfp, cosmic-cgc, and
# insider:variants's builders used to pass that directory straight to a reader that
# cannot read a directory: Hail's `FileNotFoundException: '<raw dir>' is a directory`
# from `hl.import_table`/`hl.import_vcf` (hgnc, clinvar, clingen, gencc, msigdb), or
# Python's `IsADirectoryError` from `open()`/`hfs.open()` (dbnsfp, cosmic-cgc,
# insider:variants); gevir's builder failed on any second file in the directory. A few
# tests above this point already pass a whole fixture directory
# straight to `run_builder_for_spec` (gtex-eqtl's parquet reader and both
# onek-genomes datasets each glob a directory themselves, with no resolve_raw_file-style
# helper at all) -- those were never exercising the bug this test guards, since their
# own builders already handled a directory correctly from the start.
#
# This test drives the real `hvantk reprocess` CLI instead, parametrized over every
# dataset the live plugin registry declares with no parse stage -- downloaded or
# bring-your-own alike -- so a newly added no-parse dataset is covered automatically
# rather than silently falling outside the test's scope.

_CPTAC_CANCER_TYPE_FOR_TEST = "brca"


def _install_copytree(fixture_dir: Path, raw_dir: Path) -> dict:
    """Default installer: copy the committed fixture directory verbatim.

    Correct for every no-parse dataset whose builder resolves a raw directory either by
    glob/pattern (hgnc, clinvar, clingen, gencc, gevir, gnomad-metrics, gwas-catalog,
    msigdb, dbnsfp, cosmic-cgc, insider:variants), by the same exact filename its own
    downloader writes (insider:interfaces, onek-genomes:samples), by globbing the
    directory itself with no resolver at all (alphagenome's `*.parquet` scan,
    onek-genomes:variants's `*.vcf.gz` scan), or via dict-shaped sub-paths matching the
    fixture's own names (expression-atlas:dataset, ucsc-cellbrowser:default) -- the
    fixture is already laid out the way the resolver expects, so nothing needs renaming.
    """
    shutil.copytree(fixture_dir, raw_dir)
    return {}


def _install_cptac(
    fixture_dir: Path, raw_dir: Path, *, prefix: str, matrix_name: str
) -> dict:
    """cptac:expression / cptac:phospho: rename the fixture's generic file names to the
    ``<prefix>-<cancer_type>-{matrix,metadata}.csv`` pair
    ``resolve_raw_dir_inputs`` (hvantk/skills/cptac/shared/datasets.py) requires when
    handed a raw directory instead of a ``{"expression"/"metadata"}`` mapping. The real
    downloader writes exactly that per-cancer-type pair; the committed fixture (shared
    with the direct-dict-input conformance test above) uses generic ``expression.tsv`` /
    ``phospho.tsv`` / ``metadata.tsv`` names instead, since it has no cancer type of its
    own. ``cancer_type`` is returned as a plugin-arg so the builder knows which pair to
    look for.
    """
    shutil.copytree(fixture_dir, raw_dir)
    ct = _CPTAC_CANCER_TYPE_FOR_TEST
    (raw_dir / matrix_name).rename(raw_dir / f"{prefix}-{ct}-matrix.csv")
    (raw_dir / "metadata.tsv").rename(raw_dir / f"{prefix}-{ct}-metadata.csv")
    return {"cancer_type": ct}


def _install_cptac_expression(fixture_dir: Path, raw_dir: Path) -> dict:
    return _install_cptac(
        fixture_dir, raw_dir, prefix="cptac-expression", matrix_name="expression.tsv"
    )


def _install_cptac_phospho(fixture_dir: Path, raw_dir: Path) -> dict:
    return _install_cptac(
        fixture_dir, raw_dir, prefix="cptac-phospho", matrix_name="phospho.tsv"
    )


def _install_ensembl_gene_structure(fixture_dir: Path, raw_dir: Path) -> dict:
    """ensembl-gene:structure: gzip-compress the fixture GTF under the exact pinned-
    release filename ``_resolve_gtf_path`` (hvantk/skills/ensembl_gene/structure/
    builder.py) looks for. That resolver joins a single exact name, not a glob -- and
    ``parse_gtf_structure`` picks its gzip-vs-plain reader off the ``.gz`` suffix, not
    the content -- so the fixture (a plain-text ``mini.gtf``, for human readability)
    must be both renamed AND actually gzip-compressed, not just renamed.
    """
    from hvantk.resources.ensembl_release import ENSEMBL_GTF_FILENAME

    raw_dir.mkdir(parents=True, exist_ok=True)
    gtf_bytes = (fixture_dir / "mini.gtf").read_bytes()
    with gzip.open(raw_dir / ENSEMBL_GTF_FILENAME, "wb") as fh:
        fh.write(gtf_bytes)
    return {}


def _install_ucsc_cellbrowser_ctx(fixture_dir: Path, raw_dir: Path) -> dict:
    """ucsc-cellbrowser:adult-ctx / :dev-ctx: rename the fixture's ``*-synthetic.tsv``
    files to the generic ``expression_matrix.tsv`` / ``metadata.tsv`` names
    ``_resolve_ucsc_inputs`` (hvantk/skills/ucsc_cellbrowser/builder.py) falls back to
    when neither file matches the per-dataset canonical name.
    ucsc-cellbrowser:default's fixture already uses those generic names directly (no
    installer needed); these two datasets' fixtures carry a ``-synthetic`` suffix
    instead, since they replaced real data found to have no clear licence (#426).
    """
    shutil.copytree(fixture_dir, raw_dir)
    (raw_dir / "expression_matrix-synthetic.tsv").rename(
        raw_dir / "expression_matrix.tsv"
    )
    (raw_dir / "metadata-synthetic.tsv").rename(raw_dir / "metadata.tsv")
    return {}


def _install_gtex_eqtl(fixture_dir: Path, raw_dir: Path) -> dict:
    """gtex-eqtl:eqtls: plain copy -- ``scan_tissue_files``
    (hvantk/core/utils/qtl_helpers.py) already globs a directory for parquet/TSV
    files itself, no resolver needed. ``p_threshold=0`` keeps every row: the
    builder's default (5e-8, genome-wide significance) would drop some or all of
    the fixture's few rows.
    """
    shutil.copytree(fixture_dir, raw_dir)
    return {"source": "gtex_v11", "p_threshold": 0}


def _install_pqtl(fixture_dir: Path, raw_dir: Path) -> dict:
    """pqtl:metrics: plain copy -- ``scan_tissue_files`` already globs a directory
    itself, no resolver needed. ``no_gene_map=true`` opts out of the
    gene-symbol-to-Ensembl-ID mapping that otherwise requires a ``gene_catalog``
    (``--plugin-arg hgnc_ht=<path>``), which this offline test does not build.
    """
    shutil.copytree(fixture_dir, raw_dir)
    return {"no_gene_map": "true"}


# Per-dataset fixture-install override, keyed by compound dataset name. Datasets absent
# from this mapping use `_install_copytree`. Keep additions here small and commented --
# each installer's own docstring says why that one dataset needs it.
_FIXTURE_INSTALLERS = {
    "cptac:expression": _install_cptac_expression,
    "cptac:phospho": _install_cptac_phospho,
    "ensembl-gene:structure": _install_ensembl_gene_structure,
    "ucsc-cellbrowser:adult-ctx": _install_ucsc_cellbrowser_ctx,
    "ucsc-cellbrowser:dev-ctx": _install_ucsc_cellbrowser_ctx,
    "gtex-eqtl:eqtls": _install_gtex_eqtl,
    "pqtl:metrics": _install_pqtl,
}


def _assert_artifact_nonempty(loaded) -> None:
    """Assert *loaded* (the artifact `reprocess` just wrote) holds at least one row.

    Each artifact type reports its size differently -- mirrors the per-type
    assertions the conformance tests above already use.
    """
    if isinstance(loaded, VariantMatrix):
        assert loaded.n_variants > 0
        assert loaded.n_samples > 0
    elif isinstance(loaded, ExpressionMatrix):
        assert loaded.n_obs > 0
        assert loaded.n_vars > 0
    else:
        assert loaded.count() > 0


_SKILLS_DIR = Path(__file__).resolve().parents[1] / "skills"


def _no_parse_dataset_names() -> list[str]:
    """Every dataset declared with no parse stage (downloaded or bring-your-own alike).

    Reads each `plugin.yaml` directly with `yaml.safe_load` -- the same approach
    `hvantk/tests/test_plugin_acquisition.py::_declared_acquisitions` uses -- instead
    of going through `plugin_loader.get_registry()`. *Merely constructing* the
    registry eagerly imports every plugin's builder module, Hail included
    (`PluginRegistry._build_provider`'s Pass 2 "attempted eagerly here" -- there is no
    lazier public entry point; `list_manifests()` reads `self._manifests`, which this
    same eager pass already populated as a side effect of building the registry in the
    first place). Since this parametrize list is computed at collection time, calling
    `get_registry()` here would cost every pytest invocation of the *default* suite
    that eager import, not only one gated by `-m hail`, just to list test ids. Mirrors
    the loader's own `has_parse_fn = "parse" in (ds_manifest.get("lifecycle") or {})`
    (`PluginRegistry._build_dataset_manifest`). The test body itself still goes through
    the real registry -- Hail import included -- but only for the one parametrized
    case that is actually running.
    """
    names = []
    for manifest_path in sorted(_SKILLS_DIR.glob("*/plugin.yaml")):
        manifest = yaml.safe_load(manifest_path.read_text())
        provider = manifest.get("name") or manifest_path.parent.name
        for ds in manifest.get("datasets", []):
            lifecycle = ds.get("lifecycle") or {}
            if "parse" not in lifecycle:
                names.append(f"{provider}:{ds['name']}")
    return names


@pytest.mark.hail
@pytest.mark.parametrize("dataset_name", _no_parse_dataset_names())
def test_reprocess_builds_every_no_parse_dataset_from_its_raw_directory(
    tmp_path, dataset_name
):
    """`hvantk reprocess <no-parse dataset> --raw-dir <dir> --skip-download` must build.

    Drives the real CLI (`reprocess_cmd` via `CliRunner`), not `run_builder_for_spec`
    directly like the tests above -- those pass the plugin's fixture FILE straight
    through (or, for the few builders that already glob a directory themselves, the
    whole fixture directory unmodified) and so never exercised `reprocess`'s actual
    raw-*directory* contract for a dataset with no `lifecycle.parse`, which is exactly
    where the bug lived.

    Every dataset parametrized here must build successfully offline, from its own
    committed fixture, with `--skip-download --no-check-drift` (no network: the drift
    probe is skipped, and the fixture is already on disk) -- `--skip-download` is legal
    (if redundant) for a bring-your-own dataset too. A dataset whose builder expects the
    raw directory laid out differently than its fixture is committed, or needs a
    non-default plugin arg to build offline at all, gets a small, dedicated installer in
    `_FIXTURE_INSTALLERS` instead of a silent skip.
    """
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset(dataset_name)

    # `--no-check-drift` only skips reprocess's POST-build comparison
    # (reprocess_cmd step 4); run_builder_for_spec still calls spec.drift_probe()
    # during the build itself to stamp provenance's source_fingerprint, and only
    # a DriftProbeError from that call is handled gracefully (PROBE_UNAVAILABLE_TOKEN)
    # -- any other exception a real network hiccup raises would propagate and fail
    # this test. `reg.get_dataset` caches DatasetSpec instances, and `reprocess_cmd`
    # below fetches that same cached instance (no reset in between), so patching it
    # here -- the same isolation every other test in this file applies -- keeps the
    # whole run offline and deterministic regardless of network conditions.
    object.__setattr__(
        spec, "drift_probe", lambda: {"fingerprint": f"sha256:test-{dataset_name}"}
    )

    fixture_dir = Path(spec.test_paths.fixture)
    assert fixture_dir.is_dir(), f"expected a fixture directory: {fixture_dir}"

    raw_dir = tmp_path / "raw"
    install = _FIXTURE_INSTALLERS.get(dataset_name, _install_copytree)
    plugin_args = install(fixture_dir, raw_dir)

    ext = _expected_extensions(spec)[0]
    output = tmp_path / f"out{ext}"

    args = [
        dataset_name,
        "--raw-dir",
        str(raw_dir),
        "--output",
        str(output),
        "--skip-download",
        "--no-check-drift",
    ]
    for key, value in plugin_args.items():
        args += ["--plugin-arg", f"{key}={value}"]

    runner = CliRunner()
    result = runner.invoke(reprocess_cmd, args)

    assert result.exit_code == 0, result.output
    assert output.exists(), f"{dataset_name}: expected output at {output}"
    _assert_artifact_nonempty(core_io.load(output))
