"""Per-plugin conformance tests for the Phase B builder contract.

Each test:
  1. Looks up the dataset's DatasetSpec from the live registry.
  2. Prepares minimal parsed input (using the plugin's own test fixtures).
  3. Calls run_builder_for_spec() with the spec.
  4. Verifies the saved artifact loads back via core/io with the expected type
     and provenance fields.

Each plugin gets its own test as Phase B migrates it.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

from hvantk.core import io as core_io
from hvantk.core.models import AnnotationTable, ExpressionMatrix
from hvantk.core.plugin import loader as plugin_loader
from hvantk.core.plugin.run_builder import run_builder_for_spec


# ---------- peptideatlas:phospho ----------


@pytest.fixture
def peptideatlas_phospho_parsed_tsv(tmp_path):
    """Write a minimal intermediate phospho TSV with the schema produced by parse_raw_dir."""
    cols = [
        "gene_symbol", "protein_accession", "site_position", "residue",
        "n_observations", "peptide_sequences", "source_db", "evidence_type",
    ]
    rows = [
        {"gene_symbol": "TP53", "protein_accession": "P04637",
         "site_position": "315", "residue": "S",
         "n_observations": "2", "peptide_sequences": "ABCDE;FGHIJ",
         "source_db": "peptideatlas", "evidence_type": "phospho"},
        {"gene_symbol": "BRCA1", "protein_accession": "P38398",
         "site_position": "988", "residue": "S",
         "n_observations": "1", "peptide_sequences": "KLMNO",
         "source_db": "peptideatlas", "evidence_type": "phospho"},
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
    assert spec.plugin_version  # populated by Task 1

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
    fixture = Path("hvantk/skills/clinvar/tests/testdata/raw/clinvar/clinvar_20220403_chr20.vcf.bgz")
    assert fixture.exists(), f"missing fixture: {fixture}"

    # Avoid network calls in the drift_probe (clinvar probe hits NCBI)
    object.__setattr__(spec, "drift_probe", lambda: {"fingerprint": "sha256:test-clinvar"})

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

    object.__setattr__(spec, "drift_probe", lambda: {"fingerprint": "sha256:test-gencc"})

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

    fixture = Path("hvantk/skills/clingen/tests/testdata/raw/clingen/clingen_test_sample.csv")
    assert fixture.exists()

    object.__setattr__(spec, "drift_probe", lambda: {"fingerprint": "sha256:test-clingen"})

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

    fixture = Path("hvantk/skills/gwas_catalog/tests/testdata/raw/gwas-catalog/gwas-catalog-sample.tsv")
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

    fixture = Path("hvantk/skills/gtex_eqtl/tests/testdata/raw/gtex-eqtl/Liver.v11.eQTLs.signif_pairs.parquet")
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

    fixture = Path("hvantk/skills/insider/tests/testdata/raw/insider/insider_sample.bed")
    assert fixture.exists()

    object.__setattr__(spec, "drift_probe", lambda: {"fingerprint": "sha256:test-insider"})

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

    object.__setattr__(spec, "drift_probe", lambda: {"fingerprint": "sha256:test-msigdb"})

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
        "chrom", "codon_start", "codon_end", "strand", "uniprot_id", "gene_symbol",
        "residue_pos", "amino_acid", "ptm_type", "ptm_category", "source_db",
        "evidence_type", "n_observations", "tissue_type",
    ]
    rows = [
        # Use chr17 (TP53) and chr13 (BRCA2) for plausibility
        ["17", "7676272", "7676274", "-", "P04637", "TP53", "315", "S",
         "phosphoserine", "phosphorylation", "uniprot", "experimental", "10", ""],
        ["13", "32316461", "32316463", "+", "P51587", "BRCA2", "988", "S",
         "phosphoserine", "phosphorylation", "uniprot", "experimental", "5", ""],
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

    object.__setattr__(spec, "drift_probe", lambda: {"fingerprint": "sha256:test-expr-atlas"})

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

    object.__setattr__(spec, "drift_probe", lambda: {"fingerprint": "sha256:test-cptac-expr"})

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
    expr_lines = [
        "SiteID\tS1\tS2",
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

    object.__setattr__(spec, "drift_probe", lambda: {"fingerprint": "sha256:test-cptac-phospho"})

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
        ("default",   "ucsc-cellbrowser-default-v1"),
        ("adult-ctx", "ucsc-cellbrowser-adult-ctx-v1"),
        ("dev-ctx",   "ucsc-cellbrowser-dev-ctx-v1"),
    ],
)
def test_ucsc_cellbrowser_round_trip(tmp_path, ucsc_cellbrowser_inputs, dataset_name, expected_schema):
    from hvantk.core.models import ExpressionMatrix

    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset(f"ucsc-cellbrowser:{dataset_name}")

    assert spec.artifact_type is ExpressionMatrix
    assert spec.schema_id == expected_schema

    object.__setattr__(spec, "drift_probe", lambda: {"fingerprint": f"sha256:test-ucsc-{dataset_name}"})

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


# ---------- Phase K plugins with fixtures ----------


@pytest.mark.hail
def test_gevir_metrics_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("gevir:metrics")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "gevir-metrics-v1"

    fixture = Path("hvantk/tests/testdata/raw/gevir/gevir_metrics_pmid31873297.tsv.bgz")
    assert fixture.exists()

    object.__setattr__(spec, "drift_probe", lambda: {"source_version": "test"})

    out = tmp_path / "metrics.ht"
    prov = run_builder_for_spec(
        spec, parsed_input=fixture, output_path=out, plugin_version=spec.plugin_version,
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

    fixture = Path("hvantk/tests/testdata/raw/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.chr20.tsv.bgz")
    assert fixture.exists()

    object.__setattr__(spec, "drift_probe", lambda: {"source_version": "test"})

    out = tmp_path / "metrics.ht"
    prov = run_builder_for_spec(
        spec, parsed_input=fixture, output_path=out, plugin_version=spec.plugin_version,
    )
    assert prov.plugin == "gnomad-metrics"
    assert prov.schema_id == "gnomad-metrics-v1"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.count() > 0


@pytest.mark.hail
def test_ensembl_gene_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("ensembl-gene:genes")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "ensembl-gene-v1"

    fixture = Path("hvantk/tests/testdata/raw/ensembl/ensembl_gene_biomart.tsv.bgz")
    assert fixture.exists()

    object.__setattr__(spec, "drift_probe", lambda: {"source_version": "test"})

    out = tmp_path / "genes.ht"
    prov = run_builder_for_spec(
        spec, parsed_input=fixture, output_path=out, plugin_version=spec.plugin_version,
    )
    assert prov.plugin == "ensembl-gene"

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
        spec, parsed_input=fixture, output_path=out, plugin_version=spec.plugin_version,
    )
    assert prov.plugin == "dbnsfp"

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.count() > 0
