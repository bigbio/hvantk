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
from hvantk.core.models import AnnotationTable
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
