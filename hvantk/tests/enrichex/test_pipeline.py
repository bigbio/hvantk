"""
Tests for BurdenPipeline orchestrator (Phase 3.2).
"""

import json

import pandas as pd
import pytest

from hvantk.enrichex.pipeline import BurdenConfig, BurdenPipeline, BurdenRunResult
from hvantk.utils.table_utils import leaf_name


# ---------------------------------------------------------------------------
# BurdenConfig validation tests (no Hail needed)
# ---------------------------------------------------------------------------


class TestBurdenConfig:
    """Tests for BurdenConfig validation."""

    def _write_dummy_file(self, path):
        """Create a minimal file so path-existence checks pass."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("")

    def _write_dummy_gene_sets(self, path):
        """Write a minimal GeneSetCollection JSON."""
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "gene_sets": {
                "set1": {"name": "set1", "genes": ["GENE1", "GENE2"]},
            },
            "background_genes": ["GENE1", "GENE2", "GENE3"],
        }
        path.write_text(json.dumps(data))

    def test_valid_config(self, tmp_path):
        mt_path = tmp_path / "cohort.mt"
        mt_path.mkdir()
        pheno_path = tmp_path / "pheno.ht"
        pheno_path.mkdir()
        gs_path = tmp_path / "genes.json"
        self._write_dummy_gene_sets(gs_path)

        config = BurdenConfig(
            cohort_mt_path=str(mt_path),
            phenotype_ht_path=str(pheno_path),
            gene_set_collections={"heart": str(gs_path)},
        )
        errors = config.validate()
        assert errors == []

    def test_missing_cohort_path(self, tmp_path):
        config = BurdenConfig(
            cohort_mt_path="",
            phenotype_ht_path=str(tmp_path),
            gene_set_collections={"x": str(tmp_path / "gs.json")},
        )
        errors = config.validate()
        assert any("cohort_mt_path" in e for e in errors)

    def test_missing_gene_set_collections(self, tmp_path):
        mt_path = tmp_path / "cohort.mt"
        mt_path.mkdir()
        config = BurdenConfig(
            cohort_mt_path=str(mt_path),
            phenotype_ht_path=str(tmp_path),
            gene_set_collections={},
        )
        errors = config.validate()
        assert any("gene_set_collection" in e for e in errors)

    def test_invalid_phenotype_type(self, tmp_path):
        mt_path = tmp_path / "cohort.mt"
        mt_path.mkdir()
        gs_path = tmp_path / "gs.json"
        self._write_dummy_gene_sets(gs_path)
        config = BurdenConfig(
            cohort_mt_path=str(mt_path),
            phenotype_ht_path=str(tmp_path),
            gene_set_collections={"x": str(gs_path)},
            phenotype_type="invalid",
        )
        errors = config.validate()
        assert any("phenotype_type" in e for e in errors)

    def test_invalid_correction_method(self, tmp_path):
        mt_path = tmp_path / "cohort.mt"
        mt_path.mkdir()
        gs_path = tmp_path / "gs.json"
        self._write_dummy_gene_sets(gs_path)
        config = BurdenConfig(
            cohort_mt_path=str(mt_path),
            phenotype_ht_path=str(tmp_path),
            gene_set_collections={"x": str(gs_path)},
            correction_method="foobar",
        )
        errors = config.validate()
        assert any("correction_method" in e for e in errors)

    def test_invalid_genotype_aggregation(self, tmp_path):
        mt_path = tmp_path / "cohort.mt"
        mt_path.mkdir()
        gs_path = tmp_path / "gs.json"
        self._write_dummy_gene_sets(gs_path)
        config = BurdenConfig(
            cohort_mt_path=str(mt_path),
            phenotype_ht_path=str(tmp_path),
            gene_set_collections={"x": str(gs_path)},
            genotype_aggregation="invalid_agg",
        )
        errors = config.validate()
        assert any("genotype_aggregation" in e for e in errors)

    def test_nonexistent_gene_set_path(self, tmp_path):
        mt_path = tmp_path / "cohort.mt"
        mt_path.mkdir()
        config = BurdenConfig(
            cohort_mt_path=str(mt_path),
            phenotype_ht_path=str(tmp_path),
            gene_set_collections={"x": "/nonexistent/path.json"},
        )
        errors = config.validate()
        assert any("not found" in e for e in errors)


class TestBurdenPipelineHelpers:
    """Test helper utilities used by BurdenPipeline."""

    def test_leaf_name_simple(self):
        assert leaf_name("is_case") == "is_case"

    def test_leaf_name_nested(self):
        assert leaf_name("phe.is_case") == "is_case"

    def test_leaf_name_deeply_nested(self):
        assert leaf_name("a.b.c") == "c"


class TestBurdenRunResult:
    """Test BurdenRunResult dataclass."""

    def test_creation(self):
        r = BurdenRunResult(
            variant_class="lof",
            collection_name="heart",
            results_df=pd.DataFrame({"p_value": [0.01]}),
            n_gene_sets_tested=1,
            n_significant=1,
        )
        assert r.variant_class == "lof"
        assert r.n_significant == 1


# ---------------------------------------------------------------------------
# BurdenPipeline tests
# ---------------------------------------------------------------------------


class TestBurdenPipelineShowPlan:
    """Test show_plan() — no Hail needed."""

    def test_show_plan_no_hail(self, tmp_path):
        mt_path = tmp_path / "cohort.mt"
        mt_path.mkdir()
        pheno_path = tmp_path / "pheno.ht"
        pheno_path.mkdir()
        gs_path = tmp_path / "gs.json"
        gs_path.write_text(
            json.dumps(
                {
                    "gene_sets": {"s1": {"name": "s1", "genes": ["G1"]}},
                    "background_genes": ["G1", "G2"],
                }
            )
        )

        config = BurdenConfig(
            cohort_mt_path=str(mt_path),
            phenotype_ht_path=str(pheno_path),
            gene_set_collections={"heart": str(gs_path)},
            output_dir=str(tmp_path / "out"),
        )
        pipeline = BurdenPipeline(config)
        plan = pipeline.show_plan()
        assert "BURDEN PIPELINE EXECUTION PLAN" in plan
        assert "heart" in plan

    def test_show_plan_with_variant_classes(self, tmp_path):
        from hvantk.enrichex.burden import VariantFilter

        mt_path = tmp_path / "cohort.mt"
        mt_path.mkdir()
        pheno_path = tmp_path / "pheno.ht"
        pheno_path.mkdir()
        gs_path = tmp_path / "gs.json"
        gs_path.write_text(
            json.dumps(
                {
                    "gene_sets": {"s1": {"name": "s1", "genes": ["G1"]}},
                    "background_genes": ["G1"],
                }
            )
        )

        config = BurdenConfig(
            cohort_mt_path=str(mt_path),
            phenotype_ht_path=str(pheno_path),
            gene_set_collections={"heart": str(gs_path)},
            variant_classes={
                "lof": VariantFilter(consequences=["stop_gained"]),
                "syn": VariantFilter(consequences=["synonymous_variant"]),
            },
            output_dir=str(tmp_path / "out"),
        )
        pipeline = BurdenPipeline(config)
        plan = pipeline.show_plan()
        assert "lof" in plan
        assert "syn" in plan
        assert "2 classes x 1 collections" in plan

    def test_invalid_config_raises(self, tmp_path):
        config = BurdenConfig(
            cohort_mt_path="",
            phenotype_ht_path="",
        )
        with pytest.raises(ValueError, match="validation failed"):
            BurdenPipeline(config)

    def test_valid_config_no_phenotype_ht(self, tmp_path):
        """Config is valid when phenotype_ht_path is empty (MT col fields)."""
        mt_path = tmp_path / "cohort.mt"
        mt_path.mkdir()
        gs_path = tmp_path / "gs.json"
        gs_path.write_text(
            json.dumps(
                {
                    "gene_sets": {"s1": {"name": "s1", "genes": ["G1"]}},
                    "background_genes": ["G1"],
                }
            )
        )

        config = BurdenConfig(
            cohort_mt_path=str(mt_path),
            phenotype_ht_path="",  # extract from MT cols
            phenotype_field="phe.is_case",
            gene_set_collections={"heart": str(gs_path)},
        )
        errors = config.validate()
        assert errors == []

    def test_show_plan_mt_col_fields(self, tmp_path):
        """show_plan indicates phenotype from MT column fields."""
        mt_path = tmp_path / "cohort.mt"
        mt_path.mkdir()
        gs_path = tmp_path / "gs.json"
        gs_path.write_text(
            json.dumps(
                {
                    "gene_sets": {"s1": {"name": "s1", "genes": ["G1"]}},
                    "background_genes": ["G1"],
                }
            )
        )

        config = BurdenConfig(
            cohort_mt_path=str(mt_path),
            phenotype_ht_path="",
            phenotype_field="phe.is_case",
            gene_set_collections={"heart": str(gs_path)},
            output_dir=str(tmp_path / "out"),
        )
        pipeline = BurdenPipeline(config)
        plan = pipeline.show_plan()
        assert "(from MT column fields)" in plan


# ---------------------------------------------------------------------------
# Integration tests (Hail required)
# ---------------------------------------------------------------------------


def _create_test_data(tmp_path):
    """Create and write synthetic MT, phenotype HT, and gene sets to disk."""
    import hail as hl

    # Cohort MT: 30 variants, 50 samples
    mt = hl.utils.range_matrix_table(30, 50)
    mt = mt.key_cols_by(s=hl.str(mt.col_idx))
    mt = mt.annotate_rows(
        SYMBOL=hl.if_else(
            mt.row_idx < 10,
            "GENE1",
            hl.if_else(mt.row_idx < 20, "GENE2", "GENE3"),
        ),
        gnomad_af=hl.rand_unif(0, 0.05),
        cadd_phred=hl.rand_unif(15, 30),
        consequence=hl.if_else(mt.row_idx % 3 == 0, "stop_gained", "missense_variant"),
    )
    mt = mt.annotate_entries(
        GT=hl.call(hl.int(hl.rand_bool(0.3)), hl.int(hl.rand_bool(0.1))),
        GQ=30,
        DP=20,
    )
    mt_path = str(tmp_path / "cohort.mt")
    mt.write(mt_path, overwrite=True)

    # Phenotype HT: binary phenotype + covariates
    ht = hl.utils.range_table(50)
    ht = ht.annotate(s=hl.str(ht.idx))
    ht = ht.key_by("s")
    ht = ht.annotate(
        is_case=ht.idx < 25,
        PC1=hl.rand_norm(),
    )
    pheno_path = str(tmp_path / "pheno.ht")
    ht.write(pheno_path, overwrite=True)

    # Gene set collection
    gs_data = {
        "gene_sets": {
            "set1": {"name": "set1", "genes": ["GENE1", "GENE2"]},
            "set2": {"name": "set2", "genes": ["GENE3"]},
        },
        "background_genes": ["GENE1", "GENE2", "GENE3", "GENE4"],
    }
    gs_path = tmp_path / "gene_sets.json"
    gs_path.write_text(json.dumps(gs_data))

    return mt_path, pheno_path, str(gs_path)


def _create_test_data_with_nested_cols(tmp_path):
    """Create MT with nested column structs (phenotype in col fields)."""
    import hail as hl

    # Cohort MT: 30 variants, 50 samples, phenotype in col fields
    mt = hl.utils.range_matrix_table(30, 50)
    mt = mt.key_cols_by(s=hl.str(mt.col_idx))
    mt = mt.annotate_rows(
        SYMBOL=hl.if_else(
            mt.row_idx < 10,
            "GENE1",
            hl.if_else(mt.row_idx < 20, "GENE2", "GENE3"),
        ),
        gnomad_af=hl.rand_unif(0, 0.05),
        cadd_phred=hl.rand_unif(15, 30),
        consequence=hl.if_else(mt.row_idx % 3 == 0, "stop_gained", "missense_variant"),
    )
    mt = mt.annotate_entries(
        GT=hl.call(hl.int(hl.rand_bool(0.3)), hl.int(hl.rand_bool(0.1))),
        GQ=30,
        DP=20,
    )
    # Nested phenotype struct in column fields (like CHD MT)
    mt = mt.annotate_cols(
        phe=hl.struct(
            is_case=mt.col_idx < 25,
            PC1=hl.rand_norm(),
        ),
    )
    mt_path = str(tmp_path / "cohort_nested.mt")
    mt.write(mt_path, overwrite=True)

    # Gene set collection
    gs_data = {
        "gene_sets": {
            "set1": {"name": "set1", "genes": ["GENE1", "GENE2"]},
            "set2": {"name": "set2", "genes": ["GENE3"]},
        },
        "background_genes": ["GENE1", "GENE2", "GENE3", "GENE4"],
    }
    gs_path = tmp_path / "gene_sets.json"
    gs_path.write_text(json.dumps(gs_data))

    return mt_path, str(gs_path)


@pytest.mark.hail
class TestBurdenPipelineRun:
    """Integration tests for BurdenPipeline.run()."""

    def test_run_single_collection_no_classes(self, hail_session, tmp_path):
        """Run without variant class stratification."""
        mt_path, pheno_path, gs_path = _create_test_data(tmp_path)

        config = BurdenConfig(
            cohort_mt_path=mt_path,
            phenotype_ht_path=pheno_path,
            gene_set_collections={"test": gs_path},
            min_carriers=0,
            output_dir=str(tmp_path / "out"),
            generate_report=False,
        )
        pipeline = BurdenPipeline(config)
        combined_df = pipeline.run()

        assert not combined_df.empty
        assert "gene_set_name" in combined_df.columns
        assert "p_value" in combined_df.columns
        assert "p_adjusted" in combined_df.columns
        assert "variant_class" in combined_df.columns
        assert "collection" in combined_df.columns

        # Check per-run TSV was written
        per_run_files = list((tmp_path / "out" / "per_run").glob("*.tsv"))
        assert len(per_run_files) >= 1

        # Check combined TSV
        combined_path = tmp_path / "out" / "burden_combined.tsv"
        assert combined_path.exists()

    def test_run_with_variant_classes(self, hail_session, tmp_path):
        """Run with variant class stratification."""
        from hvantk.enrichex.burden import VariantFilter

        mt_path, pheno_path, gs_path = _create_test_data(tmp_path)

        config = BurdenConfig(
            cohort_mt_path=mt_path,
            phenotype_ht_path=pheno_path,
            gene_set_collections={"test": gs_path},
            variant_classes={
                "lof": VariantFilter(
                    consequences=["stop_gained"],
                    pass_only=False,
                    min_gq=0,
                    min_dp=0,
                    min_score=None,
                ),
                "missense": VariantFilter(
                    consequences=["missense_variant"],
                    pass_only=False,
                    min_gq=0,
                    min_dp=0,
                    min_score=None,
                ),
            },
            min_carriers=0,
            output_dir=str(tmp_path / "out"),
            generate_report=False,
        )
        pipeline = BurdenPipeline(config)
        combined_df = pipeline.run()

        # Should have results from multiple variant classes
        assert set(combined_df["variant_class"].unique()) <= {"lof", "missense"}

        # Per-run TSVs
        per_run_files = list((tmp_path / "out" / "per_run").glob("burden_*.tsv"))
        assert len(per_run_files) >= 1

    def test_run_multi_collection(self, hail_session, tmp_path):
        """Run with multiple gene set collections."""
        mt_path, pheno_path, gs_path = _create_test_data(tmp_path)

        # Create a second collection
        gs2_data = {
            "gene_sets": {
                "setA": {"name": "setA", "genes": ["GENE1"]},
            },
            "background_genes": ["GENE1", "GENE2", "GENE3"],
        }
        gs2_path = tmp_path / "gene_sets2.json"
        gs2_path.write_text(json.dumps(gs2_data))

        config = BurdenConfig(
            cohort_mt_path=mt_path,
            phenotype_ht_path=pheno_path,
            gene_set_collections={
                "collection1": gs_path,
                "collection2": str(gs2_path),
            },
            min_carriers=0,
            output_dir=str(tmp_path / "out"),
            generate_report=False,
        )
        pipeline = BurdenPipeline(config)
        combined_df = pipeline.run()

        assert not combined_df.empty
        assert set(combined_df["collection"].unique()) == {
            "collection1",
            "collection2",
        }

    def test_run_empty_results_handled(self, hail_session, tmp_path):
        """Pipeline handles gene sets with no matching genes gracefully."""
        mt_path, pheno_path, _ = _create_test_data(tmp_path)

        # Gene sets with non-existent genes
        gs_data = {
            "gene_sets": {
                "empty_set": {
                    "name": "empty_set",
                    "genes": ["NONEXISTENT1", "NONEXISTENT2"],
                },
            },
            "background_genes": ["NONEXISTENT1", "NONEXISTENT2"],
        }
        gs_path = tmp_path / "empty_gs.json"
        gs_path.write_text(json.dumps(gs_data))

        config = BurdenConfig(
            cohort_mt_path=mt_path,
            phenotype_ht_path=pheno_path,
            gene_set_collections={"empty": str(gs_path)},
            min_carriers=0,
            output_dir=str(tmp_path / "out"),
            generate_report=False,
        )
        pipeline = BurdenPipeline(config)
        combined_df = pipeline.run()

        assert combined_df.empty

    def test_run_with_covariates(self, hail_session, tmp_path):
        """Pipeline passes covariates through correctly."""
        mt_path, pheno_path, gs_path = _create_test_data(tmp_path)

        config = BurdenConfig(
            cohort_mt_path=mt_path,
            phenotype_ht_path=pheno_path,
            gene_set_collections={"test": gs_path},
            covariate_fields=["PC1"],
            min_carriers=0,
            output_dir=str(tmp_path / "out"),
            generate_report=False,
        )
        pipeline = BurdenPipeline(config)
        combined_df = pipeline.run()

        assert not combined_df.empty

    def test_run_results_tracked(self, hail_session, tmp_path):
        """Pipeline tracks run results in _run_results."""
        mt_path, pheno_path, gs_path = _create_test_data(tmp_path)

        config = BurdenConfig(
            cohort_mt_path=mt_path,
            phenotype_ht_path=pheno_path,
            gene_set_collections={"test": gs_path},
            min_carriers=0,
            output_dir=str(tmp_path / "out"),
            generate_report=False,
        )
        pipeline = BurdenPipeline(config)
        pipeline.run()

        assert len(pipeline._run_results) == 1
        r = pipeline._run_results[0]
        assert isinstance(r, BurdenRunResult)
        assert r.collection_name == "test"
        assert r.variant_class == "all"
        assert r.n_gene_sets_tested > 0

    def test_run_with_mt_col_fields(self, hail_session, tmp_path):
        """Pipeline extracts phenotype from MT column fields (nested struct)."""
        mt_path, gs_path = _create_test_data_with_nested_cols(tmp_path)

        config = BurdenConfig(
            cohort_mt_path=mt_path,
            phenotype_ht_path="",  # extract from MT cols
            phenotype_field="phe.is_case",
            covariate_fields=["phe.PC1"],
            gene_set_collections={"test": gs_path},
            min_carriers=0,
            output_dir=str(tmp_path / "out"),
            generate_report=False,
        )
        pipeline = BurdenPipeline(config)
        combined_df = pipeline.run()

        assert not combined_df.empty
        assert "gene_set_name" in combined_df.columns
        assert "p_value" in combined_df.columns

        # Check that config fields were flattened to leaf names
        assert pipeline.config.phenotype_field == "is_case"
        assert pipeline.config.covariate_fields == ["PC1"]

    def test_run_with_mt_col_fields_flat(self, hail_session, tmp_path):
        """Pipeline extracts flat (non-nested) phenotype from MT column fields."""
        import hail as hl

        # MT with flat column fields
        mt = hl.utils.range_matrix_table(30, 50)
        mt = mt.key_cols_by(s=hl.str(mt.col_idx))
        mt = mt.annotate_rows(
            SYMBOL=hl.if_else(
                mt.row_idx < 10,
                "GENE1",
                hl.if_else(mt.row_idx < 20, "GENE2", "GENE3"),
            ),
            gnomad_af=hl.rand_unif(0, 0.05),
            cadd_phred=hl.rand_unif(15, 30),
            consequence=hl.if_else(
                mt.row_idx % 3 == 0, "stop_gained", "missense_variant"
            ),
        )
        mt = mt.annotate_entries(
            GT=hl.call(hl.int(hl.rand_bool(0.3)), hl.int(hl.rand_bool(0.1))),
            GQ=30,
            DP=20,
        )
        mt = mt.annotate_cols(
            is_case=mt.col_idx < 25,
        )
        mt_path = str(tmp_path / "cohort_flat.mt")
        mt.write(mt_path, overwrite=True)

        gs_data = {
            "gene_sets": {
                "set1": {"name": "set1", "genes": ["GENE1", "GENE2"]},
            },
            "background_genes": ["GENE1", "GENE2", "GENE3"],
        }
        gs_path = tmp_path / "gene_sets.json"
        gs_path.write_text(json.dumps(gs_data))

        config = BurdenConfig(
            cohort_mt_path=mt_path,
            phenotype_ht_path="",
            phenotype_field="is_case",
            gene_set_collections={"test": str(gs_path)},
            min_carriers=0,
            output_dir=str(tmp_path / "out"),
            generate_report=False,
        )
        pipeline = BurdenPipeline(config)
        combined_df = pipeline.run()

        assert not combined_df.empty
