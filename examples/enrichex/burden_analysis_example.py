#!/usr/bin/env python3
"""
Burden Analysis Example

This example demonstrates how to use EnrichEx burden testing to test if
cases have excess rare variants in gene set genes compared to controls.

Use Case: Test if Alzheimer's disease cases have increased burden of rare
damaging variants in microglia genes.

Note: This example creates synthetic test data for demonstration.
In real analysis, you would use actual cohort MatrixTables and phenotype data.
"""

import json
import logging
from pathlib import Path
import hail as hl
from hvantk.algorithms.enrichex import run_burden_analysis
from hvantk.core.utils.hail_context import init_hail

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def create_synthetic_cohort(n_samples=1000, n_variants=500):
    """
    Create a synthetic cohort MatrixTable for demonstration.

    In real analysis, you would load your actual data:
        mt = hl.read_matrix_table("real_cohort.mt")

    Parameters
    ----------
    n_samples : int
        Number of samples (500 cases, 500 controls)
    n_variants : int
        Number of variants to generate
    """
    logger.info(
        f"Creating synthetic cohort: {n_samples} samples, {n_variants} variants"
    )

    # Create base MatrixTable
    mt = hl.utils.range_matrix_table(n_variants, n_samples)

    # Add genomic coordinates
    mt = mt.annotate_rows(
        locus=hl.locus("chr17", 41000000 + mt.row_idx, reference_genome="GRCh38"),
        alleles=["A", "T"],
    )
    mt = mt.key_rows_by(mt.locus, mt.alleles)

    # Add sample IDs
    mt = mt.key_cols_by(s=hl.str(mt.col_idx))

    # Add genotypes (random for demonstration)
    # In real data, genotypes come from variant calling
    mt = mt.annotate_entries(
        GT=hl.call(
            hl.if_else(hl.rand_unif(0, 1) < 0.01, 1, 0),  # 1% het rate
            hl.if_else(hl.rand_unif(0, 1) < 0.001, 1, 0),  # 0.1% hom-alt rate
        ),
        GQ=hl.int32(hl.rand_unif(20, 99)),
        DP=hl.int32(hl.rand_unif(10, 100)),
    )

    # Add variant annotations
    # In real data, these come from VEP, dbNSFP, gnomAD, etc.
    gene_pool = [
        "TREM2",
        "CD33",
        "MS4A6A",
        "TYROBP",
        "CSF1R",  # Microglia
        "SLC17A7",
        "GRIN2A",
        "GRIN2B",
        "CAMK2A",  # Excitatory neurons
        "GAD1",
        "GAD2",
        "PVALB",
        "SST",  # Inhibitory neurons
        "GFAP",
        "AQP4",
        "SLC1A2",  # Astrocytes
        "APOE",
        "CLU",
        "BIN1",
        "PICALM",  # AD risk genes
    ]

    mt = mt.annotate_rows(
        # Gene annotation (typically from VEP)
        SYMBOL=hl.literal(gene_pool)[hl.int32(hl.rand_unif(0, len(gene_pool)))],
        # Allele frequency (typically from gnomAD)
        gnomad_af=hl.rand_unif(0.0, 0.01),
        # CADD score (typically from dbNSFP)
        cadd_phred=hl.rand_unif(15.0, 35.0),
        # VEP consequence
        most_severe_consequence=hl.literal(
            [
                "missense_variant",
                "synonymous_variant",
                "frameshift_variant",
                "stop_gained",
            ]
        )[hl.int32(hl.rand_unif(0, 4))],
        # Filter status
        filters=hl.empty_set(hl.tstr),
    )

    logger.info("Synthetic cohort created")
    return mt


def create_phenotype_table(n_samples=1000, n_cases=500):
    """
    Create phenotype table with case/control status and covariates.

    In real analysis, you would load your actual phenotypes:
        phenotypes_ht = hl.read_table("phenotypes.ht")
        phenotypes_ht = hl.import_table("phenotypes.tsv", key="sample_id")
    """
    logger.info(
        f"Creating phenotype table: {n_cases} cases, {n_samples - n_cases} controls"
    )

    # Create table with sample IDs
    ht = hl.utils.range_table(n_samples)
    ht = ht.annotate(s=hl.str(ht.idx))
    ht = ht.key_by(ht.s)

    # Add phenotype and covariates
    ht = ht.annotate(
        # Binary phenotype (case/control)
        is_case=ht.idx < n_cases,
        # Continuous phenotype (e.g., cognitive score)
        cognitive_score=hl.if_else(
            ht.idx < n_cases,
            hl.rand_norm(80.0, 10.0),  # Cases: mean=80, sd=10
            hl.rand_norm(100.0, 10.0),  # Controls: mean=100, sd=10
        ),
        # Covariates
        age=hl.rand_norm(70.0, 10.0),
        sex=hl.if_else(hl.rand_bool(0.5), "M", "F"),
        # Principal components (from ancestry PCA)
        PC1=hl.rand_norm(0, 1),
        PC2=hl.rand_norm(0, 1),
        PC3=hl.rand_norm(0, 1),
        PC4=hl.rand_norm(0, 1),
        PC5=hl.rand_norm(0, 1),
    )

    # Convert sex to numeric
    ht = ht.annotate(sex=hl.if_else(ht.sex == "M", 1, 0))

    logger.info("Phenotype table created")
    return ht


def create_gene_sets():
    """
    Define gene sets for burden testing.

    In real analysis, you would load from JSON:
        with open("gene_sets.json") as f:
            gene_sets = json.load(f)["gene_sets"]
    """
    return {
        "Microglia": [
            "TREM2",
            "CD33",
            "MS4A6A",
            "MS4A4A",
            "TYROBP",
            "CSF1R",
            "C1QA",
            "C1QB",
            "C1QC",
            "CTSS",
        ],
        "Excitatory_Neurons": [
            "SLC17A7",
            "CAMK2A",
            "GRIN2A",
            "GRIN2B",
            "NRGN",
            "SATB2",
            "TBR1",
            "CUX2",
        ],
        "Inhibitory_Neurons": [
            "GAD1",
            "GAD2",
            "SLC32A1",
            "PVALB",
            "SST",
            "VIP",
            "LAMP5",
        ],
        "Astrocytes": ["GFAP", "AQP4", "SLC1A2", "SLC1A3", "ALDOC", "GJA1"],
        "AD_Risk_Genes": ["APOE", "CLU", "CR1", "PICALM", "BIN1", "ABCA7", "SORL1"],
    }


def main():
    """Run burden analysis."""

    # Initialize Hail
    logger.info("Initializing Hail...")
    init_hail(quiet=True, min_block_size=128)

    # Create output directory
    output_dir = Path("enrichex_burden_results")
    output_dir.mkdir(exist_ok=True)
    logger.info(f"Output directory: {output_dir}")

    # Create synthetic data
    logger.info("\n" + "=" * 60)
    logger.info("CREATING SYNTHETIC DATA")
    logger.info("=" * 60)
    logger.info("Note: In real analysis, you would load actual data:")
    logger.info("  mt = hl.read_matrix_table('cohort.mt')")
    logger.info("  phenotypes_ht = hl.read_table('phenotypes.ht')")

    mt = create_synthetic_cohort(n_samples=1000, n_variants=500)
    phenotypes_ht = create_phenotype_table(n_samples=1000, n_cases=500)
    gene_sets = create_gene_sets()

    # Save gene sets
    gene_sets_full = {
        "background_genes": list(
            set([g for genes in gene_sets.values() for g in genes])
        ),
        "gene_sets": {
            name: {"name": name, "genes": genes} for name, genes in gene_sets.items()
        },
    }
    gene_sets_path = output_dir / "gene_sets.json"
    with open(gene_sets_path, "w") as f:
        json.dump(gene_sets_full, f, indent=2)
    logger.info(f"\nSaved gene sets to: {gene_sets_path}")

    # Run burden analysis - Binary phenotype
    logger.info("\n" + "=" * 60)
    logger.info("RUNNING BURDEN ANALYSIS: BINARY PHENOTYPE (CASE/CONTROL)")
    logger.info("=" * 60)

    results_binary_ht = run_burden_analysis(
        cohort_mt=mt,
        phenotype_ht=phenotypes_ht,
        gene_sets=gene_sets,
        phenotype_field="is_case",
        phenotype_type="binary",
        covariate_fields=["PC1", "PC2", "PC3", "PC4", "PC5", "age", "sex"],
        max_af=0.001,  # Rare variants (AF < 0.1%)
        min_score=25.0,  # High CADD scores
        genotype_aggregation="hets",
    )

    # Export results
    binary_results_path = output_dir / "burden_binary_results.tsv"
    results_binary_ht.export(str(binary_results_path))
    logger.info(f"\nSaved binary phenotype results to: {binary_results_path}")

    # Display results
    results_binary_df = results_binary_ht.to_pandas()
    results_binary_df = results_binary_df.sort_values("p_value")

    logger.info("\n" + "=" * 60)
    logger.info("BINARY PHENOTYPE RESULTS")
    logger.info("=" * 60)
    logger.info(f"{'Gene Set':<25} {'OR':<8} {'95% CI':<20} {'P-value':<12} {'Sig':<5}")
    logger.info("-" * 80)

    for _, row in results_binary_df.iterrows():
        ci_str = f"[{row['ci_lower']:.2f}, {row['ci_upper']:.2f}]"
        sig_marker = "***" if row["significant"] else ""
        logger.info(
            f"{row['gene_set_name']:<25} "
            f"{row['odds_ratio']:>7.2f} "
            f"{ci_str:<20} "
            f"{row['p_adjusted']:>11.2e} "
            f"{sig_marker:<5}"
        )

    # Run burden analysis - Continuous phenotype
    logger.info("\n" + "=" * 60)
    logger.info("RUNNING BURDEN ANALYSIS: CONTINUOUS PHENOTYPE")
    logger.info("=" * 60)

    results_continuous_ht = run_burden_analysis(
        cohort_mt=mt,
        phenotype_ht=phenotypes_ht,
        gene_sets=gene_sets,
        phenotype_field="cognitive_score",
        phenotype_type="continuous",
        covariate_fields=["PC1", "PC2", "PC3", "PC4", "PC5", "age", "sex"],
        max_af=0.001,
        min_score=25.0,
        genotype_aggregation="hets",
    )

    # Export results
    continuous_results_path = output_dir / "burden_continuous_results.tsv"
    results_continuous_ht.export(str(continuous_results_path))
    logger.info(f"\nSaved continuous phenotype results to: {continuous_results_path}")

    # Display results
    results_continuous_df = results_continuous_ht.to_pandas()
    results_continuous_df = results_continuous_df.sort_values("p_value")

    logger.info("\n" + "=" * 60)
    logger.info("CONTINUOUS PHENOTYPE RESULTS")
    logger.info("=" * 60)
    logger.info(f"{'Gene Set':<25} {'Beta':<10} {'SE':<10} {'P-value':<12} {'Sig':<5}")
    logger.info("-" * 70)

    for _, row in results_continuous_df.iterrows():
        sig_marker = "***" if row["significant"] else ""
        logger.info(
            f"{row['gene_set_name']:<25} "
            f"{row['beta']:>9.2f} "
            f"{row['standard_error']:>9.2f} "
            f"{row['p_adjusted']:>11.2e} "
            f"{sig_marker:<5}"
        )

    # Interpretation
    logger.info("\n" + "=" * 60)
    logger.info("INTERPRETATION GUIDE")
    logger.info("=" * 60)

    logger.info(
        """
Binary Phenotype (Case/Control):
- Odds Ratio (OR): Effect size
  - OR > 1: Cases have more burden than controls (risk factor)
  - OR < 1: Cases have less burden than controls (protective)
  - OR = 1: No difference

- Example interpretation:
  "Microglia: OR=1.57, 95% CI [1.24, 1.98], p=0.0002"
  → Cases have 57% more rare variants in microglia genes
  → Each additional rare variant increases AD risk by 57%
  → Highly significant (p=0.0002)

Continuous Phenotype (Cognitive Score):
- Beta: Effect size per variant
  - Positive beta: Higher burden → higher phenotype value
  - Negative beta: Higher burden → lower phenotype value

- Example interpretation:
  "Microglia: beta=-2.35, SE=0.65, p=0.0003"
  → Each additional rare variant in microglia genes
     decreases cognitive score by 2.35 points
  → Effect is precise (small SE)
  → Highly significant (p=0.0003)

Statistical Significance:
- p < 0.001: Very strong evidence
- p < 0.01: Strong evidence
- p < 0.05: Moderate evidence
- p > 0.05: No significant evidence

Next Steps:
1. Replicate findings in independent cohort
2. Identify specific genes driving burden
3. Test different variant filters (AF, CADD thresholds)
4. Perform stratified analyses (by sex, age, ancestry)
5. Combine with expression data for functional validation
    """
    )

    logger.info("\n" + "=" * 60)
    logger.info("ALTERNATIVE GENOTYPE AGGREGATION METHODS")
    logger.info("=" * 60)

    logger.info(
        """
You can test different genetic models by changing --genotype-aggregation:

1. Heterozygous (default): --genotype-aggregation hets
   - Counts genes with ≥1 heterozygous variant per sample
   - Most common model for complex diseases

2. Homozygous: --genotype-aggregation homs
   - Counts genes with ≥1 homozygous variant per sample
   - Use for recessive disorders

3. Compound Heterozygous: --genotype-aggregation chets
   - Counts genes with ≥2 heterozygous variants per sample
   - Tests for recessive model without requiring homozygotes

4. Homozygous + Compound Heterozygous: --genotype-aggregation homs_chets
   - Counts genes with homs OR ≥2 hets per sample
   - Combined recessive model

Example:
    hvantk enrichex burden ... --genotype-aggregation chets
    """
    )

    logger.info("\n" + "=" * 60)
    logger.info("VARIANT FILTERING STRATEGIES")
    logger.info("=" * 60)

    logger.info(
        """
Conservative (high confidence):
    --max-af 0.0001 --min-cadd 30 --consequences frameshift_variant,stop_gained
    Use for: Rare disease studies, high-penetrance variants

Moderate (balanced):
    --max-af 0.001 --min-cadd 25 --consequences missense_variant,frameshift_variant
    Use for: Complex diseases, general burden testing (DEFAULT)

Permissive (exploratory):
    --max-af 0.01 --min-cadd 20 --consequences missense_variant,synonymous_variant
    Use for: Exploratory analyses, large cohorts
    """
    )

    logger.info("\n" + "=" * 60)
    logger.info("Analysis complete!")
    logger.info("=" * 60)
    logger.info(f"\nOutput files:")
    logger.info(f"  - Gene sets: {gene_sets_path}")
    logger.info(f"  - Binary results: {binary_results_path}")
    logger.info(f"  - Continuous results: {continuous_results_path}")

    logger.info(
        """
\nTo run with your own data:

# Python API:
from hvantk.algorithms.enrichex import run_burden_analysis
import hail as hl

mt = hl.read_matrix_table("your_cohort.mt")
phenotypes_ht = hl.read_table("your_phenotypes.ht")
# ... configure and run ...

# Command-line:
hvantk enrichex burden \\
  -m your_cohort.mt \\
  -p your_phenotypes.ht \\
  -s gene_sets.json \\
  --max-af 0.001 \\
  --min-cadd 25 \\
  --covariates PC1,PC2,PC3,age,sex \\
  -o results.tsv
    """
    )


if __name__ == "__main__":
    main()
