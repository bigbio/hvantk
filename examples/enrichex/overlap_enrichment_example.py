#!/usr/bin/env python3
"""
Overlap Enrichment Example

This example demonstrates how to use EnrichEx overlap enrichment to test
if a gene list (e.g., from GWAS) is enriched in cell-type specific gene sets.

Use Case: Test if Alzheimer's disease GWAS genes are enriched in brain cell types.
"""

import logging
from pathlib import Path
from hvantk.algorithms.enrichex import (
    GeneSet,
    GeneSetCollection,
    compute_overlap_enrichment_pandas,
)
from hvantk.core.hail_context import init_hail

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def create_example_gene_sets() -> GeneSetCollection:
    """Create example brain cell-type gene sets."""

    # Microglia markers (immune cells in brain)
    microglia = GeneSet(
        name="Microglia",
        genes={
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
            "HLA-DRA",
            "HLA-DRB1",
            "AIF1",
            "CX3CR1",
            "P2RY12",
        },
        source="Lake et al. Nat Biotech 2018",
        metadata={"tissue": "brain", "technology": "snRNA-seq", "species": "human"},
    )

    # Excitatory neuron markers
    excitatory = GeneSet(
        name="Excitatory_Neurons",
        genes={
            "SLC17A7",
            "CAMK2A",
            "GRIN2A",
            "GRIN2B",
            "NRGN",
            "SATB2",
            "TBR1",
            "CUX2",
            "RORB",
            "FEZF2",
            "FOXP2",
            "BCL11B",
            "CRYM",
            "LPL",
            "PRSS12",
        },
        source="Lake et al. Nat Biotech 2018",
        metadata={"cell_type": "neuron", "subtype": "excitatory"},
    )

    # Inhibitory neuron markers
    inhibitory = GeneSet(
        name="Inhibitory_Neurons",
        genes={
            "GAD1",
            "GAD2",
            "SLC32A1",
            "PVALB",
            "SST",
            "VIP",
            "LAMP5",
            "ADARB2",
            "LHX6",
            "NKX2-1",
            "SOX6",
            "CNR1",
            "CCK",
            "TAC1",
            "RELN",
        },
        source="Lake et al. Nat Biotech 2018",
        metadata={"cell_type": "neuron", "subtype": "inhibitory"},
    )

    # Astrocyte markers
    astrocytes = GeneSet(
        name="Astrocytes",
        genes={
            "GFAP",
            "AQP4",
            "SLC1A2",
            "SLC1A3",
            "ALDOC",
            "GJA1",
            "SLC14A1",
            "BMPR1B",
            "MFGE8",
            "AGT",
            "ATP1A2",
            "GLUL",
            "S100B",
            "ALDH1L1",
            "FGFR3",
        },
        source="Lake et al. Nat Biotech 2018",
        metadata={"cell_type": "glia", "subtype": "astrocyte"},
    )

    # Oligodendrocyte markers
    oligodendrocytes = GeneSet(
        name="Oligodendrocytes",
        genes={
            "MBP",
            "MOG",
            "MAG",
            "PLP1",
            "MOBP",
            "OLIG1",
            "OLIG2",
            "SOX10",
            "CLDN11",
            "CNP",
            "ERMN",
            "GJC2",
            "ENPP2",
            "UGT8",
            "FA2H",
        },
        source="Lake et al. Nat Biotech 2018",
        metadata={"cell_type": "glia", "subtype": "oligodendrocyte"},
    )

    # Endothelial cells
    endothelial = GeneSet(
        name="Endothelial_Cells",
        genes={
            "CLDN5",
            "FLT1",
            "PECAM1",
            "VWF",
            "NOSTRIN",
            "ITM2A",
            "SLC2A1",
            "ABCB1",
            "ABCG2",
            "SLC38A5",
            "ADGRL4",
            "KDR",
            "ESAM",
            "CDH5",
            "TIE1",
        },
        source="Lake et al. Nat Biotech 2018",
        metadata={"cell_type": "vascular", "subtype": "endothelial"},
    )

    # Define background: all protein-coding genes expressed in brain
    # For this example, we'll use a union of all gene sets + some additional genes
    all_genes = set()
    for gs in [
        microglia,
        excitatory,
        inhibitory,
        astrocytes,
        oligodendrocytes,
        endothelial,
    ]:
        all_genes.update(gs.genes)

    # Add some additional common genes to background
    all_genes.update(
        [
            "GAPDH",
            "ACTB",
            "TUBB",
            "APOE",
            "APP",
            "MAPT",
            "SNCA",
            "PSEN1",
            "PSEN2",
            "CLU",
            "CR1",
            "PICALM",
            "BIN1",
            "ABCA7",
            "INPP5D",
            "MEF2C",
            "HLA-DRB5",
            "PTK2B",
            "SORL1",
            "SLC24A4",
            "CASS4",
            "FERMT2",
        ]
    )

    # Create collection
    collection = GeneSetCollection(
        gene_sets={
            "Microglia": microglia,
            "Excitatory_Neurons": excitatory,
            "Inhibitory_Neurons": inhibitory,
            "Astrocytes": astrocytes,
            "Oligodendrocytes": oligodendrocytes,
            "Endothelial_Cells": endothelial,
        },
        background_genes=all_genes,
        source_description="Brain cell-type markers from Lake et al. 2018",
        metadata={
            "reference": "Lake et al. Nature Biotechnology 2018",
            "pmid": "29227469",
            "tissue": "brain",
            "technology": "snRNA-seq",
        },
    )

    return collection


def get_alzheimers_gwas_genes():
    """
    Return list of genes associated with Alzheimer's disease from GWAS.

    These are well-established AD risk genes from large-scale GWAS studies.
    """
    return [
        # Strongest associations
        "APOE",  # APOE epsilon4 allele - strongest genetic risk factor
        "TREM2",  # Microglia receptor - rare variants
        # Immune/microglia genes
        "CD33",  # Microglia - inhibits Aβ clearance
        "MS4A6A",  # Microglia - immune regulation
        "MS4A4A",  # Microglia - immune regulation
        "ABCA7",  # Microglia - lipid transport
        "CR1",  # Complement receptor
        "INPP5D",  # Microglia - PI3K pathway
        "HLA-DRB1",  # MHC class II
        "HLA-DRB5",  # MHC class II
        # Endocytosis/trafficking
        "PICALM",  # Clathrin-mediated endocytosis
        "BIN1",  # Endocytosis
        "SORL1",  # APP trafficking
        "RIN3",  # Endosomal trafficking
        # Lipid metabolism
        "CLU",  # Clusterin/apolipoprotein J
        "PTK2B",  # Focal adhesion kinase
        # Other pathways
        "FERMT2",  # Cell adhesion
        "SLC24A4",  # Calcium homeostasis
        "CASS4",  # Cytoskeleton
        "MEF2C",  # Transcription factor
    ]


def main():
    """Run overlap enrichment analysis."""

    # Initialize Hail
    logger.info("Initializing Hail...")
    init_hail(quiet=True, min_block_size=128)

    # Create output directory
    output_dir = Path("enrichex_overlap_results")
    output_dir.mkdir(exist_ok=True)
    logger.info(f"Output directory: {output_dir}")

    # Create gene sets
    logger.info("Creating brain cell-type gene sets...")
    gene_sets = create_example_gene_sets()

    # Save gene sets for future use
    gene_sets_path = output_dir / "brain_cell_types.json"
    gene_sets.save(gene_sets_path)
    logger.info(f"Saved gene sets to: {gene_sets_path}")

    # Get Alzheimer's GWAS genes
    query_genes = get_alzheimers_gwas_genes()
    logger.info(f"\nQuery: {len(query_genes)} Alzheimer's disease GWAS genes")
    logger.info(f"Gene sets: {len(gene_sets)} brain cell types")
    logger.info(f"Background: {len(gene_sets.background_genes)} genes")

    # Compute enrichment
    logger.info("\n" + "=" * 60)
    logger.info("Running overlap enrichment analysis...")
    logger.info("=" * 60)

    results_df = compute_overlap_enrichment_pandas(
        query_genes=query_genes,
        gene_set_collection=gene_sets,
        correction_method="benjamini-hochberg",
    )

    # Save all results
    results_path = output_dir / "ad_enrichment_all_results.tsv"
    results_df.to_csv(results_path, sep="\t", index=False)
    logger.info(f"\nSaved all results to: {results_path}")

    # Display significant results
    significant = results_df[results_df["significant"]].sort_values("p_adjusted")

    logger.info("\n" + "=" * 60)
    logger.info("SIGNIFICANT ENRICHMENTS (FDR < 0.05)")
    logger.info("=" * 60)

    if len(significant) == 0:
        logger.info("No significant enrichments found.")
    else:
        for _, row in significant.iterrows():
            logger.info(f"\n{row['gene_set_name']}:")
            logger.info(
                f"  Overlap: {row['n_overlap']}/{row['n_query']} query genes, "
                f"{row['n_overlap']}/{row['n_gene_set']} gene set genes"
            )
            logger.info(
                f"  Odds Ratio: {row['odds_ratio']:.2f} "
                f"(95% CI: [{row['ci_lower']:.2f}, {row['ci_upper']:.2f}])"
            )
            logger.info(
                f"  P-value: {row['p_value']:.2e} (adjusted: {row['p_adjusted']:.2e})"
            )

            # Show overlapping genes
            overlap_genes = row["overlap_genes"]
            if isinstance(overlap_genes, str):
                genes = [g.strip() for g in overlap_genes.split(",")]
                logger.info(f"  Genes: {', '.join(genes)}")

    # Display all results summary
    logger.info("\n" + "=" * 60)
    logger.info("ALL RESULTS SUMMARY")
    logger.info("=" * 60)
    logger.info(f"{'Cell Type':<25} {'OR':<8} {'P-adj':<12} {'Sig':<5}")
    logger.info("-" * 60)
    for _, row in results_df.sort_values("p_adjusted").iterrows():
        sig_marker = "***" if row["significant"] else ""
        logger.info(
            f"{row['gene_set_name']:<25} "
            f"{row['odds_ratio']:>7.2f} "
            f"{row['p_adjusted']:>11.2e} "
            f"{sig_marker:<5}"
        )

    # Biological interpretation
    logger.info("\n" + "=" * 60)
    logger.info("BIOLOGICAL INTERPRETATION")
    logger.info("=" * 60)

    if len(significant[significant["gene_set_name"] == "Microglia"]) > 0:
        logger.info(
            """
Microglia Enrichment Detected!

This finding supports the immune/neuroinflammatory hypothesis of Alzheimer's
disease. Key microglia genes in the overlap include:
- TREM2: Regulates microglial activation and Aβ clearance
- CD33: Inhibits Aβ phagocytosis
- MS4A6A/MS4A4A: Immune regulation

Clinical Implications:
- Microglia-targeted therapies may be beneficial
- Anti-inflammatory approaches warrant investigation
- Focus on innate immune pathway modulation

Next Steps:
1. Validate findings in independent GWAS cohorts
2. Test for rare variant burden in microglia genes (see burden_analysis_example.py)
3. Perform functional studies in microglia cell lines
4. Screen for microglia-modulating compounds
        """
        )
    else:
        logger.info(
            """
No significant microglia enrichment detected in this example.

This could indicate:
1. Sample size limitations
2. Disease heterogeneity
3. Involvement of other cell types
4. Need for rare variant analysis (not just common GWAS hits)

Consider:
- Running burden analysis on sequencing data
- Testing additional gene sets
- Combining with expression data
        """
        )

    logger.info("\n" + "=" * 60)
    logger.info("Analysis complete!")
    logger.info("=" * 60)
    logger.info("\nOutput files:")
    logger.info(f"  - Gene sets: {gene_sets_path}")
    logger.info(f"  - Results: {results_path}")


if __name__ == "__main__":
    main()
