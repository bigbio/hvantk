#!/usr/bin/env python3
"""
Create Gene Sets from Marker Files

This example demonstrates how to create gene set collections from various
input formats commonly used in single-cell and bulk RNA-seq analysis.

Common sources:
- Seurat FindAllMarkers() output
- Scanpy rank_genes_groups() output
- DESeq2 differential expression results
- Manual curation from literature
"""

import logging
from pathlib import Path
import pandas as pd
from hvantk.enrichex import (
    GeneSet,
    GeneSetCollection,
    load_marker_genes,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def example_seurat_format():
    """
    Example: Load gene sets from Seurat FindAllMarkers output.

    Seurat marker format (TSV):
        gene        cluster     p_val       avg_log2FC  pct.1   pct.2
        TREM2       Microglia   1.23e-50    2.5         0.95    0.10
        CD33        Microglia   4.56e-45    2.1         0.88    0.15
        ...
    """
    logger.info("\n" + "=" * 60)
    logger.info("Example 1: Seurat FindAllMarkers Format")
    logger.info("=" * 60)

    # Create example data
    seurat_data = pd.DataFrame(
        {
            "gene": [
                "TREM2",
                "CD33",
                "MS4A6A",
                "TYROBP",
                "CSF1R",
                "SLC17A7",
                "GRIN2A",
                "GRIN2B",
                "CAMK2A",
                "NRGN",
                "GAD1",
                "GAD2",
                "PVALB",
                "SST",
                "VIP",
            ],
            "cluster": [
                "Microglia",
                "Microglia",
                "Microglia",
                "Microglia",
                "Microglia",
                "Excitatory",
                "Excitatory",
                "Excitatory",
                "Excitatory",
                "Excitatory",
                "Inhibitory",
                "Inhibitory",
                "Inhibitory",
                "Inhibitory",
                "Inhibitory",
            ],
            "p_val": [
                1e-50,
                1e-45,
                1e-40,
                1e-38,
                1e-35,
                1e-48,
                1e-42,
                1e-40,
                1e-38,
                1e-35,
                1e-45,
                1e-44,
                1e-40,
                1e-38,
                1e-35,
            ],
            "avg_log2FC": [
                2.5,
                2.1,
                2.0,
                1.9,
                1.8,
                3.0,
                2.8,
                2.5,
                2.3,
                2.0,
                2.8,
                2.7,
                2.4,
                2.2,
                2.0,
            ],
            "pct.1": [
                0.95,
                0.88,
                0.85,
                0.82,
                0.80,
                0.98,
                0.95,
                0.92,
                0.90,
                0.85,
                0.96,
                0.95,
                0.90,
                0.88,
                0.85,
            ],
            "pct.2": [
                0.10,
                0.15,
                0.18,
                0.20,
                0.22,
                0.05,
                0.08,
                0.10,
                0.12,
                0.15,
                0.08,
                0.10,
                0.12,
                0.15,
                0.18,
            ],
        }
    )

    # Save to TSV
    output_dir = Path("enrichex_gene_sets")
    output_dir.mkdir(exist_ok=True)
    seurat_path = output_dir / "seurat_markers.tsv"
    seurat_data.to_csv(seurat_path, sep="\t", index=False)
    logger.info(f"Created example Seurat file: {seurat_path}")

    # Load as gene sets
    gene_sets = load_marker_genes(
        marker_file=str(seurat_path),
        gene_column="gene",
        cluster_column="cluster",
    )

    logger.info(f"\nLoaded {len(gene_sets)} gene sets:")
    for gs in gene_sets:
        logger.info(f"  {gs.name}: {gs.n_genes} genes")

    # Save as JSON for EnrichEx
    json_path = output_dir / "seurat_gene_sets.json"
    gene_sets.save(json_path)
    logger.info(f"\nSaved gene sets to: {json_path}")

    return gene_sets


def example_manual_creation():
    """
    Example: Manually create gene sets from literature or databases.

    Use this when you have curated gene lists from:
    - Literature review
    - MSigDB pathways
    - GO terms
    - KEGG pathways
    - Custom functional groupings
    """
    logger.info("\n" + "=" * 60)
    logger.info("Example 2: Manual Gene Set Creation")
    logger.info("=" * 60)

    # Define gene sets manually
    ad_immune = GeneSet(
        name="AD_Immune_Genes",
        genes={
            "TREM2",
            "CD33",
            "MS4A6A",
            "MS4A4A",
            "TYROBP",
            "CR1",
            "CLU",
            "ABCA7",
            "INPP5D",
            "HLA-DRB1",
        },
        source="Lambert et al. Nat Genet 2013; Jansen et al. Nat Genet 2019",
        metadata={"category": "immune", "disease": "Alzheimers", "evidence": "GWAS"},
    )

    ad_lipid = GeneSet(
        name="AD_Lipid_Genes",
        genes={
            "APOE",
            "APOC1",
            "ABCA7",
            "CLU",
            "SORL1",
            "ABCA1",
            "LIPC",
            "CETP",
            "APOB",
            "LDLR",
        },
        source="Lambert et al. Nat Genet 2013",
        metadata={
            "category": "lipid_metabolism",
            "disease": "Alzheimers",
            "evidence": "GWAS",
        },
    )

    ad_endocytosis = GeneSet(
        name="AD_Endocytosis_Genes",
        genes={
            "PICALM",
            "BIN1",
            "CD2AP",
            "EPHA1",
            "RIN3",
            "SORL1",
            "PTK2B",
            "SH3KBP1",
            "GRB2",
            "AP2A2",
        },
        source="Lambert et al. Nat Genet 2013",
        metadata={
            "category": "endocytosis",
            "disease": "Alzheimers",
            "evidence": "GWAS",
        },
    )

    # Define background (all genes of interest)
    all_genes = set()
    for gs in [ad_immune, ad_lipid, ad_endocytosis]:
        all_genes.update(gs.genes)

    # Add additional background genes
    all_genes.update(
        [
            "ACTB",
            "GAPDH",
            "TUBB",
            "APP",
            "MAPT",
            "SNCA",
            "PSEN1",
            "PSEN2",
            "MEF2C",
            "CASS4",
            "FERMT2",
        ]
    )

    # Create collection
    gene_sets = GeneSetCollection(
        gene_sets={
            "AD_Immune": ad_immune,
            "AD_Lipid": ad_lipid,
            "AD_Endocytosis": ad_endocytosis,
        },
        background_genes=all_genes,
        source_description="Alzheimer's disease GWAS genes grouped by function",
        metadata={
            "references": [
                "Lambert et al. Nat Genet 2013 PMID:24162737",
                "Jansen et al. Nat Genet 2019 PMID:30617256",
            ],
            "date_created": "2026-01-28",
        },
    )

    # Save
    output_dir = Path("enrichex_gene_sets")
    output_dir.mkdir(exist_ok=True)
    json_path = output_dir / "ad_functional_gene_sets.json"
    gene_sets.save(json_path)

    logger.info(f"\nCreated {len(gene_sets)} functional gene sets:")
    for gs in gene_sets:
        logger.info(f"  {gs.name}: {gs.n_genes} genes")

    logger.info(f"\nSaved gene sets to: {json_path}")

    return gene_sets


def example_msigdb_format():
    """
    Example: Convert MSigDB GMT format to EnrichEx JSON.

    MSigDB GMT format:
        HALLMARK_INFLAMMATORY_RESPONSE   http://...   GENE1  GENE2  GENE3  ...
        HALLMARK_APOPTOSIS                http://...   GENE4  GENE5  GENE6  ...

    Download from: https://www.gsea-msigdb.org/gsea/msigdb/
    """
    logger.info("\n" + "=" * 60)
    logger.info("Example 3: MSigDB GMT Format")
    logger.info("=" * 60)

    # Create example GMT data
    output_dir = Path("enrichex_gene_sets")
    output_dir.mkdir(exist_ok=True)
    gmt_path = output_dir / "example.gmt"

    gmt_content = """HALLMARK_INFLAMMATORY_RESPONSE\thttp://www.gsea-msigdb.org/gsea/msigdb\tTREM2\tCD33\tC1QA\tIL1B\tTNF\tIL6
HALLMARK_APOPTOSIS\thttp://www.gsea-msigdb.org/gsea/msigdb\tBCL2\tBCL2L1\tBAX\tBID\tCASP3\tCASP9
HALLMARK_OXIDATIVE_PHOSPHORYLATION\thttp://www.gsea-msigdb.org/gsea/msigdb\tMT-ND1\tMT-ND2\tMT-CO1\tMT-CO2\tMT-ATP6
"""

    with open(gmt_path, "w") as f:
        f.write(gmt_content)

    logger.info(f"Created example GMT file: {gmt_path}")

    # Parse GMT file
    gene_sets_dict = {}
    all_genes = set()

    with open(gmt_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            url = parts[1]
            genes = set(parts[2:])

            gene_sets_dict[name] = GeneSet(
                name=name, genes=genes, source="MSigDB", metadata={"url": url}
            )
            all_genes.update(genes)

    # Create collection
    gene_sets = GeneSetCollection(
        gene_sets=gene_sets_dict,
        background_genes=all_genes,
        source_description="MSigDB Hallmark Gene Sets",
        metadata={
            "source": "MSigDB",
            "url": "https://www.gsea-msigdb.org/gsea/msigdb/",
        },
    )

    # Save
    json_path = output_dir / "msigdb_gene_sets.json"
    gene_sets.save(json_path)

    logger.info(f"\nParsed {len(gene_sets)} MSigDB gene sets:")
    for gs in gene_sets:
        logger.info(f"  {gs.name}: {gs.n_genes} genes")

    logger.info(f"\nSaved gene sets to: {json_path}")

    return gene_sets


def example_from_dataframe():
    """
    Example: Create gene sets from pandas DataFrame.

    Use this for custom data manipulation before creating gene sets.
    """
    logger.info("\n" + "=" * 60)
    logger.info("Example 4: From pandas DataFrame (Custom)")
    logger.info("=" * 60)

    # Create example marker data
    marker_data = pd.DataFrame(
        {
            "gene": [
                "TREM2",
                "CD33",
                "MS4A6A",
                "SLC17A7",
                "GRIN2A",
                "GAD1",
                "GAD2",
                "GFAP",
                "AQP4",
                "MBP",
            ],
            "cell_type": [
                "Microglia",
                "Microglia",
                "Microglia",
                "Excitatory",
                "Excitatory",
                "Inhibitory",
                "Inhibitory",
                "Astrocyte",
                "Astrocyte",
                "Oligodendrocyte",
            ],
            "log2fc": [3.2, 2.8, 2.5, 3.5, 3.1, 3.0, 2.9, 3.3, 3.0, 3.4],
            "p_value": [
                1e-50,
                1e-45,
                1e-40,
                1e-48,
                1e-42,
                1e-45,
                1e-44,
                1e-46,
                1e-43,
                1e-47,
            ],
            "pct_in": [0.95, 0.88, 0.85, 0.98, 0.95, 0.96, 0.95, 0.97, 0.93, 0.98],
            "pct_out": [0.10, 0.15, 0.18, 0.05, 0.08, 0.08, 0.10, 0.06, 0.09, 0.05],
        }
    )

    logger.info("Example marker data:")
    logger.info(marker_data.to_string())

    # Apply custom filters
    logger.info("\nApplying filters:")
    logger.info("  - log2FC > 2.0")
    logger.info("  - p_value < 1e-40")
    logger.info("  - pct_in > 0.90")

    filtered = marker_data[
        (marker_data["log2fc"] > 2.0)
        & (marker_data["p_value"] < 1e-40)
        & (marker_data["pct_in"] > 0.90)
    ]

    logger.info(f"\nFiltered: {len(filtered)}/{len(marker_data)} genes")

    # Group by cell type
    gene_sets_dict = {}
    all_genes = set(filtered["gene"])

    for cell_type in filtered["cell_type"].unique():
        genes = set(filtered[filtered["cell_type"] == cell_type]["gene"])
        gene_sets_dict[cell_type] = GeneSet(
            name=cell_type,
            genes=genes,
            source="Custom filtered markers",
            metadata={
                "log2fc_threshold": 2.0,
                "p_value_threshold": 1e-40,
                "pct_in_threshold": 0.90,
            },
        )

    # Create collection
    gene_sets = GeneSetCollection(
        gene_sets=gene_sets_dict,
        background_genes=all_genes,
        source_description="Stringently filtered cell-type markers",
    )

    # Save
    output_dir = Path("enrichex_gene_sets")
    output_dir.mkdir(exist_ok=True)
    json_path = output_dir / "filtered_gene_sets.json"
    gene_sets.save(json_path)

    logger.info(f"\nCreated {len(gene_sets)} gene sets:")
    for gs in gene_sets:
        logger.info(f"  {gs.name}: {gs.n_genes} genes")

    logger.info(f"\nSaved gene sets to: {json_path}")

    return gene_sets


def main():
    """Run all examples."""

    logger.info("=" * 60)
    logger.info("EnrichEx: Gene Set Creation Examples")
    logger.info("=" * 60)

    # Run examples
    example_seurat_format()
    example_manual_creation()
    example_msigdb_format()
    example_from_dataframe()

    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("SUMMARY")
    logger.info("=" * 60)

    logger.info(
        """
Gene sets created in multiple formats:
1. From Seurat FindAllMarkers output
2. Manual curation from literature
3. MSigDB GMT format conversion
4. Custom pandas DataFrame filtering

Output directory: enrichex_gene_sets/
  - seurat_gene_sets.json
  - ad_functional_gene_sets.json
  - msigdb_gene_sets.json
  - filtered_gene_sets.json

Use these gene sets with EnrichEx:

# Overlap enrichment
hvantk enrichex overlap \\
  -g my_genes.txt \\
  -s enrichex_gene_sets/seurat_gene_sets.json \\
  -o enrichment_results.tsv

# Burden testing
hvantk enrichex burden \\
  -m cohort.mt \\
  -p phenotypes.ht \\
  -s enrichex_gene_sets/seurat_gene_sets.json \\
  -o burden_results.tsv

# Python API
from hvantk.enrichex import GeneSetCollection
gene_sets = GeneSetCollection.load("enrichex_gene_sets/seurat_gene_sets.json")
    """
    )

    logger.info("\n" + "=" * 60)
    logger.info("BEST PRACTICES")
    logger.info("=" * 60)

    logger.info(
        """
1. Gene Set Size:
   - Minimum: 20 genes (for statistical power)
   - Maximum: 500 genes (avoid overly broad sets)
   - Optimal: 50-200 genes per set

2. Background Definition:
   - Use all genes expressed in your tissue/cell type
   - Match to your variant calling strategy
   - Typically 15,000-20,000 protein-coding genes

3. Gene Identifiers:
   - Use consistent ID system (symbols vs Ensembl IDs)
   - Ensure query genes and gene sets match
   - Check for gene name updates (e.g., MARCH1 → MARCHF1)

4. Quality Control:
   - Remove low-confidence markers (high p-values)
   - Filter by effect size (log2FC > 1.0)
   - Check for cell-type specificity (pct.in vs pct.out)

5. Source Documentation:
   - Always cite original source
   - Include PMID or URL in metadata
   - Note filtering criteria used

6. Version Control:
   - Date your gene sets
   - Track which version was used in analyses
   - Archive gene sets with results
    """
    )

    logger.info("\n" + "=" * 60)
    logger.info("ADDITIONAL RESOURCES")
    logger.info("=" * 60)

    logger.info(
        """
Gene Set Databases:
- MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/
- Gene Ontology: http://geneontology.org/
- KEGG: https://www.genome.jp/kegg/pathway.html
- Reactome: https://reactome.org/
- PanglaoDB (cell types): https://panglaodb.se/
- CellMarker: http://bio-bigdata.hrbmu.edu.cn/CellMarker/

Single-Cell Marker Tools:
- Seurat: https://satijalab.org/seurat/
- Scanpy: https://scanpy.readthedocs.io/
- scRNA-tools: https://www.scrna-tools.org/
    """
    )


if __name__ == "__main__":
    main()
