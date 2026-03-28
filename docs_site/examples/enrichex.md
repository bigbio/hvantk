# EnrichEx Examples

EnrichEx tests whether a gene list is enriched in curated gene sets (overlap) or whether cases carry excess rare variants in gene set genes (burden).

See the [EnrichEx reference](../tools/enrichex.md) for full CLI options and output format details.

## Overlap Enrichment

Test gene list enrichment using Fisher's exact test:

### CLI

```bash
hvantk enrichex overlap \
  -g my_genes.txt \
  -s gene_sets.json \
  -o results/overlap_results.tsv \
  --generate-report
```

### Python API

```python
from hvantk.enrichex import GeneSetCollection, compute_overlap_enrichment

query_genes = ["BRCA1", "TP53", "EGFR"]  # or load from my_genes.txt
gene_sets = GeneSetCollection.load("gene_sets.json")
results = compute_overlap_enrichment(
    query_genes=query_genes,
    gene_set_collection=gene_sets,
    correction_method="benjamini-hochberg",
)
```

## Burden Testing

Test rare variant burden in gene sets using Hail regression:

### CLI

```bash
hvantk enrichex burden \
  -m cohort.mt \
  -p phenotypes.ht \
  -s gene_sets.json \
  --covariates PC1,PC2,PC3,age,sex \
  --max-af 0.001 \
  --min-cadd 25 \
  -o results/burden_results.tsv \
  --generate-report
```

## Input Requirements

**Gene list** (overlap): plain text, one gene symbol per line.

**Gene sets** (both): JSON with `gene_sets` and optional `background_genes`:

```json
{
  "gene_sets": {
    "Microglia": {
      "name": "Microglia",
      "genes": ["APOE", "TREM2", "CD33", "MS4A6A"]
    }
  },
  "background_genes": ["APOE", "TREM2", "CD33", "..."]
}
```

**Cohort MatrixTable** (burden): row annotations `SYMBOL`, `gnomad_af`, `cadd_phred`; entry field `GT`.

**Phenotype table** (burden): sample IDs + phenotype field + optional covariates.

## Expected Outputs

- TSV results with p-values, odds ratios, confidence intervals, and significance flags
- PNG plots (dot plot for overlap, forest plot for burden)
- Self-contained HTML reports (with `--generate-report`)

## Gene Set Sources

- **Cell-type markers**: Lake et al. 2018, PanglaoDB, CellMarker
- **Pathway databases**: MSigDB, Gene Ontology, KEGG
- **GWAS genes**: GWAS Catalog

## Runnable Scripts

See the [`examples/enrichex/`](https://github.com/bigbio/hvantk/tree/main/examples/enrichex/) directory for Python API examples and pre-generated results.
