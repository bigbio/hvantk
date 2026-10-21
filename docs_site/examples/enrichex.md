# EnrichEx examples

Worked examples for gene-set overlap enrichment and case-control burden testing. The
[EnrichEx reference](../tools/enrichex.md) explains the statistics, every option and
the output columns.

## Overlap on the bundled gene sets

`examples/enrichex/synthetic_gene_sets.json` holds six brain cell-type sets with a
35-gene background. Overlap needs no Hail. From the repository root:

```bash
mkdir -p results
printf "APOE\nCD33\nMS4A6A\nTREM2\nGRIN2A\nNRGN\nAQP4\nGFAP\n" > query_genes.txt
hvantk enrichex overlap \
  -g query_genes.txt \
  -s examples/enrichex/synthetic_gene_sets.json \
  -o results/overlap_results.tsv \
  --generate-report
```

`results/overlap_results.tsv`, selected columns, values rounded:

```text
gene_set_name       n_overlap  odds_ratio  ci_lower  ci_upper  p_value  p_adjusted  overlap_genes
Excitatory Neurons  2          8.67        0.67      112.04    0.124    0.744       GRIN2A,NRGN
Oligodendrocytes    0          0.00        0.00      inf       0.315    0.805
Microglia           4          2.38        0.47      11.92     0.402    0.805       APOE,CD33,MS4A6A,TREM2
Endothelial         0          0.00        0.00      inf       0.553    0.830
Astrocytes          2          1.17        0.19      7.34      1.000    1.000       AQP4,GFAP
Inhibitory Neurons  0          0.00        0.00      inf       1.000    1.000
```

No set is significant with eight query genes and a 35-gene background. The table also
shows that the test is two-sided: Oligodendrocytes, with no overlap at all, ranks above
Microglia. The command also writes `results/enrichex_overlap.png` and
`results/enrichex_overlap_report.html`; copies of all three files are in
[`examples/enrichex/results/`](https://github.com/bigbio/hvantk/tree/main/examples/enrichex/results).

The same test from Python:

```python
from hvantk.algorithms.enrichex import (
    GeneSetCollection,
    compute_overlap_enrichment_pandas,
)

gene_sets = GeneSetCollection.load("examples/enrichex/synthetic_gene_sets.json")
query = ["APOE", "CD33", "MS4A6A", "TREM2", "GRIN2A", "NRGN", "AQP4", "GFAP"]

df = compute_overlap_enrichment_pandas(
    query, gene_sets, correction_method="benjamini-hochberg"
)
print(df[["gene_set_name", "n_overlap", "odds_ratio", "p_value", "p_adjusted"]])
```

For your own analysis, build the gene-set JSON with an explicit background (see
[Gene sets](../tools/genesets.md)): the background decides what the test means.

## Burden testing

Burden testing runs Hail on a cohort MatrixTable with a `GT` entry and, by default, the
row fields `SYMBOL`, `gnomad_af`, `cadd_phred` and `consequence`. The phenotype table is
keyed by the same sample IDs as the MatrixTable columns.

```bash
hvantk enrichex burden \
  -m cohort.mt \
  -p phenotypes.ht \
  -s gene_sets.json \
  --covariates PC1,PC2,PC3,age,sex \
  --max-af 0.001 \
  --variant-classes lof,missense_constrained,synonymous \
  --generate-report \
  -o results/burden/
```

This writes `burden_results.tsv` (all classes), one `burden_<class>.tsv` per variant
class, and the report to `results/burden/`. Without a cohort of your own, run
`synthetic_burden_demo.py` (below): it simulates a cohort with planted signal and runs
the same test through `BurdenPipeline`.

## Scripts

The scripts in
[`examples/enrichex/`](https://github.com/bigbio/hvantk/tree/main/examples/enrichex)
write to the current directory unless noted.

| Script | Needs Hail | Shows | Writes |
| --- | --- | --- | --- |
| `overlap_enrichment_example.py` | no | Overlap of Alzheimer's disease GWAS genes with brain cell-type sets, Python API | `enrichex_overlap_results/` |
| `create_gene_sets_example.py` | no | Gene-set JSON from a Seurat marker table, by hand, from GMT, and from a filtered DataFrame | `enrichex_gene_sets/` |
| `plot_example.py` | no | Dot and bar plots from `results/overlap_results.tsv` | `plots/synthetic_enrichex/` |
| `report_example.py` | no | HTML report from `results/overlap_results.tsv` | `reports/` (`--output` to change) |
| `burden_analysis_example.py` | yes | Binary and continuous burden tests on a random synthetic cohort, with Benjamini-Hochberg correction | `enrichex_burden_results/` |
| `synthetic_burden_demo.py` | yes | `BurdenPipeline` with variant classes on a cohort with planted signal; heatmap, volcano and forest plots; a type-I error check on the synonymous class | A temporary directory, deleted at exit |
| `chd_burden_run.py` | yes, except `--dry-run` | `BurdenPipeline` on a MatrixTable with custom field names and the phenotype in MatrixTable column fields | `--output-dir` |

Plots and reports need the `enrichex` extra (`pip install 'hvantk[enrichex]'`).
