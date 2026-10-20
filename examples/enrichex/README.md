# EnrichEx examples

Scripts for gene-set overlap enrichment and case-control burden testing with
`hvantk enrichex`. The [EnrichEx reference](https://bigbio.github.io/hvantk/tools/enrichex/)
explains the statistics, options and output columns; the
[EnrichEx examples page](https://bigbio.github.io/hvantk/examples/enrichex/) walks through
the overlap run below.

## Files

```text
examples/enrichex/
├── synthetic_gene_sets.json        # six brain cell-type gene sets, 35-gene background
├── overlap_enrichment_example.py   # overlap, Python API (no Hail)
├── create_gene_sets_example.py     # build gene-set JSON files (no Hail)
├── plot_example.py                 # plots from results/overlap_results.tsv (no Hail)
├── report_example.py               # HTML report from results/overlap_results.tsv (no Hail)
├── burden_analysis_example.py      # burden test on a random synthetic cohort (Hail)
├── synthetic_burden_demo.py        # BurdenPipeline on a cohort with planted signal (Hail)
├── chd_burden_run.py               # BurdenPipeline with custom MatrixTable field names (Hail)
└── results/                        # sample outputs read by plot_example.py and report_example.py
    ├── overlap_results.tsv
    ├── enrichex_overlap.png
    ├── enrichex_overlap_report.html
    ├── burden_results.tsv
    ├── enrichex_burden.png
    └── enrichex_burden_report.html
```

The `overlap_*` files in `results/` are the output of the overlap command below. The
`burden_*` files come from an earlier synthetic burden run; a current `burden` run adds
the `gene_set_size`, `n_genes_found` and `gene_coverage_pct` columns.

## Overlap from the command line

Overlap needs no Hail. From this directory:

```bash
mkdir -p out
printf "APOE\nCD33\nMS4A6A\nTREM2\nGRIN2A\nNRGN\nAQP4\nGFAP\n" > out/query_genes.txt
hvantk enrichex overlap \
  -g out/query_genes.txt \
  -s synthetic_gene_sets.json \
  -o out/overlap_results.tsv \
  --generate-report
```

This writes `out/overlap_results.tsv`, `out/enrichex_overlap.png` and
`out/enrichex_overlap_report.html`, the same tables as in `results/`. No gene set is
significant: eight query genes against a 35-gene background give the test little power.

## Burden from the command line

Burden needs Hail and a cohort MatrixTable with a `GT` entry and, by default, the row
fields `SYMBOL`, `gnomad_af`, `cadd_phred` and `consequence`, plus a phenotype table
keyed by the MatrixTable's sample IDs:

```bash
hvantk enrichex burden \
  -m cohort.mt \
  -p phenotypes.ht \
  -s synthetic_gene_sets.json \
  --covariates PC1,PC2,PC3,age,sex \
  --max-af 0.001 \
  -o out/burden/ \
  --generate-report
```

## Scripts

Run them with hvantk importable (installed, or `PYTHONPATH` set to the repository root).
They write to the current directory unless noted.

| Script | Needs Hail | Writes |
| --- | --- | --- |
| `python overlap_enrichment_example.py` | no | `enrichex_overlap_results/brain_cell_types.json`, `ad_enrichment_all_results.tsv` |
| `python create_gene_sets_example.py` | no | `enrichex_gene_sets/` (Seurat, hand-made, GMT and filtered collections) |
| `python plot_example.py` | no | `plots/synthetic_enrichex/enrichex_dotplot.png`, `enrichex_barplot.png` |
| `python report_example.py` | no | `reports/synthetic_enrichex_report.html` and `reports/enrichex_overlap.png` (`--output` to change) |
| `python burden_analysis_example.py` | yes | `enrichex_burden_results/gene_sets.json`, `burden_binary_results.tsv`, `burden_continuous_results.tsv` |
| `python synthetic_burden_demo.py` | yes | A temporary directory, deleted at exit; the log lists the files |
| `python chd_burden_run.py --cohort-mt ... --gene-sets NAME:PATH --dry-run` | only without `--dry-run` | `--output-dir` (default `results/celltype_burden`) |

`burden_analysis_example.py` uses random genotypes and phenotypes, so no gene set should
be significant. `synthetic_burden_demo.py` plants signal in two gene sets and checks that
the synonymous class stays calibrated. Plots and reports need the `enrichex` extra
(`pip install 'hvantk[enrichex]'`).
