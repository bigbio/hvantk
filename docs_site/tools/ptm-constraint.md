# PTM Constraint: Stratified Allele-Frequency Depletion at PTM Codons

`hvantk ptm constraint` is a stratified depletion analysis that compares gnomAD
allele-frequency distributions between PTM-proximal and non-PTM variants,
grouped by tissue, cell type, or any categorical field found in an expression
dataset. It is the CLI incarnation of the multi-dataset EDA consolidated in
`local/notebooks/ptm-eda/PTM-EDA_Consolidated_Report.md` (Notebooks E–I).

> **Not a per-variant scorer.** For per-site PTM flags use `hvantk ptm annotate`.
> This command produces *stratified group-level* statistics.

## What it does

1. Reads a PTM-annotated variant Hail Table (output of `hvantk ptm annotate`).
2. Loads an expression dataset from one of three backends:
   - Hail MatrixTable (GTEx-style bulk RNA/proteomics)
   - AnnData `.h5ad` (scverse / single-cell)
   - Pre-aggregated wide tabular file (`.parquet`/`.pkl`/`.tsv`)
3. Aggregates expression to a per-group median (or mean / median-nonzero).
4. Computes per-gene τ (Yanai 2005 tissue specificity) via the maintained
   [`tspex`](https://pypi.org/project/tspex/) package.
5. Runs five statistical tests mirroring the EDA:
   - per-group ranking (Mann-Whitney)
   - τ quartile stratification
   - τ × LOEUF 2×2 factorial
   - PTM category × group heatmap
   - within-gene paired Wilcoxon

Outputs: five TSVs, four PNG panels, one self-contained HTML report, and a
`summary.json` for machine consumption.

## Quick start

```bash
# GTEx tissues (Hail MT backend)
hvantk ptm constraint \
  --variants-ht clinvar_ptm.ht \
  --expression-source hail-mt \
  --expression-path /data/GTEX_v7_TPMs.mt \
  --grouping SMTSD \
  --output-dir results/ptm-gtex/

# Farah developmental scRNA-seq (AnnData backend)
hvantk ptm constraint \
  --variants-ht clinvar_ptm.ht \
  --expression-source anndata \
  --expression-path /data/farah_2024.h5ad \
  --grouping major_cell_class \
  --output-dir results/ptm-farah/

# Pre-computed gene-by-group matrix (tabular backend)
hvantk ptm constraint \
  --variants-ht clinvar_ptm.ht \
  --expression-source tabular \
  --expression-path /data/gene_celltype.tsv \
  --grouping cell_class \
  --output-dir results/ptm-custom/
```

## Inputs

| Flag | Required | Description |
|------|:--------:|-------------|
| `--variants-ht` | ✅ | PTM-annotated variant Hail Table (from `hvantk ptm annotate`). Must contain `is_ptm_site`, `is_ptm_proximal`, `ptm_types`, the configured gene / AF / LOEUF / label fields. |
| `--expression-source` | ✅ | One of `hail-mt`, `anndata`, `tabular`. |
| `--expression-path` | ✅ | Path to the expression data. |
| `--grouping` | ✅ | Metadata field to stratify by. |
| `--output-dir` | ✅ | Output directory. |
| `--label-filter` |  | `TN` (default, benign) / `TP` (pathogenic) / `all`. |
| `--label-field` |  | Column in the variants HT holding the label. Default `rf_label`. |
| `--gene-field` |  | Variant gene identifier column. Default `gene_symbol`. |
| `--af-field` |  | Allele-frequency column. Default `gnomad_af_genomes`. |
| `--loeuf-field` |  | LOEUF column. Default `loeuf`. |
| `--ptm-category-field` |  | Field holding PTM category set. Default `ptm_types`. |
| `--gene-id-mapping` |  | Optional TSV mapping expression gene IDs → symbols. |
| `--expression-metric` |  | `median` (default) / `mean` / `median_nonzero`. |
| `--min-cells-per-group` |  | Drop small groups. Default 50. |
| `--min-variants-per-group` |  | Skip groups with fewer PTM + non-PTM variants. Default 20. |
| `--expressed-threshold` |  | Minimum per-gene max expression. Default 1.0. |

## Outputs

```
<output-dir>/
├── summary.json                     # top-line effect sizes + config echo
├── per_group_ranking.tsv            # Mann-Whitney per group, sorted
├── tau_quartile.tsv                 # τ-quartile stratified
├── loeuf_factorial.tsv              # 2x2 factorial (τ bin × LOEUF bin)
├── category_group_heatmap.tsv       # PTM category × group, long-form
├── within_gene_paired.tsv           # per-gene median AF PTM vs non-PTM
├── plots/
│   ├── per_group_ranking.png
│   ├── tau_quartile.png
│   ├── loeuf_factorial.png
│   └── category_group_heatmap.png
└── report.html                      # self-contained, base64-embedded plots
```

## Decisions and caveats

- **τ is always computed in pandas** via `tspex` — even for Hail MT input, a
  gene × group cache is materialised first. Keeps metric computation consistent
  across backends.
- **Primary group** per gene is the `argmax` of aggregate expression across
  groups. Genes expressed broadly lose information here; this is a deliberate
  trade-off, matching the EDA.
- **Gene identifier alignment** is the user's responsibility. If the expression
  matrix uses Ensembl IDs but the variants HT uses symbols, pass
  `--gene-id-mapping <tsv>`.
- **No trinucleotide-context control** — PTM codons (Ser/Thr/Tyr/Lys) may have
  different baseline mutation rates. This CLI does not correct for that; the
  synonymous-variant control at PTM codons is a separate (future) analysis.
- **Large Hail MTs** are aggregated via `group_cols_by → aggregate → entries →
  export TSV → pandas` to avoid the JVM-heap hit of `mt.to_pandas()`.

## Related

- [`hvantk ptm annotate`](ptm.md) — prerequisite per-variant PTM flagging.
- [`hvantk ptm landscape`](ptm.md) — ClinVar-based enrichment (complementary).
- [`hvantk ptm population`](ptm.md) — genome-wide gnomAD AF at PTM sites.
