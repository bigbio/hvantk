# QTL cascade walkthrough

This walkthrough runs the eQTL × pQTL cascade for one tissue, liver, from raw files to
the HTML report, reads the results, extends the run to more tissues, and then
colocalises two GWAS loci with `gwas-coloc`. Options, methods and output columns are on
the [QTL cascade](../tools/qtlcascade.md) reference page.

## Get the raw files

Both QTL sources are bring-your-own (see [Data sources](../guide/data-sources.md#qtl-data)).
Place them like this:

```text
data/
├── gtex_v11/liver_allpairs/              # GTEx v11 liver all-pairs parquet, file names starting "Liver."
└── fang_pqtl/
    └── Liver.allpairs_nobsGE72.txt.gz    # Fang et al. (2025) liver cis-pQTL, every pair
```

The builders read the tissue from the file name before the first dot, so
`Liver.allpairs_nobsGE72.txt.gz` is tissue `Liver`. The HGNC and gnomAD files are
downloaded in the next step.

Colocalisation needs every variant-gene pair around each gene, so it needs the GTEx
all-pairs files. With only the significant pairs (`Liver.v11.eQTLs.signif_pairs.parquet`)
you can still run the cascade; leave out the two all-pairs options of `run` to skip
coloc.

## Build the tables

```bash
# HGNC lookup, to map the pQTL gene symbols to Ensembl IDs
hvantk reprocess hgnc:lookup --raw-dir data/hgnc/ --output data/hgnc.ht

# eQTL, every variant-gene pair (p_threshold=0), so the table also serves coloc
hvantk reprocess gtex-eqtl:eqtls \
  --raw-dir data/gtex_v11/liver_allpairs/ \
  --output data/eqtl_liver_allpairs.ht \
  --plugin-arg source=gtex_v11 \
  --plugin-arg tissue=Liver \
  --plugin-arg p_threshold=0

# pQTL: the Fang file already holds every pair, so no p_threshold
hvantk reprocess pqtl:metrics \
  --raw-dir data/fang_pqtl/ \
  --output data/pqtl_liver.ht \
  --plugin-arg tissue=Liver \
  --plugin-arg hgnc_ht=data/hgnc.ht

# LOEUF overlay: gnomAD v2.1.1 constraint by gene
hvantk reprocess gnomad-metrics:metrics \
  --raw-dir data/gnomad_metrics/ \
  --output data/gnomad_metrics.ht
```

- The pQTL build logs `Mapped N of M pQTL gene symbol(s)`. A symbol left unmapped keeps
  the symbol as `gene_id` and never joins the eQTL table.
- Pass `tissue=` whenever the raw directory holds files of other tissues or other
  `.txt.gz` files: without it the pQTL builder imports every `.txt.gz` and `.tsv.gz`
  file it finds.
- `reprocess` runs a drift check after each build, which needs network access; add
  `--no-check-drift` on a node without it.
- All-pairs tables are large (the Fang liver file alone has about 34 million rows).
  Build them as cluster jobs ([Running on HPC](../guide/hpc-migration.md)).

## Check the plan

```bash
hvantk qtlcascade run \
  --eqtl-ht data/eqtl_liver_allpairs.ht \
  --pqtl-ht data/pqtl_liver.ht \
  --eqtl-allpairs data/eqtl_liver_allpairs.ht \
  --pqtl-allpairs data/pqtl_liver.ht \
  --constraint-ht data/gnomad_metrics.ht \
  --tissues Liver \
  -o results/qtl_liver \
  --dry-run
```

```text
======================================================================
QTL CASCADE PIPELINE — EXECUTION PLAN
======================================================================

Tissues:          Liver
eQTL table:       data/eqtl_liver_allpairs.ht
pQTL table:       data/pqtl_liver.ht
Coloc:            yes
Constraint:       data/gnomad_metrics.ht
Disease genes:    (none)
Output directory: results/qtl_liver

Stages:
  build_cascade
  run_coloc
  build_gene_summary
  generate_outputs

Total runs: 1
======================================================================
```

Each table is passed twice: as a cascade input, which `run` filters at p ≤ 5e-8, and as
a coloc input, which keeps every pair.

## Run the cascade

Run the same command without `--dry-run`. It ends with:

```text
Completed: 1 tissue(s)
  Liver: <n> genes, <m> colocalised

Outputs: results/qtl_liver
```

`<n>` counts every gene in the gene summary, whatever its class; `<m>` counts the genes
with H4 > 0.8. The files are laid out as in the reference
[quick start](../tools/qtlcascade.md#quick-start), under `results/qtl_liver/`.

## Read the results

Start with `results/qtl_liver/qtlcascade_report.html`: class counts, the genes with the
most concordant pairs, the coloc summary and the plots.

To list the genes with a complete cascade (at least one `eqtl_mediated` pair) whose eQTL
and pQTL colocalise:

```python
import pandas as pd

genes = pd.read_csv("results/qtl_liver/per_tissue/Liver/gene_summary.tsv", sep="\t")
hits = genes[genes["has_complete_cascade"] & (genes["coloc_max_h4"] > 0.8)]
columns = ["gene_id", "n_concordant", "n_discordant", "coloc_max_h4", "oe_lof_upper"]
print(hits.sort_values("coloc_max_h4", ascending=False)[columns].head(20))
```

A low `oe_lof_upper` (LOEUF) marks a gene that tolerates loss-of-function variants
poorly.

To see which hypothesis wins for each tested gene:

```python
coloc = pd.read_csv("results/qtl_liver/per_tissue/Liver/coloc_results.tsv", sep="\t")
coloc["best"] = coloc[["H0", "H1", "H2", "H3", "H4"]].idxmax(axis=1)
print(coloc["best"].value_counts())
```

Genes where H3 wins have an eQTL and a pQTL at the same variant but distinct causal
variants: the cascade there is an LD artifact. See
[Interpreting results](../tools/qtlcascade.md#interpreting-results) for the other
cases.

## Rerun one stage

To check whether the coloc calls depend on the window, colocalise the same genes with a
narrower one. The `awk` line selects what `run` tests, the genes with an
`eqtl_mediated` or `discordant` pair:

```bash
awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) col[$i] = i; next}
     $col["n_concordant"] + $col["n_discordant"] > 0 {print $col["gene_id"]}' \
  results/qtl_liver/per_tissue/Liver/gene_summary.tsv > results/coloc_genes.txt

hvantk qtlcascade coloc \
  --eqtl-allpairs data/eqtl_liver_allpairs.ht \
  --pqtl-allpairs data/pqtl_liver.ht \
  --cascade-genes results/coloc_genes.txt \
  --tissue Liver \
  --window-kb 250 \
  -o results/coloc_liver_250kb.tsv
```

`hvantk qtlcascade report` rebuilds the HTML report from the TSVs and plots if you need
it again ([Report](../tools/qtlcascade.md#report-report)).

## Add tissues

Liver, lung and thyroid carry the same names in GTEx and in the Fang data. Build the eQTL
and pQTL tables without `--plugin-arg tissue=`, from directories that hold only the QTL
files you want (the builders import every matching file they find), then run the tissues
together:

```bash
hvantk qtlcascade run \
  --eqtl-ht data/eqtl_allpairs.ht \
  --pqtl-ht data/pqtl.ht \
  --eqtl-allpairs data/eqtl_allpairs.ht \
  --pqtl-allpairs data/pqtl.ht \
  --constraint-ht data/gnomad_metrics.ht \
  --tissues Liver,Lung,Thyroid \
  -o results/qtl_three_tissues
```

Each tissue gets its own `per_tissue/<tissue>/` directory and plots,
`plots/cross_tissue_heatmap.png` compares the class counts across tissues, and
`qtlcascade_report.html` covers all three.

GTEx splits heart and colon into sub-tissues (`Heart_Left_Ventricle`,
`Heart_Atrial_Appendage`, `Colon_Transverse`, `Colon_Sigmoid`), so `--tissues Heart`
would find no eQTL rows. Pair such a tissue explicitly: build one table on each side and
run them without `--tissues`.

```bash
hvantk reprocess gtex-eqtl:eqtls \
  --raw-dir data/gtex_v11/heart_lv_allpairs/ \
  --output data/eqtl_heart_lv_allpairs.ht \
  --plugin-arg tissue=Heart_Left_Ventricle \
  --plugin-arg p_threshold=0

hvantk reprocess pqtl:metrics \
  --raw-dir data/fang_pqtl/ \
  --output data/pqtl_heart.ht \
  --plugin-arg tissue=Heart \
  --plugin-arg hgnc_ht=data/hgnc.ht

hvantk qtlcascade run \
  --eqtl-ht data/eqtl_heart_lv_allpairs.ht \
  --pqtl-ht data/pqtl_heart.ht \
  --eqtl-allpairs data/eqtl_heart_lv_allpairs.ht \
  --pqtl-allpairs data/pqtl_heart.ht \
  -o results/qtl_heart_lv
```

The outputs then go to `per_tissue/all/`, the coloc `tissue` column reads `unknown`, and
the log warns that the inputs carry a `tissue` field, which is harmless when each table
holds one tissue. To include such a pair in a multi-tissue run instead, give the GTEx
file, or a symlink to it, a name that starts with the Fang tissue (`Heart.`) before
building.

## Colocalise a GWAS locus with `gwas-coloc`

`gwas-coloc` does not use the tables above. It streams FinnGen and eQTL Catalogue
statistics on every run and saves the remote tabix indexes in the current directory, so
start it from a scratch directory on a node with network access.

### A confirmed effector: atrial fibrillation and MYOZ1

MYOZ1 at 10q22 is a literature-supported candidate effector gene for atrial fibrillation
and serves as the positive control. The command tests the FinnGen endpoint `I9_AF`
against GTEx heart atrial appendage eQTLs (`QTD000251`) and confirms MYOZ1
(`ENSG00000177791`):

```bash
hvantk qtlcascade gwas-coloc \
  --endpoint I9_AF --chrom 10 --lead 73600000 \
  --eqtl QTD000251 \
  --gene ENSG00000177791 \
  --gwas-n 261395 --eqtl-n 372 \
  -o results/myoz1
```

Output from one run (PP4 values are rounded to three decimals, and the remote data can
change between releases):

```text
Region chr10:73100000-74100000  |  GWAS min-p 2.7e-14  |  51 genes tested
Top effector: ENSG00000177791 (PP4=0.812)
Fine-map: credible sets GWAS=1 eQTL=1; coloc.susie PP4=0.649
VERDICT: CONFIRMED (ABF + fine-mapping agree)
Report: results/myoz1/report_I9_AF_10_73600000.json
```

MYOZ1 ranks first by ABF PP4, each trait has a credible set, and coloc.susie PP4 is
above 0.5, so the call is `CONFIRMED`.

### ABF only

```bash
hvantk qtlcascade gwas-coloc \
  --endpoint I9_AF --chrom 10 --lead 73600000 \
  --eqtl QTD000251 --no-fine-map \
  -o results/myoz1_abf
```

`--no-fine-map` needs neither sample sizes nor the LD reference, but the GWAS and eQTL
statistics are still streamed. The verdict stops at
`SUGGESTIVE (ABF only; fine-mapping not run)`.

### A refuted look-alike: septal defects and NSF at 17q21

```bash
hvantk qtlcascade gwas-coloc \
  --endpoint Q17_SEPTA_DEFEC --chrom 17 --lead 46890164 \
  --eqtl QTD000136 \
  --gene ENSG00000073969 \
  --gwas-n 412181 --eqtl-n 213 \
  -o results/nsf_17q21
```

```text
Top effector: ENSG00000073969 (PP4=0.810)
Fine-map: credible sets GWAS=0 eQTL=0; coloc.susie PP4=0.000
VERDICT: REFUTED (no fine-mappable signal — single-variant-ABF artifact)
```

By ABF, NSF (`ENSG00000073969`) looks like MYOZ1. Fine-mapping finds no credible set in
either trait, so coloc.susie PP4 is 0 and the call is `REFUTED`: the ABF hit is a
single-variant-ABF artifact, the kind a strong GWAS signal next to a weak eQTL produces.

### Read the report

```python
import json

with open("results/myoz1/report_I9_AF_10_73600000.json") as fh:
    report = json.load(fh)

print(report["verdict"])
print(report["results"]["top_effector"], report["results"]["goi_PP4_abf"])
print(report["fine_map"]["credible_sets_gwas"], report["fine_map"]["coloc_susie_PP4"])
```

`coloc_I9_AF_10_73600000.tsv` in the same directory ranks every tested gene by PP4. To
run the same analysis from Python, see
[GWAS colocalisation from Python](../tools/qtlcascade.md#gwas-colocalisation-from-python).

## See also

- [QTL cascade](../tools/qtlcascade.md): options, methods, output columns and the Python
  API.
- [Data sources](../guide/data-sources.md#qtl-data): where to get the eQTL and pQTL
  files.
- [Building datasets](../guide/usage.md#build-a-dataset): `hvantk reprocess`.
- [Running on HPC](../guide/hpc-migration.md): building the all-pairs tables as cluster
  jobs.
