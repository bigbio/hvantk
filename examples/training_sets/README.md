# ClinVar x heart ML training set

`clinvar_heart.yaml` declares a labelled feature dataset for training variant
classifiers on heart genes. The rows are ClinVar GRCh38 missense SNVs labelled 1
(pathogenic, likely pathogenic) or 0 (benign, likely benign), with at least one review
star. Each row carries:

- dbNSFP 4.9a rankscores for 44 predictors and conservation scores;
- gnomAD allele frequencies;
- INSIDER interface and UniProt PTM flags;
- gene constraint, interaction-partner counts and gene structure;
- heart expression: bulk developmental, fetal single-cell and adult single-cell.

Every feature column has a leakage tier, L1 to L4. Training on `tier <= L1` gives a
model that saw nothing derived from the label or from allele frequency. See
[the command's page](../../docs_site/tools/training-set.md) for the tier rule.

## Inputs

Build these datasets first; `hvantk reprocess` builds most of them:

| `--input` / argument | Dataset |
| --- | --- |
| `clinvar` | `clinvar:variants` |
| `dbnsfp` | `dbnsfp:variants` (4.9a, rankscores stored as floats) |
| `insider_variants` | `insider:variants` |
| `ptm` | `uniprot-ptm:sites` |
| `hgnc` | `hgnc:lookup` |
| `--spine` | `hvantk annotate spine` output |
| `layer1` inputs | `gnomad-metrics:metrics`, `gevir:metrics`, `insider:interfaces`, and three expression summaries (`.h5ad`) |
| `--gene-list is_chd_gene=` | a text file of gene symbols, one per line or space-separated |

## Build

```bash
# 1. Bulk heart expression by developmental stage
hvantk expression summarize -m e-mtab-6814.h5ad \
    --group-by developmental_stage --filter-by organism_part=heart \
    --min-cells 1 -o heart_bulk_stage_summary.h5ad

# 2. One prepared table per gene axis, then the gene matrix
for axis in gnomad_constraint gevir insider_interfaces heart_bulk farah hca; do
  hvantk annotate prepare --spec clinvar_heart.yaml --axis "$axis" \
      --input "<built source for $axis>" --spine spine.ht --hgnc hgnc.ht \
      --output "prepared/$axis.ht"
done
hvantk annotate compose --spec clinvar_heart.yaml --spine spine.ht \
    --prepared gnomad_constraint=prepared/gnomad_constraint.ht \
    --prepared gevir=prepared/gevir.ht \
    --prepared insider_interfaces=prepared/insider_interfaces.ht \
    --prepared heart_bulk=prepared/heart_bulk.ht \
    --prepared farah=prepared/farah.ht --prepared hca=prepared/hca.ht \
    --output gene_matrix.ht

# 3. The training set
hvantk training-set build --spec clinvar_heart.yaml \
    --input clinvar=clinvar.ht --input dbnsfp=dbnsfp.ht \
    --input insider_variants=insider_variants.ht --input ptm=ptm_sites.ht \
    --input hgnc=hgnc.ht --gene-matrix gene_matrix.ht \
    --gene-list is_chd_gene=chd_genes.txt \
    --output out/clinvar_heart
```

`--legacy-names columns.tsv` takes a two-column TSV with the header
`column	legacy_name`. It fills the dictionary's `legacy_name` column, for mapping to an
older dataset's column names.

## Outputs

| File | Holds |
| --- | --- |
| `clinvar_heart.ht` | The dataset as a Hail table, keyed by `locus` and `alleles` |
| `clinvar_heart.tsv.gz` | The same rows as a gzipped TSV, byte-identical across runs |
| `clinvar_heart.dictionary.tsv` | One row per column: legacy name, feature or not, tier, level, source, source version, description, non-null rate, true rate (booleans) |
| `clinvar_heart.manifest.json` | Row counts at each filter, label balance overall and per gene list, review stars, how genes were chosen, tier counts, excluded columns, gene-list mapping and input provenance |
