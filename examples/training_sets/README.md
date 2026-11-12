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

Build these datasets first; `hvantk reprocess` builds most of them. Download ClinVar at a
dated release (`hvantk download clinvar --version YYYYMMDD`), so a rebuild reads the same
variants. On a Slurm cluster, the `ml-features` workflow skill (`skills/ml-features/`)
runs every step below as a job, with the input checks and a receipt.

| `--input` / argument | Dataset |
| --- | --- |
| `clinvar` | `clinvar:variants` |
| `dbnsfp` | `dbnsfp:variants` (4.9a, rankscores stored as floats) |
| `insider_variants` | `insider:variants` |
| `ptm` | `uniprot-ptm:sites` |
| `hgnc` | `hgnc:lookup` |
| gene index (`--spine` of `annotate prepare` and `compose`) | `hvantk annotate spine` output |
| `layer1` inputs | `gnomad-metrics:metrics`, `gevir:metrics`, `insider:interfaces`, and three expression summaries (`.h5ad`) |
| `--gene-list is_chd_gene=` | a text file of gene symbols, one per line or space-separated |

## Build

```bash
# 1. Heart expression summaries, one per matrix axis
hvantk expression summarize -m e-mtab-6814.h5ad \
    --group-by developmental_stage --filter-by organism_part=heart \
    --min-cells 1 -o heart_bulk_stage_summary.h5ad
hvantk expression summarize-ucsc -e farah/exprMatrix.tsv.gz -m farah/meta.tsv \
    --group-by Region --group-by Age --group-by major_cell_class --min-cells 10 \
    -o farah_region_age_class_summary.h5ad
hvantk expression summarize -m heart_cell_atlas_global.h5ad --group-by cell_type \
    -o hca_cell_type_summary.h5ad

# 2. One prepared table per gene axis (--hgnc for the axes not keyed by gene_id),
#    then the gene matrix
hvantk annotate prepare --spec clinvar_heart.yaml --axis gnomad_constraint \
    --input gnomad_metrics.ht --spine gene-index.ht --output prepared/gnomad_constraint.ht
hvantk annotate prepare --spec clinvar_heart.yaml --axis gevir \
    --input gevir.ht --spine gene-index.ht --output prepared/gevir.ht
hvantk annotate prepare --spec clinvar_heart.yaml --axis insider_interfaces \
    --input insider_interfaces.ht --spine gene-index.ht --hgnc hgnc.ht \
    --output prepared/insider_interfaces.ht
hvantk annotate prepare --spec clinvar_heart.yaml --axis heart_bulk \
    --input heart_bulk_stage_summary.h5ad --spine gene-index.ht --output prepared/heart_bulk.ht
hvantk annotate prepare --spec clinvar_heart.yaml --axis farah \
    --input farah_region_age_class_summary.h5ad --spine gene-index.ht --hgnc hgnc.ht \
    --output prepared/farah.ht
hvantk annotate prepare --spec clinvar_heart.yaml --axis hca \
    --input hca_cell_type_summary.h5ad --spine gene-index.ht --hgnc hgnc.ht \
    --output prepared/hca.ht
hvantk annotate compose --spec clinvar_heart.yaml --spine gene-index.ht \
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

`--legacy-names columns.tsv` takes a two-column TSV whose header row is `column` and
`legacy_name`, separated by a tab. It fills the dictionary's `legacy_name` column, for mapping to an
older dataset's column names.

## Outputs

| File | Holds |
| --- | --- |
| `clinvar_heart.ht` | The dataset as a Hail table, keyed by `locus` and `alleles` |
| `clinvar_heart.tsv.gz` | The same rows as a gzipped TSV, byte-identical across runs |
| `clinvar_heart.dictionary.tsv` | One row per column: legacy name, feature or not, tier, level, source, source version, description, non-null rate, true rate (booleans) |
| `clinvar_heart.manifest.json` | Row counts at each filter, label balance overall and per gene list, review stars, how genes were chosen, tier counts, excluded columns, gene-list mapping and input provenance |
