# PTM example: ClinVar and gnomAD

This example builds the PTM sites table, classifies ClinVar variants, tests whether
pathogenic variants are enriched at PTM sites, compares predictor performance at PTM
and non-PTM variants, summarizes allele frequencies by PTM class, and renders a report.
[PTM](../tools/ptm.md) documents the method, options and outputs of each command.

## Before you start

You need hvantk with Hail ([Installation](../getting-started/installation.md)) and, for
the plots and the report, the `viz` extra (`poetry install --extras viz`).

| Input | How to get it |
| --- | --- |
| ClinVar table | `hvantk reprocess clinvar:variants --raw-dir data/clinvar/ --output data/clinvar.ht` |
| Variant table with allele frequencies (step 5) | Bring your own, for example a gnomAD sites release; hvantk has no gnomAD variant dataset |
| dbNSFP table (step 4) | Build it as described in [Data sources](../guide/data-sources.md) |
| UniProt PTM sites and the Ensembl GTF | Downloaded by `hvantk ptm build` |

## Step 1: Build the PTM sites table

```bash
hvantk ptm build \
  --output-dir data/ptm/ \
  --output-ht data/ptm/ptm_sites.ht
```

`build` downloads today's UniProt sites and the Ensembl GTF into `data/ptm/`, maps every
site to its GRCh38 codon and writes the table. It prints the sites mapped per source,
the overall mapping rate, how many proteins each transcript-resolution rule handled, and
the paths of the mapped TSV and the table.

UniProt has no versioned archive, so keep the downloaded
`data/ptm/uniprot-ptm-human-<date>.tsv`: passing it back with `--ptm-tsv` rebuilds the
same table later.

To include PeptideAtlas phosphosites, run these two commands instead of the one above,
passing the TSV path that the download prints:

```bash
hvantk download peptideatlas-phospho -o data/ptm/
hvantk ptm build \
  --output-dir data/ptm/ \
  --output-ht data/ptm/ptm_sites.ht \
  --peptideatlas-tsv data/ptm/peptideatlas-phospho-<build_date>-<build_id>.tsv
```

## Step 2: Classify ClinVar variants

```bash
hvantk ptm annotate \
  --variants-ht data/clinvar.ht \
  --ptm-ht data/ptm/ptm_sites.ht \
  -o data/clinvar_ptm.ht
```

The output is the ClinVar table plus `is_ptm_site`, `is_ptm_proximal`, `ptm_types`,
`ptm_distance` and `ptm_evidence`, and the command prints the number of site and
proximal variants. [How variants are classified](../tools/ptm.md#how-variants-are-classified)
gives the rules.

## Step 3: Pathogenic variants at PTM sites

```bash
hvantk ptm landscape \
  --clinvar-ht data/clinvar.ht \
  --ptm-ht data/ptm/ptm_sites.ht \
  -o results/landscape/ \
  --save-plots
```

`landscape` classifies the ClinVar table itself, so it does not need step 2. It writes
`results/landscape/landscape_summary.json` and, with `--save-plots`,
`landscape_summary.png`, `overlap_by_category.png` and `distance_distribution.png`. In
the JSON, `enrichment` holds the odds ratio of P/LP against B/LB variants for PTM (site
or proximal) against non-PTM positions, with its 95% confidence interval and Fisher's
exact p-value. `category_enrichment` repeats the test for each PTM category among the
PTM variants.

## Step 4: Predictor performance at PTM and non-PTM variants

Export the two variant lists, then run PS-ROC once on each:

```bash
hvantk ptm export-strata --annotated-ht data/clinvar_ptm.ht -o strata/

hvantk psroc \
  --variants strata/ptm_variants.txt \
  --clinvar-ht data/clinvar.ht \
  --dbnsfp-ht data/dbnsfp.ht \
  --scores "CADD_phred,REVEL_score,MetaLR_score" \
  --output-dir results/psroc_ptm/

hvantk psroc \
  --variants strata/non_ptm_variants.txt \
  --clinvar-ht data/clinvar.ht \
  --dbnsfp-ht data/dbnsfp.ht \
  --scores "CADD_phred,REVEL_score,MetaLR_score" \
  --output-dir results/psroc_non_ptm/
```

Compare the AUCs in `results/psroc_ptm/psroc_metrics.json` and
`results/psroc_non_ptm/psroc_metrics.json`. [PS-ROC](../tools/psroc.md) describes the
metrics and options.

## Step 5: Allele frequencies by PTM class

```bash
hvantk ptm population \
  --gnomad-ht data/gnomad_sites.ht \
  --ptm-ht data/ptm/ptm_sites.ht \
  --af-field AF \
  -o results/population/ \
  --save-plots
```

`--af-field` must name a top-level numeric field. gnomAD release tables keep allele
frequencies in the `freq` array, so copy the entry you want into a top-level `AF`
field first. A row with a missing value is skipped rather than an error. The non-PTM
class holds every other variant in the table; to compare against coding variants only,
filter the table before this step. The
command writes `population_summary.json`, with the counts and mean allele frequency of
each class, and, with `--save-plots`, `population_af.png`.

## Step 6: Report

```bash
hvantk ptm report \
  -o results/ptm_report.html \
  --landscape-json results/landscape/landscape_summary.json \
  --population-json results/population/population_summary.json \
  --title "PTM-variant report" \
  --description "ClinVar and gnomAD variants at UniProt PTM sites"
```

The page embeds its plots and also writes them as PNG files next to
`ptm_report.html`.

## See also

- [PTM](../tools/ptm.md): method, options and outputs of every `hvantk ptm` command, and
  the [Python API](../tools/ptm.md#python-api) for these steps
- [PTM constraint](../tools/ptm-constraint.md): compare PTM and non-PTM allele frequencies
  within the tissue or cell type where each gene is most expressed
- [PS-ROC](../tools/psroc.md): the predictor evaluation used in step 4
- [Data sources](../guide/data-sources.md): building the ClinVar and dbNSFP tables, and
  the PTM site downloads
