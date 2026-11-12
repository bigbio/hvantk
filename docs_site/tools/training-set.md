# ML training sets

`hvantk training-set build` turns ClinVar, variant-level scores and the
[gene annotation matrix](annotate.md) into one labelled feature table for training
variant classifiers. Every feature column carries a **leakage tier**, so a model can be
trained on only the columns that saw nothing derived from the label or from allele
frequency.

## When to use it

Use it when you need ClinVar variants labelled pathogenic/benign together with predictor
scores, allele frequencies, protein-site flags and gene-level features. Typical uses are
training a classifier or benchmarking predictors without circularity. One feature-spec
YAML declares the whole dataset. The same file drives
[`annotate prepare` / `compose`](annotate.md) for the gene columns, so the gene layer and
the variant layer cannot drift apart.

## Inputs

- A feature spec (`--spec`) with `label`, `variant`, `gene` and `layer1` sections; see
  [The spec](#the-spec).
- One built dataset per input the spec names (`--input NAME=PATH`): the label source
  (`clinvar:variants`), each variant source (e.g. `dbnsfp:variants`,
  `insider:variants`, `uniprot-ptm:sites`) and the gene section's `hgnc:lookup`.
- The gene matrix (`--gene-matrix`): the `annotate compose` output built from the same
  spec, with the manifest compose wrote next to it.
- Optionally, gene lists (`--gene-list COLUMN=FILE`, symbols separated by whitespace) and
  a legacy-name map (`--legacy-names`).

## Worked example

[`examples/training_sets/`](https://github.com/bigbio/hvantk/tree/main/examples/training_sets)
builds a ClinVar x heart dataset from `clinvar_heart.yaml`:

```bash
hvantk training-set build --spec clinvar_heart.yaml \
    --input clinvar=clinvar.ht --input dbnsfp=dbnsfp.ht \
    --input insider_variants=insider_variants.ht --input ptm=ptm_sites.ht \
    --input hgnc=hgnc.ht --gene-matrix gene_matrix.ht \
    --gene-list is_chd_gene=chd_genes.txt \
    --output out/clinvar_heart
```

The example's README has the `expression summarize`, `annotate prepare` and
`annotate compose` steps that build `gene_matrix.ht`.

## The spec

The training-set sections sit beside `layer1` in the feature spec; `annotate` ignores
them.

```yaml
label:
  input: clinvar                 # the --input name
  source: "clinvar:variants"
  provenance: [ClinVar]          # what the label derives from (drives L2)
  pathogenic: [Pathogenic, Likely_pathogenic, Pathogenic/Likely_pathogenic]
  benign: [Benign, Likely_benign, Benign/Likely_benign]
  min_review_stars: 1            # default 1
  consequence: missense_variant  # default
  snv_only: true                 # default
variant:
  - input: dbnsfp
    source: "dbnsfp:variants"    # its plugin's `scores: trained_on` drive the tiers
    join: allele                 # same locus and alleles
    gene_field: Ensembl_geneid   # breaks multi-gene ties (at most one source)
    columns:
      REVEL_rankscore: {field: REVEL_rankscore}
      CADD_raw_rankscore: {field: CADD_raw_rankscore, tier: L2}
      gnomad_af_popmax:
        fields: [gnomAD.gnomAD_exomes_POPMAX_AF, gnomAD.gnomAD_genomes_POPMAX_AF]
        cast: float64            # parse a string field
        reduce: max
        tier: L4
  - input: insider_variants
    source: "insider:variants"
    join: interval               # locus inside an interval: one boolean column
    columns: {ppi_site: {tier: L1}}
  - input: ptm
    source: "uniprot-ptm:sites"
    join: codon                  # locus inside a site's codon: one boolean column
    columns: {ptm_site: {tier: L1}}
gene:
  hgnc_input: hgnc
  spine_columns: {cds_length: {tier: L1}}
gene_lists:
  is_chd_gene: {description: "..."}
```

Each `layer1` axis also takes a `tier:`; an axis without one is left out of the training
set. The full schema is `hvantk/resources/schemas/feature_spec.schema.json`.

## How rows, labels and genes are derived

- **Rows**: biallelic SNVs (with `snv_only`) whose ClinVar `MC` lists `consequence`.
- **Label**: the primary classification is the text of `CLNSIG` before the first `|`,
  so `Pathogenic|risk_factor` is pathogenic. It gives 1 when in `pathogenic`, 0 when in
  `benign`; any other row is dropped.
- **Review stars**: from `CLNREVSTAT` re-joined with `,` (Hail splits it at every
  comma). `practice_guideline` 4, `reviewed_by_expert_panel` 3,
  `criteria_provided,_multiple_submitters,_no_conflicts` 2, other `criteria_provided`
  statuses 1, anything else 0. Rows below `min_review_stars` are dropped.
- **Gene**: ClinVar's `GENEINFO` genes, in order, are mapped Entrez -> HGNC -> Ensembl.
  A single gene is used as is. When several are listed, the first one the `gene_field`
  source also assigns to the variant is used, falling back to the first listed gene.
- **Gene columns** are joined from the matrix by `gene_id`. A gene that is not on the
  spine has every gene column missing, its `{axis}_present` flags included. Missing values
  are never filled.

## The tier rule

Nested training arms are `L1`, `<= L2`, `<= L3` and all, so a higher tier keeps a column
out of more arms. Each column's tier is the **higher** of two:

- **Derived from its declared `trained_on`** (the plugin's `scores:` block): **L2** if it
  overlaps the label's `provenance`, **L3** if it contains `population_af` (the higher one
  when both apply), otherwise **L1**, including for an empty `trained_on`. Curated disease
  databases (ClinVar, ClinGen, GenCC, HGMD, OMIM) count as one source here, the same
  equivalence `hvantk rerank` uses.
- **Declared in the spec** with `tier:`: **L4** for a population allele frequency or a
  score that takes one as input (`trained_on` describes training data, not inputs),
  **L3** for gene constraint, **L2** for e.g. predictors trained on simulated vs observed
  alleles, **L1** for expression, structure and protein context.

A column with neither is **left out** and listed under `excluded` in the manifest.

A declaration can raise a tier but never lower it, so it cannot hide circularity: a column
trained on ClinVar is at least L2 whatever the spec says. A column built from several
scores (`fields` + `reduce`) takes the union of their `trained_on`, and the score columns
of a `layer1` axis that aggregates per-variant scores (`key: variant`) take their source
scores' `trained_on`; one undeclared input makes the column undeclared. `{axis}_present`
flags are L1.

## Options reference

| Option | Required | Description |
| --- | --- | --- |
| `--spec` | yes | Feature-spec YAML |
| `--input` | yes | Repeatable `NAME=PATH`, exactly the names the spec declares |
| `--gene-matrix` | yes | `annotate compose` output from the same spec |
| `--gene-matrix-manifest` | no | Compose manifest (default `<gene-matrix>.manifest.json`) |
| `--gene-list` | if the spec has `gene_lists` | Repeatable `COLUMN=FILE` |
| `--legacy-names` | no | TSV with header `column`, `legacy_name` |
| `--output` | yes | Output prefix |

Name mismatches with the spec, an empty gene list, a missing field path (in a spec column,
the `gene_field`, a codon source's `chrom`/`codon_start`/`codon_end` or HGNC's
`entrez_id`/`ensembl_gene_id`), or a gene matrix that lacks a declared column all stop the
command with a usage error before Hail runs any job.

## Outputs

| File | Contents |
| --- | --- |
| `<prefix>.ht` | The dataset, keyed by `locus`, `alleles`, with a provenance sidecar |
| `<prefix>.tsv.gz` | The same rows as a [deterministic gzipped TSV](../guide/usage.md#output-extension-by-dataset-type) |
| `<prefix>.dictionary.tsv` | Per column: `column`, `legacy_name`, `training_feature`, `tier`, `level` (`variant`/`protein`/`gene`), `source`, `source_version`, `description`, `non_null_rate`, `true_rate` (booleans) |
| `<prefix>.manifest.json` | Counts at each filter, label balance overall and per gene list, review stars, rows per input, how each gene was chosen, tier counts, excluded columns, gene-list mapping (symbols in, symbols mapped, genes, every unmapped symbol), inputs and their provenance |

## Caveats

- Only `clinvar:variants` labels are supported.
- The gene columns' `source_version` is empty: `annotate compose` records no provenance.
- A `--gene-list` symbol that does not map onto the spine is listed in the manifest, not
  an error.

## See also

- [Gene annotation matrix](annotate.md): the gene layer this command joins.
- [Re-rank](rerank.md): the same `trained_on` declarations drive its clean and all arms.
