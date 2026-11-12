# The spec, the tiers, and choosing training columns

One feature-spec YAML declares the whole dataset, and the same file drives
`annotate prepare`/`compose` for the gene matrix. The command's page,
`docs_site/tools/training-set.md`, is the reference for every key and for the full tier
rule; this page covers the decisions a run needs.

## Adapting the example spec

Start from a copy of `examples/training_sets/clinvar_heart.yaml` and change only what the
task asks for.

- **Another gene list.** Rename the `gene_lists` entry (its key is the output column,
  e.g. `is_kidney_gene`) and pass `--gene-list is_kidney_gene=FILE`. The file is gene
  symbols separated by whitespace. A gene list flags rows; it does not filter them: the
  dataset still holds every ClinVar gene, and a model for the listed genes trains on the
  rows where the flag is true. The flag is not a training feature (`training_feature =
  no` in the dictionary), and the manifest reports the label balance inside the list and
  every symbol that did not map.
- **Fewer or other gene axes.** Delete or add `layer1` entries, then rebuild the gene
  matrix from the edited spec: `training-set build` checks the matrix against the spec
  and stops with `axis '<axis>' is not in the gene matrix's manifest; was the matrix
  composed from this spec?` or `gene matrix lacks column(s): ...` when they differ. A
  matrix composed from the old spec is never reused.
- **Other variant columns.** `field: <path>` takes one field (dotted paths reach into
  structs, e.g. `gnomAD.gnomAD_exomes_AF`); `fields: [...]` with `reduce: max` takes the
  maximum of several; `cast: float64` parses a field stored as text. A path that does not
  exist stops the build before any job with `field '<path>': '<segment>' is not a field
  of ...`. An `interval` or `codon` source gives exactly one boolean column.
- **The label.** `pathogenic` and `benign` list the primary-CLNSIG values for 1 and 0;
  `min_review_stars` (default 1), `consequence` (default `missense_variant`) and
  `snv_only` (default true) filter the rows. Only `clinvar:variants` can be the label
  source. Changing any of these changes the rows, so record the spec's SHA-256 (the
  manifest does: `spec.sha256`).
- **Old column names.** `--legacy-names` takes a TSV with the header `column` and
  `legacy_name`; it only fills the dictionary's `legacy_name` column.

## Tiers: what to declare

Every training feature gets a tier, L1 to L4: the **higher** of the tier derived from
the plugin's `trained_on` declarations and the tier the spec declares with `tier:`. A
declaration can raise a tier, never lower it. A column with neither is left out of the
dataset and listed under `excluded` in the manifest.

Declare a `tier:` for:

| Column | Tier |
| --- | --- |
| a population allele frequency, or a score that takes one as an input | `L4` |
| gene constraint (depletion of variation in a population) | `L3` |
| a predictor trained on simulated versus observed alleles | `L2` |
| expression, structure, protein context, conservation | `L1` |

`trained_on` describes training data, not inputs, so the plugin cannot know that a score
uses allele frequency as a feature: that is what the spec's `L4` declaration is for.
Every `layer1` axis needs a `tier:`, or its columns are excluded. A variant source that
is not a registered dataset has no `trained_on` at all: the build logs
`source '<source>' is not a registered dataset; its columns need a tier in the spec`,
and each of its columns needs a `tier:`.

After the build, read the manifest's `tiers` and `excluded`. Anything in `excluded` was
dropped on purpose by the rule above; give it a tier in the spec, or leave it out, and
rebuild.

## Choosing training columns

The tiers nest into four training arms. Each keeps a column only when its tier is at or
below the arm's limit:

| Arm | Keeps | Leakage it avoids |
| --- | --- | --- |
| `L1` | `L1` | anything trained on the label's sources, on allele frequency or using it |
| `<= L2` | `L1`, `L2` | allele frequency and constraint |
| `<= L3` | `L1` to `L3` | allele frequency as a value or an input |
| all | `L1` to `L4` | nothing: a reference arm, not a clean model |

Pick the columns from the dictionary, never by name patterns:

```python
import pandas as pd

d = pd.read_csv("OUT/clinvar_heart.dictionary.tsv", sep="\t", dtype=str).fillna("")
arm = {"L1"}  # or {"L1", "L2"}, {"L1", "L2", "L3"}, {"L1", "L2", "L3", "L4"}
features = d.loc[(d.training_feature == "yes") & d.tier.isin(arm), "column"].tolist()
```

`training_feature = no` marks the label, the identifiers (`locus`, `alleles`, `vid`,
`gene_id`, `gene_symbol`), the ClinVar fields the label was derived from
(`clinvar_clnsig`, `clinvar_review_status`, `clinvar_review_stars`, `consequence`) and
the gene-list flags: never train on them.

Two kinds of leakage no column choice removes:

- **Gene identity.** Gene-level columns repeat across every variant of a gene, and the
  label balance differs a lot between genes, so a random split by variant lets a model
  learn each gene's label rate. Split train and test by `gene_id`.
- **The labels themselves.** ClinVar submitters classify partly from allele frequency
  and from computational and conservation evidence (the ACMG criteria BA1, BS1, PM2,
  PP3, BP4), so even an `L1` predictor agrees with the labels partly by construction.
  Comparing the same model, on the same gene-grouped folds, across the four arms
  measures how much the higher tiers add.
