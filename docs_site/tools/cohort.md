# External cohorts

`hvantk cohort` brings a cohort's own gene-level results alongside the Layer-1
annotation matrix (built by [`hvantk annotate`](annotate.md)), without imputing
or mixing cohorts. A cohort is declared by a manifest, not by flags: a gene key, a
cohort-derived prior statistic (and its direction), and optionally labels and feature
axes.

## When to use it

Use it to fold an external cohort's burden results, or any other gene-level table, into
the same `gene_id`-keyed universe as your annotation matrix -- as a precursor to
[`hvantk rerank`](rerank.md), or to inspect a cohort's own coverage on its own.

## Inputs

- A cohort manifest YAML -- see [Manifest reference](#manifest-reference) below.
- The cohort's own table: TSV or comma-separated (sniffed from the header line; gzip/bgzf
  transparently decompressed), with one row per gene and no duplicate gene keys. It must
  carry the manifest's `key_column` plus every column the manifest declares (the prior
  column and every `cohort_axes` column).
- Optionally, a labels pointer: a `GeneSetCollection` JSON (the format
  [`hvantk genesets`](genesets.md) and its Python API produce).
- For `hvantk cohort burden`: a clean, split bi-allelic case/control `MatrixTable` with
  genotypes (`GT`), a gene column, a variant-route column, and a boolean case/control
  column -- sample and variant QC are upstream of this command, not its job.
- `attach` additionally needs the Layer-1 matrix (`--layer1`, a `.ht` from
  `hvantk annotate compose`) and, when the cohort's `key` is `hgnc_id` or `symbol` or
  the manifest declares labels, the `hgnc:lookup` table (`--hgnc`).

## Worked example

```bash
# Check the manifest against its table (and labels) before anything reads it
hvantk cohort validate --cohort cohort.yaml

# Join the cohort's declared columns onto the Layer-1 matrix
hvantk cohort attach --cohort cohort.yaml --layer1 out/features.ht --output out/attached.ht

# Originate a gene-level prior from a clean case/control MatrixTable instead
hvantk cohort burden \
  --mt data/cc_clean.mt --gene-col gene_id --route-col route --arm-col is_case \
  --key gene_id --carrier-mode het --output out/burden.tsv --emit-manifest out/cohort.yaml
```

`validate` first is the point: it fails on a key column that isn't present, or a
declared axis column missing from the table, rather than letting a silently-empty join
propagate into the matrix.

## Options reference

### `cohort validate`

| Option | Required | Description |
| --- | --- | --- |
| `--cohort TEXT` | yes | Path to the cohort manifest YAML. |

Prints the manifest's key, prior, and axes, loads the labels (if declared), and exits
nonzero with a specific message on the first thing that doesn't match the table.

### `cohort attach`

| Option | Required | Description |
| --- | --- | --- |
| `--cohort TEXT` | yes | Path to the cohort manifest YAML. |
| `--layer1 TEXT` | yes | Path to the composed Layer-1 matrix (`.ht`). |
| `--output TEXT` | yes | Output path for the attached table (`.ht`). |
| `--report TEXT` | no | Output path for the report JSON. Default `<output>.report.json`. |
| `--hgnc TEXT` | only when the cohort `key` is `hgnc_id`/`symbol`, or labels are declared | Path to the `hgnc:lookup` table. |

No `--overwrite` flag: `attach` always overwrites an existing `--output`.

### `cohort burden`

| Option | Required | Description |
| --- | --- | --- |
| `--mt TEXT` | yes | Clean case/control genotype `MatrixTable` (`.mt`). |
| `--gene-col TEXT` | yes | Row field holding the gene key. |
| `--route-col TEXT` | yes | Row field holding the variant-class "route" (`str` or `array<str>`; a variant tagged with several routes contributes to each). |
| `--arm-col TEXT` | yes | Column field: boolean case/control, `True` = case. Must be defined and non-constant for every sample. |
| `--key [gene_id\|hgnc_id\|symbol]` | yes | Identifier space `--gene-col` values live in. |
| `--carrier-mode [het\|hom\|chet\|homs_chet]` | no | Default `het`. See below. |
| `--score-field TEXT` | no | Row field with a per-variant score (e.g. REVEL); reduced into `mean_score_case`. |
| `--mtc [bh\|bonferroni]` | no | Multiple-testing correction on the gene-level min-p; adds `p_adj`. |
| `--output TEXT` | yes | Output gene table (TSV). |
| `--emit-manifest TEXT` | no | Also write a ready-to-use `cohort.yaml` stub pointing at `--output`. |

`--carrier-mode` sets the gene-level "qualifies" rule: `het` -- at least one
heterozygous carrier; `hom` -- at least one homozygous-variant carrier; `chet` -- at
least two heterozygous carriers (a compound-het proxy); `homs_chet` -- either of the
last two.

For each gene and route, `burden` counts the cases and the controls that qualify and
runs a two-sided Fisher's exact test on that 2×2 table (qualifying or not, case or
control). The gene keeps the route with the smallest p-value as `minp`, with that
route's odds ratio. `minp` is not adjusted for the number of routes tried; `--mtc`
corrects across genes only.

## Outputs

**`attach`** -- the Layer-1 columns, the cohort's declared columns (never imputed: a
gene absent from the cohort table keeps missing values), and `cohort_tested` (true iff
the gene had a row in the cohort table), plus `label` when the manifest declares
labels. The report JSON:

```json
{
  "cohort": "my-cohort",
  "n_genes": 20000,
  "n_tested": 180,
  "key": "symbol",
  "key_column": "gene",
  "mapping_rate": 0.97,
  "prior": {"column": "minp", "direction": "lower_is_better", "non_null_rate": 0.009, "positive_rate": 1.0},
  "axes": {"architecture": {"n_case_var": {"non_null_rate": 0.009, "positive_rate": 0.95}, "...": {}}}
}
```

A `labels` block (`gene_set`, `set_name`, `n_in`, `n_mapped`, `mapping_rate`) is added
when the manifest declares labels. `positive_rate` is `null` for a declared column Hail
typed as non-numeric (one dirty cell, e.g. `<0.001` or `Inf`, is enough -- a warning
names the column).

**`burden`** -- a TSV, one row per gene (its winning, min-p route): `gene`, `route`,
`minp`, `odds_ratio`, `n_case_var`, `conc`, `driver_af`, `mean_score_case`,
`frac_case_private`, plus `p_adj` when `--mtc` is set. `conc` is the within-route
case/control concordance fraction; `driver_af` is the max-case-carrier variant's
control-carrier frequency (an upper-bound leak signal, not a population allele
frequency). With `--emit-manifest`, the stub's `cohort_axes: [{axis: architecture,
columns: [n_case_var, conc, driver_af]}]` is exactly what
[`hvantk rerank`](rerank.md)'s `CaseControlArchitectureAudit` looks for.

No `--overwrite` flag on either command: both always overwrite an existing output.

## Manifest reference

Validated against
[`cohort_manifest.schema.json`](https://github.com/bigbio/hvantk/blob/main/hvantk/resources/schemas/cohort_manifest.schema.json):

```yaml
name: my-cohort
key: symbol              # or gene_id / hgnc_id
key_column: gene         # the column in `table` holding that key
table: cohort.tsv
prior:
  column: minp
  direction: lower_is_better
cohort_axes:
  - axis: architecture
    columns: [n_case_var, conc, driver_af]
```

| Field | Required | Notes |
| --- | --- | --- |
| `name` | yes | |
| `key` | yes | One of `gene_id`, `hgnc_id`, `symbol` -- the identifier *space*. |
| `key_column` | no | The column in `table` that actually holds those values. Defaults to `key`; set it explicitly whenever the table's column is named differently (a symbol-keyed table commonly names its column `gene`, not `symbol`). |
| `table` | yes | Path to the cohort's own table. |
| `min_mapping_rate` | no | Default `0.9`. |
| `prior.column` / `prior.direction` | yes | `direction` is `lower_is_better` or `higher_is_better` -- never inferred, since guessing silently inverts the ranking. |
| `labels.gene_set` | only if declaring labels | Path to a `GeneSetCollection` JSON. |
| `labels.set_name` | no | Required only when the collection holds more than one set. |
| `labels.min_mapping_rate` | no | Default `0.9`. |
| `cohort_axes[].axis` / `.columns` | no | Repeatable; each `axis` label must be unique, and every declared column (the prior's and every axis's) must be unique across the whole manifest. |

## Bridging to `hvantk rerank`

`rerank` does not read a cohort manifest's table directly; it takes its own YAML
config that *points* at the manifest and adds the feature axes and labels:

```yaml
cohort: cohort.yaml        # the manifest above
features:
  - {name: constraint, path: constraint.tsv}
labels:
  path: labels.txt
```

`cohort`, `features`, and `labels` are the only top-level keys `rerank` accepts; an
older config with a top-level `prior:` block is rejected with a hint to delete it --
the prior now comes from the cohort manifest's own `prior:` block. The manifest's
`prior.column` is carried through to `rerank`'s output as `prior_stat`. See
[`examples/rerank`](https://github.com/bigbio/hvantk/tree/main/examples/rerank) for a
complete runnable configuration, and [Re-rank](rerank.md) for the rest of its options.

## Python API

```python
from hvantk.algorithms.cohort.spec import load_cohort
from hvantk.algorithms.cohort.checks import (
    read_header, check_key_column_exists, check_declared_columns_exist, detect_delimiter,
)
from hvantk.algorithms.cohort.attach import attach

manifest = load_cohort("cohort.yaml")
header = read_header(manifest.table)
check_key_column_exists(manifest, header)        # raises ValueError, not click.UsageError
check_declared_columns_exist(manifest, header)

import hail as hl
cohort_ht = hl.import_table(manifest.table, delimiter=detect_delimiter(manifest.table),
                             impute=True, missing=["NA", ""], key=manifest.key_column)
attached, report = attach(hl.read_table("out/features.ht"), cohort_ht, manifest)
```

`attach()` itself does **not** run the `checks` functions, and does not reject a blank
or `NA` key cell -- the CLI (`cohort attach`) runs both before calling it. A Python
caller of `attach()` directly should run them too, in this order, or a malformed
manifest or table surfaces later as a raw Hail error instead of a clear message.

For `hvantk cohort burden`'s pipeline:
`hvantk.algorithms.burden.pipeline.run_from_mt(mt, *, gene_col, route_col, arm_col, key,
carrier_mode="het", score_field=None, mtc=None)` returns the same `pandas.DataFrame`
the CLI writes; `hvantk.algorithms.burden.manifest.render_cohort_manifest(*, name, key,
table, key_column="gene", prior_column="minp")` renders the same manifest stub text.

## Caveats

- The cohort table must have exactly one row per gene key; `attach` raises, naming the
  duplicated keys, rather than silently picking one row.
- A manifest column that collides with `attach`'s own reserved output names
  (`cohort_tested`, and `label` when labels are declared) is rejected before any Hail
  call, as is a declared column that already exists in the Layer-1 matrix -- `attach`
  never silently overwrites either side.
- `burden`'s `--mt` must already satisfy the clean-input contract (split bi-allelic
  rows; the given columns present; the arm column defined and not all-one-value for
  every sample); it fails fast with a specific message rather than computing on dirty
  input.

## See also

- [Gene annotation matrix](annotate.md) -- build the Layer-1 matrix a cohort attaches onto
- [Re-rank](rerank.md) -- re-rank genes using this cohort plus features and labels
- [Gene sets](genesets.md) -- build the `GeneSetCollection` JSON a cohort's labels point at
