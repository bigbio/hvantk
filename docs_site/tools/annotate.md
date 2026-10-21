# Gene annotation matrix

`hvantk annotate` turns already-built datasets into one gene x feature matrix ("Layer
1"). It is three commands run in order -- `spine`, `prepare`, `compose` -- because each
stage is independently re-runnable and the expensive middle stage (`prepare`) is
per-source, so it is the one worth parallelizing across a cluster.

## When to use it

Use it once you have a canonical gene table and the per-source datasets you want to
join onto it (built separately with [`hvantk reprocess`](../guide/usage.md#build-a-dataset)),
and you want them combined into a single `gene_id`-keyed table, for example as the
feature axes of [`hvantk rerank`](rerank.md) or any other analysis that expects one row
per gene.

## Inputs

- `ensembl-gene:structure` -- the canonical gene table (`--gene-table`), keyed on `gene_id`.
- `hgnc:lookup` -- the HGNC gene catalog (`--hgnc`), used both to build the spine and,
  in `prepare`, to resolve any source key that isn't already `gene_id`.
- A feature-spec YAML (`--spec`) naming each source, its join key, and the columns to
  keep -- see [Feature spec reference](#feature-spec-reference) below.
- Per axis, the built source itself: an `AnnotationTable` (`.ht`) for a `gene_id` /
  `hgnc_id` / `symbol` / `uniprot_id` / `variant` entry, or an `ExpressionMatrix`
  (`.h5ad`) for a `matrix` entry.

## Worked example

### 1. Build the spine

```bash
hvantk annotate spine \
  --gene-table out/ensembl_structure.ht \
  --hgnc out/hgnc_lookup.ht \
  --output out/spine.ht
```

`--biotype` defaults to `protein_coding`; pass `--biotype all` to keep every biotype.
The `hgnc_id` join is a left join on purpose: a protein-coding Ensembl gene with no
HGNC record stays in the spine with a missing `hgnc_id` rather than being dropped.

### 2. Declare the sources

```yaml
name: chd-features
layer1:
  - {axis: constraint, source: gnomad-metrics:metrics, key: gene_id, columns: [mis_z]}
  - {axis: gevir,      source: gevir:metrics,          key: gene_id, columns: [gevir_pct]}
```

`layer1:` is a literal key, not a label you choose -- it names the gene-level evidence
axes. Each `axis` label must be unique: it addresses the entry in the next two stages
(`--axis constraint`, `--prepared constraint=...`) and names the `{axis}_present` column
`compose` adds.

### 3. Prepare each axis, then compose

```bash
hvantk annotate prepare --spec features.yaml --axis constraint \
  --input out/gnomad_metrics.ht --spine out/spine.ht --output out/constraint.ht

hvantk annotate prepare --spec features.yaml --axis gevir \
  --input out/gevir_metrics.ht --spine out/spine.ht --output out/gevir.ht

hvantk annotate compose --spec features.yaml --spine out/spine.ht \
  --prepared constraint=out/constraint.ht \
  --prepared gevir=out/gevir.ht \
  --output out/features.ht
```

`compose` requires exactly one `--prepared` per axis the spec declares -- missing one
is a `click.UsageError` before any Hail call; naming one the spec doesn't declare is
the same. The manifest JSON lands next to the output
(`out/features.ht.manifest.json` unless `--manifest` says otherwise).

## Options reference

### `annotate spine`

| Option | Required | Description |
| --- | --- | --- |
| `--gene-table TEXT` | yes | Path to the `ensembl-gene:structure` table. |
| `--hgnc TEXT` | yes | Path to the `hgnc:lookup` table. |
| `--output TEXT` | yes | Output path for the spine table. |
| `--biotype TEXT` | no | Gene biotype to keep. Default `protein_coding`; pass `all` to keep every biotype. |

### `annotate prepare`

| Option | Required | Description |
| --- | --- | --- |
| `--spec TEXT` | yes | Path to the feature-spec YAML. |
| `--axis TEXT` | yes | Which `layer1` entry (by its `axis` label) to prepare. |
| `--input TEXT` | yes | The built source: an `AnnotationTable` (`.ht`), or an `AnnData` (`.h5ad`) for a `matrix` entry. |
| `--spine TEXT` | yes | Path to the gene spine table (`.ht`). |
| `--output TEXT` | yes | Output path for the prepared table (`.ht`). |
| `--hgnc TEXT` | only when the rule below requires it | Path to the `hgnc:lookup` table. |

### `annotate compose`

| Option | Required | Description |
| --- | --- | --- |
| `--spec TEXT` | yes | Path to the feature-spec YAML. |
| `--spine TEXT` | yes | Path to the gene spine table (`.ht`). |
| `--prepared TEXT` | yes, repeatable | `axis=path`, one per spec `layer1` axis. |
| `--output TEXT` | yes | Output path for the composed matrix (`.ht`). |
| `--manifest TEXT` | no | Output path for the manifest JSON. Default `<output>.manifest.json`. |

None of the three commands has an `--overwrite` flag: each always overwrites an
existing `--output` (`ht.write(..., overwrite=True)`), with no confirmation prompt.

## The `--hgnc` rule

`prepare --hgnc` is required exactly when the axis can't resolve onto `gene_id` on its
own:

- `key: hgnc_id`, `key: symbol`, or `key: uniprot_id` -- always needs `--hgnc`.
- `key: variant` with an `aggregate.to` other than `gene_id` -- needs `--hgnc` to carry
  the aggregated id the rest of the way onto `gene_id`.
- `key: variant` aggregating straight to `gene_id`, or a plain `key: gene_id` entry --
  does not need it.
- A `matrix` block (expression-atlas-shaped sources) is always `key: symbol` (the schema
  enforces this when loading the spec), so it always needs `--hgnc`.

`hvantk annotate prepare --help` states the rule in its own words; the logic lives in
`needs_hgnc` in `hvantk/tools/annotation/annotate_cli.py`.

## Feature spec reference

Validated against
[`feature_spec.schema.json`](https://github.com/bigbio/hvantk/blob/main/hvantk/resources/schemas/feature_spec.schema.json).
Two keys are functional: `name` and `layer1` (a non-empty list). Each `layer1` entry:

| Field | Required | Notes |
| --- | --- | --- |
| `axis` | yes | Unique label; see above. |
| `source` | yes | A `provider:dataset` string, used in log lines and error messages; `--input` on `prepare` is what actually locates the file, not this field. |
| `key` | yes | One of `gene_id`, `hgnc_id`, `symbol`, `uniprot_id`, `variant`. |
| `columns` | yes | Output column names, non-empty. For a `variant` entry these must exactly match what `aggregate` produces (see below); for a `matrix` entry they are checked against, but not limited to, what the reducer emits. |
| `min_mapping_rate` | no | Default `0.9`. `prepare` raises if the axis maps below this fraction of its input rows onto the spine. |
| `collapse` | no | One of `max`, `min`, `mean`, `sum`. Required whenever a source maps many-to-one onto a `gene_id` (e.g. several UniProt accessions for one gene); omitting it is an error, not a silent drop. |
| `origin`, `ablate_separately` | no | Accepted by the schema and stored on the parsed entry; nothing in `prepare` or `compose` currently reads either one. |
| `aggregate` | required iff `key: variant` | See below. |
| `matrix` | only for an expression-atlas-shaped source | Requires `key: symbol`. `group_axis`, `atlas` required; `tissue_tag`, `stats` (`mean`/`fraction_expressed`/`frac`), `drop_groups`, `specificity` optional. |

A `variant` entry's `aggregate` block (`by`, `to`, `scores` required): groups rows by
`by` (a `;`-delimited gene-string column, dbNSFP's shape -- a single value still works),
reduces each named score to a scalar per variant (`reduce: max`, the default, collapses
a per-transcript dict; `reduce: identity` passes an already-scalar column through), then
aggregates to `to` (one of `gene_id`, `hgnc_id`, `symbol`, `uniprot_id`) with the
declared `stats` (`mean`, `max`, `min`, `count`, or `frac_gt_<threshold>`). Each
`{score}_{stat}` pair becomes one output column (e.g. `revel_mean`), plus a row-count
column named by `count_name` (default `n_possible_missense`). `filter` accepts only the
named predicate `"missense"` -- this is a closed vocabulary, not a free-form expression.

The schema also accepts top-level `spine` and `exclude` keys, but `load_spec` (in
`hvantk/algorithms/annotation/spec.py`) never reads either -- they currently have no
effect on the build.

## Outputs

- **Spine** (`annotate spine`): keyed on `gene_id`; columns `gene_name`, `chromosome`,
  `gene_start`, `gene_end`, `gene_biotype`, `cds_length`, `n_coding_exons`,
  `n_transcripts`, `mane_select`, `hgnc_id`, `gene_group`.
- **Prepared** (`annotate prepare`): keyed on `gene_id`; exactly the axis's declared
  `columns`, restricted to genes on the spine.
- **Composed** (`annotate compose`): the spine's columns, every axis's declared columns,
  and one `{axis}_present` boolean per axis (true iff that gene had a row in the axis's
  prepared table -- a value can be missing *and* `{axis}_present` false, distinguishing
  "absent from source" from any other reason a value is null).
- **Manifest** (`<output>.manifest.json`): `{"n_genes": N, "axes": {"<axis>": {"<column>":
  {"non_null_rate": ..., "positive_rate": ...}}}}` -- non-null rate is over all spine
  genes, positive rate is the fraction of non-null values that are `> 0`.

## Python API

```python
import hail as hl
from hvantk.core.utils.hail_context import init_hail
from hvantk.algorithms.annotation.spec import load_spec
from hvantk.algorithms.annotation.spine import build_spine, spine_mapping_rate
from hvantk.algorithms.annotation.prepare import prepare_source
from hvantk.algorithms.annotation.mapping import enforce_rate
from hvantk.algorithms.annotation.compose import compose

init_hail()

spine = build_spine(
    hl.read_table("out/ensembl_structure.ht"),
    hl.read_table("out/hgnc_lookup.ht"),
    biotype="protein_coding",
)
print(spine_mapping_rate(spine))  # fraction of spine genes with an hgnc_id

spec = load_spec("features.yaml")  # the two-axis spec above
spine_ids = spine.gene_id.collect()

prepared_by_axis = {}
for axis, source_ht in [
    ("constraint", hl.read_table("out/gnomad_metrics.ht")),
    ("gevir", hl.read_table("out/gevir_metrics.ht")),
]:
    entry = spec.entry(axis)
    prepared, report = prepare_source(source_ht, spine_ids, entry)
    enforce_rate(report, entry.min_mapping_rate)  # raises MappingRateError below the threshold
    prepared_by_axis[axis] = prepared

composed, manifest = compose(spine, prepared_by_axis, spec)  # compose needs every axis
```

This mirrors the CLI one-to-one: the CLI is a thin wrapper over these same functions.
`prepare_variant_source` and `prepare_matrix_source` (same module) take the same
`(source, spine_gene_ids, entry, *, hgnc=None)` shape for `key: variant` and `matrix`
entries respectively.

## Caveats

- `compose` raises if two axes declare the same output column name -- checked once in
  pure Python against the spec, and again against the prepared tables' actual columns
  (a `matrix` axis's per-group vector columns aren't known until then).
- Two `matrix` axes may not share an `atlas` label: their vector columns are named
  `{atlas}_{group}` with no axis component, so a shared atlas would collide. This is
  checked at spec-load time, before any Hail call.
- Scheduling is deliberately external: all three commands run in-process, and hvantk
  never imports a scheduler -- an sbatch wrapper submits them.

## See also

- [Building datasets](../guide/usage.md) -- build the per-source datasets this consumes
- [External cohorts](cohort.md) -- attach a cohort's own results onto this matrix
- [Re-rank](rerank.md) -- re-rank genes using this matrix plus a cohort and labels
