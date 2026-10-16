# Adding annotation axes onto the gene index

`annotate prepare` and `annotate compose` (`hvantk/tools/annotation/annotate_cli.py`)
turn the gene index plus a feature-spec YAML into one gene x feature matrix -- the
Layer-1 matrix `cohort attach --layer1` can take, in place of the gene index alone.
`hvantk annotate prepare --help` and `hvantk annotate compose --help` list every flag,
but do not fully document the rule for when `--hgnc` is required (its own help text
names only the `hgnc_id`/`symbol` cases, not `uniprot_id` or the `variant`-aggregate
case the code also enforces -- see "`annotate prepare`" below for the real rule). This
page covers the spec format (`hvantk/algorithms/annotation/spec.py::load_spec`,
`hvantk/resources/schemas/feature_spec.schema.json`) and the two commands' own rules,
verified against that code, not against `--help` alone.

## The feature spec

`load_spec` reads exactly two top-level keys into the object the two commands use:
`name` (a label) and `layer1` (a list of source entries; at least one). Each entry:

```yaml
name: my-features
layer1:
  - {axis: constraint, source: gnomad-metrics:metrics, key: gene_id, columns: [mis_z]}
  - {axis: gevir,      source: gevir:metrics,          key: gene_id, columns: [gevir_pct]}
```

Required per entry: `axis` (a unique label -- it is how later stages address the entry,
`--axis constraint`, `--prepared constraint=…`, and it names the `{axis}_present` column
`compose` emits), `source` (a label, informational), `key` (one of `gene_id`, `hgnc_id`,
`symbol`, `uniprot_id`, `variant`) and `columns` (the source columns to keep).

Optional per entry: `min_mapping_rate` (0-1, default `0.9` -- `prepare` raises if the
entry's actual mapping rate onto `gene_id` falls below it); `collapse`
(`max`/`min`/`mean`/`sum` -- required whenever a source key maps many rows onto one
`gene_id`, e.g. `insider:interfaces` or `uniprot-ptm:sites`, otherwise `prepare` raises
rather than silently picking one); `aggregate` (required when, and only when, `key:
variant` -- aggregates per-variant scores up to a gene; see below); `matrix` (for a
`genes x groups` expression or specificity source -- requires `key: symbol`, since the
reduced table it produces is symbol-keyed); `origin` and `ablate_separately` (parsed
onto the entry by `load_spec` -- a free-form string and a bool, respectively -- but as
of this checkout neither is read by `prepare` or `compose` themselves; accepted spec
metadata, not yet a lever either command acts on).

`aggregate`'s own keys: `by` (the source column carrying the gene identifier(s) to
aggregate onto -- dbNSFP's is a `;`-delimited string, exploded to one row per gene),
`to` (the gene-key space the aggregate targets: `gene_id`, `hgnc_id`, `symbol` or
`uniprot_id`), `scores` (a map of output name -> `{column, stats}`, `stats` being
`mean`/`max`/`min`/`count`/`frac_gt_<threshold>`), and optionally `filter` (a named
predicate), `reduce` (`max` or `identity` -- how to collapse a per-transcript dict
column, such as dbNSFP's, to one scalar per variant BEFORE the group-by; default
`max`; use `identity` for a source whose score columns are already plain scalars,
e.g. a PTM site's `n_observations` or a GTEx eQTL pair's `slope`, where `max` would
raise), and `count_name` (the emitted row-count column's name; default
`n_possible_missense`, true only for dbNSFP -- a source counting something else must
set this, or two aggregate axes collide on one column name at compose time).

A duplicate `axis` label, or an entry whose `key: variant` presence disagrees with
whether it carries an `aggregate` block, fails to load before any Hail call runs.

## `annotate prepare`

Run once per axis:

```bash
hvantk annotate prepare --spec features.yaml --axis constraint \
  --input gnomad_metrics.ht --spine gene-index.ht --output constraint.ht
```

(`--spine` is the flag's name; it takes the gene index's path either way.) `--input` is
a built `.ht` for an ordinary entry, or a built AnnData `.h5ad` for an entry carrying a
`matrix` block (a `genes x groups` expression or specificity source) -- `prepare`
dispatches on whether the entry declares `matrix`, not on the file extension. `prepare`
maps the source's key onto the gene index's `gene_id`, keeps only the entry's declared
columns, and enforces `min_mapping_rate`. Pass `--hgnc <hgnc:lookup .ht>` whenever the
entry's `key` is `hgnc_id`, `symbol` or `uniprot_id`, or a `variant` entry whose
`aggregate.to` is anything other than `gene_id` -- omitting it raises `entry '<source>'
needs --hgnc <hgnc:lookup .ht>`. A `gene_id`-keyed entry never needs it.

## `annotate compose`

```bash
hvantk annotate compose --spec features.yaml --spine gene-index.ht \
  --prepared constraint=constraint.ht \
  --prepared gevir=gevir.ht \
  --output features.ht
```

`compose` left-joins every prepared axis onto the gene index's `gene_id` and is strict:
it requires exactly one `--prepared axis=path` per `layer1` axis the spec declares, and
rejects both a missing axis and a `--prepared` naming an axis the spec does not declare
-- a partial matrix is never produced silently. The output carries the gene index's own
columns, every axis's declared columns, and an `{axis}_present` boolean per axis so a
downstream consumer can tell "missing because the gene had no row in that source" apart
from any other reason a value is null. A manifest JSON lands next to the output
(`<output>.manifest.json` unless `--manifest` says otherwise), recording, per axis and
column, the non-null rate over the gene index and the positive rate among non-null
values.

## What `cohort attach --layer1` can take

`cohort attach --layer1 <path>` (`hvantk/tools/cohort/cohort_cli.py`) only needs a
`gene_id`-keyed Hail table with a `gene_id` field to join a cohort's results onto -- it
does not require the composed matrix specifically. The gene index alone, with no axes
added at all, satisfies that; `--layer1` accepts either the gene index or the composed
matrix, and the composed matrix is a superset of the gene index (its own columns plus
every added axis), never a different shape.
