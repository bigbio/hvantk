# The cohort manifest

A cohort manifest (`cohort.yaml`) declares how one external or self-generated
case/control cohort is presented to hvantk: its gene key, its cohort-derived prior
statistic (and its direction), and, optionally, where its labels live and which of its
columns form which feature axis. Parsed and validated by
`hvantk.algorithms.cohort.spec.load_cohort` against
`hvantk/resources/schemas/cohort_manifest.schema.json` (`additionalProperties: false`
at every level -- an unknown key is a schema error, not a silent no-op). `hvantk cohort
validate` and `hvantk cohort attach` both load it the same way, so they can never
disagree about its shape.

## Fields

| field | required | default | meaning |
| --- | --- | --- | --- |
| `name` | yes | -- | a label; shows up in printed messages and the attach report's `cohort` field. |
| `key` | yes | -- | the identifier SPACE the gene key lives in: `gene_id` (Ensembl), `hgnc_id`, or `symbol`. Decides which mapper `attach` dispatches to. |
| `table` | yes | -- | path to the cohort's gene-level table (TSV or CSV; `.gz`/`.bgz` ok; the delimiter is sniffed from the header line, not the extension). See "Where `table` is resolved from" below. |
| `prior.column` | yes | -- | the column holding the statistic that defines the ranking. |
| `prior.direction` | yes | -- | `lower_is_better` or `higher_is_better`. Never inferred -- a p-value and a Z-score or Bayes factor point opposite ways, and guessing silently inverts the ranking. |
| `key_column` | no | `key`'s own value | the COLUMN NAME in the table that actually holds the `key`-space identifiers. Real tables routinely differ: a symbol-keyed table commonly names that column `gene`, not `symbol` -- write `key: symbol` / `key_column: gene` for that case. |
| `min_mapping_rate` | no | `0.9` | `attach` raises `MappingRateError` if fewer than this fraction of the table's gene keys map onto the spine. |
| `labels.gene_set` | no (block optional) | -- | path to a `GeneSetCollection` JSON (never a bare column) -- see `hvantk genesets prepare`/`clingen`/`gencc`. |
| `labels.set_name` | no | -- | which named set inside the collection, if it holds more than one. |
| `labels.min_mapping_rate` | no | `0.9` | mapping-rate floor for the label set specifically (labels are always symbol-keyed, independent of the cohort's own `key`). |
| `cohort_axes` | no | `[]` | a list of `{axis, columns}` groups -- feature axes a downstream consumer (`rerank-evaluate`'s architecture audit) reads by axis name. |

`declared_columns()` is `prior.column` followed by every axis's columns, in
declaration order; every one of those columns must appear in the table (`validate`
checks this) and must be unique across the whole manifest -- a column cannot be both
the prior and an axis member, and no two axes can share a column name. This is
enforced at construction time, including for a manifest built directly in Python, not
only one parsed from YAML.

## `key` vs `key_column`

These answer two different questions, and conflating them is the exact bug the pair of
fields exists to prevent:

- `key` is the identifier SPACE -- what `attach` needs to know to pick the right
  mapper (`gene_id`/`hgnc_id`/`symbol`). The enum stays closed to those three
  spine-mappable spaces.
- `key_column` is the COLUMN NAME that actually holds those values in the cohort's own
  table.

When they match (a manifest whose identifier column happens to be named after its own
id-space -- e.g. `key: symbol` over a column literally named `symbol`), nothing extra
is needed. They do NOT match for burden's own `--emit-manifest` output even by
default: `render_cohort_manifest` always sets `key_column: "gene"` regardless of what
`key` is -- a `gene_id`-keyed manifest from burden still has `key: gene_id` /
`key_column: gene` as two different literal strings, because the burden TSV's gene
column is always named `gene` (`hvantk/algorithms/burden/fet.py` renames it that way),
never `gene_id`. When the two differ, `attach` silently renames the column internally
before resolving identifiers; get `key_column` wrong and the failure is
`cohort table ... has no column 'x' (key_column)`, which already names both fields so
you can tell which one is wrong.

## Why the `architecture` axis matters to `rerank-evaluate`

`hvantk cohort burden`'s `--emit-manifest` always emits a `cohort_axes` entry named
`architecture` holding `n_case_var`, `conc`, `driver_af` -- the per-gene reductions of
the burden op's winning route (how many case-carrying variants, how concentrated the
case carriers are on the top variant, and that top variant's control-carrier
frequency). `rerank-evaluate`'s advisory `CaseControlArchitectureAudit` flags a gene
whose signal looks like it is really driven by one recurrent, not-actually-rare
variant.

The real wiring condition is slightly looser than "an axis named `architecture`": the
rerank CLI (`hvantk/tools/rerank/rerank_cli.py`, `has_architecture_columns`,
`hvantk/algorithms/rerank/audit.py`) wires the audit whenever the **union of every
`cohort_axes` entry's columns** -- whatever each axis is individually named -- is a
superset of `n_case_var`, `conc`, `driver_af`; it never looks at axis names, only at
which columns are declared somewhere. In practice keep them together under
`architecture` anyway, because the CLI's own near-miss warning is keyed on that name:
if an axis is literally named `architecture` but is MISSING one of the three columns
(e.g. a hand-edited manifest kept only `n_case_var`/`conc`), and no other axis
supplies it either, the CLI logs a warning naming exactly which column is missing --
but still leaves the audit unwired (`NoAudit`), silently, rather than raising. Splitting
the three columns across two differently-named axes also wires the audit (the union
test does not care), but gets no such warning if something is missing, since nothing
is then named `architecture` to warn about.

## `labels`

A pointer to a `GeneSetCollection` JSON, never a bare column -- label sets are always
symbol-keyed (`hvantk genesets prepare` rejects Ensembl or Entrez ids outright), so
`attach` needs the HGNC catalog (`--hgnc`) to resolve them onto the spine regardless of
the cohort's own `key`. Declare `labels` only once you are about to rerank; a
credibility or veto report has no use for them.

## Two examples

Burden's own output, as `--emit-manifest` writes it (edit only `name`):

```yaml
name: my-cohort
key: gene_id
key_column: gene
table: /abs/path/gene_burden.tsv
prior:
  column: minp
  direction: lower_is_better
cohort_axes:
  - axis: architecture
    columns: [n_case_var, conc, driver_af]
```

An external, symbol-keyed gene table with its own p-value column and no
variant-level architecture information:

```yaml
name: my-gene-table
key: symbol
key_column: gene_symbol
table: /abs/path/my_gene_table.tsv
min_mapping_rate: 0.85
prior:
  column: p_value
  direction: lower_is_better
labels:
  gene_set: /abs/path/my-genes.geneset.json
```

## How `validate` and `attach` use it

`validate` is pure Python and Hail-free: it loads the manifest, reads only the table's
*header* line (sniffing the delimiter, decompressing `.gz`/`.bgz` transparently), and
checks that `key_column` and every declared column are present. It never reads a data
row, so it cannot catch a duplicate key or a blank identifier -- that is `attach`'s
job, at Hail-join time (see the skill's Traps). `attach` additionally requires
`--hgnc` whenever `key` is `hgnc_id`/`symbol` or `labels` is declared, checks for
column-name collisions against the Layer-1 table and against its own reserved output
names (`cohort_tested`, and `label` when labels are declared) before running anything,
rejects a table with more than one row per gene key, and enforces `min_mapping_rate`
(raising `MappingRateError`) **before** the actual index-join onto the Layer-1
matrix's `gene_id` -- a cohort that fails the rate check never gets as far as touching
Layer-1 data at all.

## Where `table` is resolved from

`table` is kept **verbatim** by every writer and every reader: `--emit-manifest`
writes it as exactly the `--output` path you passed to `cohort burden` (see
`hvantk/algorithms/burden/manifest.py`'s `render_cohort_manifest` docstring), and both
`cohort validate` and `cohort attach` open it relative to the **process's current
working directory** at the moment they run (`hvantk/tools/cohort/cohort_cli.py`'s
`_read_header` passes `manifest.table` straight through to a plain file open) --
**never relative to the manifest file's own location**. A manifest and its table that
sit side by side in the same directory will still break if `validate`/`attach` is run
from a different working directory than the one `table` was written relative to (a
container job's working directory is not guaranteed to match where you authored the
manifest by hand). Give `--output`/`table` an **absolute path** and this never comes
up; if you must use a relative one, run every command from the same directory you
wrote the manifest in.
