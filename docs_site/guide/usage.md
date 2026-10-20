# Building datasets

This guide covers the data layer: discovering what hvantk can build, building a dataset
with `hvantk reprocess`, checking a source for drift, and browsing the dataset catalog.
The analysis tools that consume built datasets (annotation matrices, cohorts,
expression, gene sets and the rest) each have a page under Tools; [See also](#see-also)
links them.

If you haven't installed hvantk yet, see [Installation](../getting-started/installation.md).

## Discover datasets

Every in-tree data source ships as a plugin under `hvantk/skills/<provider>/`, declared
by a `plugin.yaml` manifest: its datasets, lifecycle (download / parse / build), drift
probe, and dataset contract. Every analysis command ships the same way, as a
`tool.yaml` manifest. Two command groups inspect these registries.

### `hvantk plugins`

```bash
# List every loaded provider and how many datasets each declares
hvantk plugins list

# Full detail for one provider: its datasets, skill paths, drift fingerprints
hvantk plugins describe clinvar

# Name anything that failed to load, and why
hvantk plugins errors

# Validate a plugin.yaml against the schema, offline
hvantk plugins validate hvantk/skills/clinvar/plugin.yaml
```

`list` and `describe` are informational and exit 0 even when something failed to load
(`describe` exits 1 only if the named provider itself is unknown or failed to load as a
whole); `errors` is the one whose exit code reports load failures (1 if it lists
anything, 0 if clean); `validate` exits 1 on a schema violation.

### `hvantk tools`

```bash
hvantk tools list
hvantk tools describe annotate
hvantk tools errors
hvantk tools validate hvantk/tools/annotation/annotate_cli.tool.yaml
```

Same shape as `plugins`: `list` and `describe` exit 0; `errors` exits 1 if it lists
anything, 0 if clean; `validate` exits 1 on a schema violation.

## Build a dataset

`hvantk reprocess <provider>:<dataset>` chains a plugin's lifecycle stages (download,
parse, build), then runs the dataset's drift probe:

```bash
hvantk reprocess clinvar:variants \
  --raw-dir data/clinvar/ \
  --output out/clinvar.ht
```

A dataset declared `acquisition.mode: byo` (see [Data Sources](data-sources.md)) treats
`--skip-download` as implicit and instead requires `--raw-dir` to already hold files --
unless `--skip-parse` is combined with an explicit `--intermediate`, in which case the
raw directory is never read and does not need to be populated.

### Options

| Option | Required | Description |
| --- | --- | --- |
| `DATASET` (argument) | yes | `<provider>:<dataset>`, e.g. `clinvar:variants`. |
| `--raw-dir DIRECTORY` | yes | Where to write (or, with `--skip-download`, read) raw files. |
| `--output PATH` | yes | Final builder output path. Extension must match the dataset's native format (see below); a mismatch fails before any stage runs. |
| `--intermediate PATH` | no | Path for the parsed intermediate. If the plugin declares `lifecycle.parse` and you omit this, it defaults to `<raw-dir>/<provider>_<dataset>.intermediate` and the default is logged. |
| `--skip-download` | no | Skip the download stage; read from `--raw-dir` instead. |
| `--skip-parse` | no | Skip the parse stage; the build reads `--intermediate` if given, else `--raw-dir`. |
| `--skip-build` | no | Skip the build stage. Nothing is written to `--output`. |
| `--plugin-arg KEY=VALUE` | no, repeatable | Forwarded to the plugin's download, parse, and build functions (e.g. `--plugin-arg cancer_type=brca` for `cptac:phospho`, `--plugin-arg reference_genome=GRCh38` for `clinvar`). |
| `--check-drift` / `--no-check-drift` | no | Run the post-build drift probe (default: on) or skip it. |
| `-v`, `--verbose` | no | Surface the builder's own INFO logs (e.g. per-chunk progress) on stderr. |
| `-q`, `--quiet` | no | Suppress the per-stage progress lines; errors still print. |

`--plugin-arg` values are coerced before they reach the plugin: `true`/`false`
(case-insensitive) become `bool`; an integer string becomes `int`; a decimal or
scientific-notation string becomes `float`; a value containing a comma becomes a list
of its comma-separated parts, each coerced the same way; anything else stays `str`.

### Output extension by dataset type

`reprocess` writes each dataset in its own native format -- there is no cross-format
conversion -- so `--output`'s extension must match the dataset's declared
`artifact_type` and `backend`:

| Dataset type | Backend | Extension | Example |
| --- | --- | --- | --- |
| `AnnotationTable` | `hail` | `.ht` | `clinvar:variants` |
| `AnnotationTable` | `pandas` | `.parquet` | `peptideatlas:phospho` |
| `ExpressionMatrix` | `anndata` | `.h5ad` | `cptac:phospho` |
| `VariantMatrix` | `hail` | `.mt` | `onek-genomes:variants` |
| `GeneSet` | -- | `.geneset.json` | none built by an in-tree plugin yet |

### Existing output and provenance

`reprocess` always overwrites an existing `--output`: there is no `--overwrite` flag and
no confirmation prompt. (A few individual tool subcommands -- `hvantk genesets` and
`hvantk expression` -- do refuse an existing output unless you pass their own
`--overwrite`; that is a property of those commands, not of `reprocess`.)

Every build writes a sidecar `<output>.provenance.json` next to the output, recording:

```json
{
  "plugin": "clinvar",
  "dataset": "clinvar:variants",
  "plugin_version": "0.3.0",
  "source_fingerprint": "sha256:...",
  "schema_id": "clinvar-variants-v3",
  "build_timestamp": "2026-10-20T12:00:00+00:00",
  "builder_commit": null,
  "parents": [],
  "build_parameters": {}
}
```

`source_fingerprint` is the drift probe's fingerprint at build time. `parents` is
filled only for a dataset derived from other hvantk datasets. `build_parameters` records
the options a builder chose to stamp, and is `null` only in a sidecar written before the
field existed. `builder_commit` is always `null` for `reprocess` builds today. Loading the
dataset later (`hvantk.core.io.load`) re-attaches this `Provenance`.

### Drift check at the end of a build

With `--check-drift` (the default), `reprocess` runs the dataset's drift probe once the
build finishes and prints `drift: <status>` -- this never changes `reprocess`'s own exit
code, even for `drifted`. Use `--no-check-drift` to skip the probe entirely (e.g.
offline), or run `hvantk drift <dataset>` separately when you need drift status to
affect a script's exit code.

### Downloading without building

`hvantk download <source>` runs one provider's downloader on its own, without the
parse/build/drift stages `reprocess` chains it into. Not every dataset has one -- see
[Data Sources](data-sources.md) for which do, and for the ones that are manual-download
only.

## Check a source for drift

`hvantk drift <provider>:<dataset>` (or `--all`, optionally filtered with `--domain`)
compares the dataset's live drift-probe output against the fingerprint file committed
in the plugin's own tree, as named by its manifest's `tests.drift_fingerprint` -- for
example `tests/drift_fingerprint.json` for a single-dataset plugin, or, for a
multi-dataset one, either a per-dataset subfolder (`phospho/tests/drift_fingerprint.json`)
or a suffixed filename (`tests/drift_fingerprint_variants.json`):

```bash
hvantk drift clinvar:variants
hvantk drift --all --domain genomics --json
```

| Status | Meaning |
| --- | --- |
| `clean` | The live probe output matches the committed fingerprint. |
| `drifted` | The probe ran, but the source has changed since the fingerprint was captured; printed with a diff. |
| `probe_failed` | The probe raised, timed out, or returned something invalid (e.g. unreachable network). |
| `stub` | Documentation-only source with no live probe (not an error; exits clean, with a stderr warning). |

Exit code is the worst outcome across the datasets checked:

| Exit code | Meaning |
| --- | --- |
| 0 | Clean (and, with `--all`, no load errors). |
| 1 | At least one dataset drifted. |
| 2 | A probe failed, a dataset's plugin failed to load, or (with `--all`) nothing was registered to check. |
| 3 | The dataset name is unknown. |
| 4 | An invalid combination of flags. |

`--regenerate`, `--accepted-in`, `--no-ledger-row`, `--ledger`, and `--mark-rebuilt`
rewrite the committed fingerprint and the rebuild ledger; they are for maintainers
promoting a drift result, not for routine checks -- see
[Contributing](../contributing.md).

## Browse the catalog

`hvantk catalog` reads the dataset metadata every plugin declares in its own
`catalog/datasets.json` (not every plugin has one):

```bash
hvantk catalog list --omics-type genomics --limit 10
hvantk catalog search "clinical significance" --omics-type genomics
hvantk catalog show ClinVar_latest --format json
hvantk catalog stats
```

| Subcommand | Key options |
| --- | --- |
| `list` | `--omics-type {transcriptomics,proteomics,genomics,epigenomics}`, `--data-source TEXT` (substring, case-insensitive), `--organism TEXT` (substring), `--limit INTEGER` (default 50) |
| `search QUERY` | Matches title, description, and accession. `--omics-type`, `--limit` (default 50) |
| `show ACCESSION` | `--format {yaml,json}` (default `yaml`) |
| `stats` | No options; aggregate counts by omics type, organism, and data source |

## Utilities

```bash
# Convert standard gzip to BGZF for parallel Hail import
hvantk utils convert-bgz data/scores.tsv.gz --threads 4

# Verify a BGZF file's block integrity
hvantk utils validate-bgzf data/clinvar.vcf.gz --quick

# Smoke-test the install: Hail init, a tiny synthetic MatrixTable, describe()
hvantk utils check-install
```

| Command | Options | Notes |
| --- | --- | --- |
| `convert-bgz INPUT_FILE` | `-o, --output PATH` (default: replaces `.gz` with `.bgz`); `-t, --threads INTEGER` (default 4) | No-ops if the input is already BGZF; exits 1 if the input isn't gzip at all. |
| `validate-bgzf PATH` | `--quick` (header-only; skips CRC32/ISIZE) | `PATH` may be a directory, in which case every `*.vcf.gz` inside it is checked. Exits 1 if any file fails. |
| `check-install` | none | The one `utils` command that starts Hail/Spark: it initializes Hail, builds a `balding_nichols_model` MatrixTable (3 populations, 50 samples, 200 variants), and prints its schema. |

Hail reads standard gzip (`.gz`) single-threaded; BGZF (`.bgz`) parallelizes, which
matters for large raw files before `hvantk reprocess` imports them.

## Tips and troubleshooting

- A plugin's on-disk directory name and its manifest `name:` sometimes differ by `_`
  versus `-` (e.g. the `gwas_catalog` directory declares `name: gwas-catalog`). If
  `hvantk plugins describe <provider>` or `hvantk drift <provider>:<dataset>` says
  "unknown provider" or "unknown dataset", run `hvantk plugins errors` first -- a
  provider that *failed to load* can look identical to a typo.
- A build that fails immediately with a message about `--output`'s extension is
  `reprocess`'s own fast-fail check, before any stage runs -- see the extension table
  above, not a Hail error.
- A `stub` result from `hvantk drift` means the source has no live, programmatic probe
  (documentation-only acquisition); it is expected, not a failure.

## See also

- [Data Sources](data-sources.md) -- available datasets and how to acquire them
- [Gene annotation matrix](../tools/annotate.md) -- build a gene x feature matrix from built datasets
- [External cohorts](../tools/cohort.md) -- bring a cohort's own results alongside it
- [Expression analysis](../tools/expression.md) -- summarize and extract markers from expression AnnData
- [Gene sets](../tools/genesets.md) -- extract and prepare `GeneSetCollection` files
- [Architecture](../architecture.md) -- system design and extension points
- [Contributing](../contributing.md) -- development workflow and contribution guidelines
