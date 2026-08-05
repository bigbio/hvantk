---
name: hvantk:conventions
description: Conventions every hvantk plugin assumes. Read first.
status: provisional
---

# hvantk Plugin Conventions

These conventions apply to every per-resource plugin under `hvantk/skills/`. Per-resource skills MAY assume everything below without restating it.

## 1. Repository map

- `hvantk/skills/<provider>/` — the plugin folder. Contains `plugin.yaml`, `builder.py`, `cli.py`, `drift_probe.py`, `SKILL.md`, and `tests/`.
- `hvantk/skills/<provider>/<dataset>/` — for providers that ship more than one dataset (e.g., `cptac/expression/`, `cptac/phospho/`). One `plugin.yaml` per provider declares all datasets; each dataset folder owns its builder, drift probe, CLI, and tests.
- `hvantk/skills/<provider>/shared/` — code reused across two or more datasets in the same provider (e.g., the shared CPTAC dataset class).
- `hvantk/skills/_conventions/SKILL.md` — this file. The shared contract.
- `hvantk/core/utils/hail_helpers.py` — shared Hail Table helpers: `create_table_base()` (import → transform → checkpoint → optional TSV export) and `cleanup_temp_file()` (best-effort temp cleanup). Other shared helpers: `contig_recoding()` in `hvantk/core/utils/genome.py`; `parse_gtex_variant_id()` / `strip_ensembl_version()` in `hvantk/core/utils/qtl_helpers.py`. There is no `hvantk/core/builders/table.py` and no `_create_table_base`.
- `hvantk/core/plugin/api.py` — plugin runtime dataclasses (`Provider`, `DatasetSpec`, `DatasetManifest`, `TestPaths`) plus `PluginLoadError`, `PluginNameCollision`, `DriftProbeError`.
- `hvantk/core/plugin/loader.py` — discovery (filesystem + Python entry points), manifest validation, and lazy callable resolution. The module-level `get_registry()` returns a `PluginRegistry`; `registry.get_dataset("<provider>:<dataset>")` yields the executable `DatasetSpec`. There is no `registry.py`, no `TABLE_BUILDERS` / `MATRIX_BUILDERS`, and no `create_table_adapter()` / `create_matrix_adapter()`.
- `hvantk/core/plugin/run_builder.py` — `run_builder_for_spec(...)` is the sole dispatch path: it runs the drift probe, builds a `BuildContext`, invokes the builder, validates the returned artifact's type/`schema_id` against `plugin.yaml`, then saves and returns provenance.
- `hvantk/tools/` — top-level CLI (`hvantk plugins`, `hvantk drift`, `hvantk reprocess`, `hvantk catalog`). The `reprocess` command lives in `hvantk/tools/plugins/reprocess_cli.py`. Per-provider downloader CLI lives in the plugin's own `cli.py` and is wired by `plugin.yaml`'s `cli:` block.
- `hvantk/skills/<provider>/catalog/datasets.json` — per-plugin dataset catalog (URLs, version cadence, license, per-accession metadata). Aggregated by `hvantk.resources.unified_registry.HvantkRegistry` and surfaced via `hvantk catalog {list,show,stats,search}`.

When in doubt, READ existing code under these paths before inferring shape.

## 2. Authoritative spec sources

Each plugin's `catalog/datasets.json` (under `hvantk/skills/<provider>/catalog/`) is the source of truth for that provider's metadata: URLs, version strings, license, citation, release cadence. NEVER restate this content in a skill. Reference the catalog file instead, or query it via `hvantk catalog show <accession>` / `hvantk catalog stats`.

Every provider MUST ship a `plugin.yaml` with `api_version: 2`. Datasets are addressed by compound key `provider:dataset` (e.g., `clinvar:variants`); the loader registers each `DatasetManifest` under that key automatically and resolves its callables lazily on first `get_dataset("provider:dataset")`. The manifest schema is `hvantk/core/plugin/manifest.schema.json`.

Every per-resource `SKILL.md` MUST cover these nine sections, in order, with these exact headings:

1. `## 1. Status & scope`
2. `## 2. Source identity`
3. `## 3. Backend choice + reasoning`
4. `## 4. Raw format & gotchas`
5. `## 5. Output contract`
6. `## 6. hvantk integration points`
7. `## 7. Workflow steps`
8. `## 8. Update playbook`
9. `## 9. Validation contract`

Optional sections (only if they add information not covered above): `## 10. Cross-reference notes`, `## 11. Performance notes`.

## 3. Keying conventions per data domain

- Variants → key `(locus, alleles)`, Hail Table or MatrixTable
- Genes → key `gene_id`, Hail Table
- Proteins → key `interval` or `protein_id`
- Expression (dense, Hail-friendly) → MatrixTable, rows = genes, cols = samples/cells
- Expression (sparse, single-cell) → anndata h5ad
- Lookup / mapping (e.g., HGNC) → small Hail Table or pandas DataFrame, depending on use

## 4. Required helpers

- `create_table_base()` — `hvantk/core/utils/hail_helpers.py`. Optional scaffold for the small set of builders that follow the import → transform → checkpoint → optional TSV-export pattern. Its `import_func` accepts any `Callable[[], hl.Table]` — `hl.import_table` (TSV), `hl.import_vcf().rows()`, or `hl.import_lines` for line-oriented formats like GMT. Most builders build the Table inline instead of using it.
- `cleanup_temp_file()` — `hvantk/core/utils/hail_helpers.py`. Best-effort cleanup of local / Hadoop / S3 / GS temp files. This is the only shared temp helper.
- `init_hail()` — `hvantk/core/utils/hail_context.py`. Idempotent Hail init. Tests use the session-scoped `hail_session` fixture from `conftest.py`.
- AnnData helpers — `annotate_column_summary_ad` in `hvantk/core/models/anndata_utils.py`; `save_anndata` in `hvantk/core/io/anndata_io.py`. Provenance is stamped on the returned Artifact via `ctx.provenance(schema_id=...)` — builders no longer write a separate `hvantk_metadata` dict.
- Plugin runtime — `hvantk/core/plugin/api.py` defines `Provider`, `DatasetSpec`, `DatasetManifest`, and `DriftProbeError`. Tests/CLI consume the populated registry via `get_registry()` in `hvantk/core/plugin/loader.py`.

**Builder contract (current):** plugin builders are functions
`(parsed_input, ctx: BuildContext, **params) -> Artifact` that return
an `AnnotationTable`, `ExpressionMatrix`, or `GeneSet` (see
`hvantk/core/models/`). There is no `(input_path, output_path, overwrite,
export_tsv)` signature — output path and persistence are owned by the
orchestrator, not the builder. The platform invokes builders via
`hvantk.core.plugin.run_builder.run_builder_for_spec(...)`, which runs the
drift probe, constructs the `BuildContext`, calls the builder, validates the
returned artifact's type and `schema_id` against `plugin.yaml`, stamps
source-fingerprint provenance, and saves the artifact. The builder function
name is whatever `plugin.yaml`'s `builder.function` declares — a
`build_<...>` name such as `build_clinvar`, `build_hgnc_gene_lookup`, or
`build_ucsc_cellbrowser`. The old `create_<x>_tb` / `build_<x>_ad` names do
not exist.

NEVER paste these helpers' source into a skill. Reference them by path.

## 5. Builder pattern

- Function naming: `build_<source>` (the exact name is declared in `plugin.yaml`'s `builder.function`).
- Signature shape: `(parsed_input, ctx, **params) -> Artifact`. `parsed_input` is whatever `lifecycle.parse` returned (often a raw path or directory); `ctx` is the platform-supplied `BuildContext`. Common `params`: `reference_genome: str`, plus dataset-specific flags forwarded from `--plugin-arg`.
- The builder returns an `AnnotationTable`, `ExpressionMatrix`, or `GeneSet` wrapper (from `hvantk/core/models/`), stamping provenance via `ctx.provenance(schema_id=...)`. The builder does NOT take an `output_path` / `overwrite` kwarg and does NOT checkpoint itself — `run_builder_for_spec` saves the returned artifact.
- Location: `hvantk/skills/<provider>/builder.py` for single-dataset providers, `hvantk/skills/<provider>/<dataset>/builder.py` for multi-dataset providers.

## 6. Registry registration via plugin.yaml

The plugin loader discovers every `plugin.yaml` under `hvantk/skills/` (plus `hvantk.providers` entry points) and registers each `datasets[]` entry under the compound key `provider:dataset`. There is no `TABLE_BUILDERS` / `MATRIX_BUILDERS` dict and no `create_table_adapter()` to hand-edit — the manifest IS the registration. Example:

```yaml
api_version: 2
name: hgnc
version: 0.1.0
datasets:
  - name: lookup
    domain: mapping
    backend: hail
    artifact_type: AnnotationTable
    schema_id: hgnc-lookup-v1
    builder:
      module: hvantk.skills.hgnc.builder
      function: build_hgnc_gene_lookup
```

`get_registry().get_dataset("hgnc:lookup")` returns the executable `DatasetSpec` (callables resolved lazily); top-level builds run through `run_builder_for_spec`. No `registry.py` edit and no `_apply_plugin_registrations` step is involved.

### Optional: `scores:` — per-predictor training provenance

A dataset that ships *predictor* columns (a pathogenicity score, a trained metric) may declare what each one was fit on:

```yaml
    scores:
      phyloP100way_vertebrate_rankscore: {trained_on: []}
      CADD_raw_rankscore:                {trained_on: [simulated]}
      REVEL_rankscore:                   {trained_on: [HGMD, ClinVar]}
```

The plugin author declares this **once**, because circularity is a property of the data, not of any one analysis — a predictor trained on ClinVar will correlate with a ClinGen-derived label partly because it was trained on genes like those, and no statistical filter can detect that (filtering *rewards* it). Downstream consumers inherit the declaration instead of re-reviewing dozens of scores per cohort; `hvantk`'s rerank uses it to split a run into `clean` and `all` arms and report the gap.

Three states, and the difference matters:

| `trained_on` | meaning |
|---|---|
| `[]` | fit on nothing label-derived (e.g. pure conservation) — never conflicts |
| `[SourceA, …]` | fit on those sources; conflicts with a label derived from any of them |
| score omitted | **unknown**, treated conservatively — usable, but never counted as clean |

Source names are free strings; consumers compare them through an equivalence map (`ClinVar`, `ClinGen`, `GenCC`, `HGMD`, `OMIM` are all curated disease databases and conflict with one another). Declare the sources as the authors describe them and let the consumer's map do the grouping. `scores:` is optional — a dataset with no trained predictors should omit it entirely rather than declare an empty block.

## 7. CLI command pattern

Per-provider CLI lives in `hvantk/skills/<provider>/cli.py` (single-dataset) or `hvantk/skills/<provider>/<dataset>/cli.py` (multi-dataset). The manifest's `cli:` block registers Click commands at top-level discovery time:

```yaml
cli:
  - command: hgnc-download
    module: hvantk.skills.hgnc.cli
    function: download_cmd
```

Top-level data-build invocations go through:

```bash
hvantk reprocess <provider>:<dataset> --raw-dir <dir> --output <out> [--plugin-arg KEY=VALUE]
```

This resolves the manifest via `get_registry().get_dataset(...)`, runs `lifecycle.download` / `lifecycle.parse` / `builder` in sequence, and stamps provenance. Build-time kwargs flow through `--plugin-arg KEY=VALUE`. Provider names use hyphens (e.g., `gtex-eqtl`, `gnomad-metrics`, `ucsc-cellbrowser`, `cosmic-cgc`, `uniprot-ptm`). The legacy `mktable` / `mkmatrix` commands have been retired.

## 8. Test pattern

- Tests live next to the code: `hvantk/skills/<provider>/tests/` (single-dataset) or `hvantk/skills/<provider>/<dataset>/tests/` (multi-dataset).
- Round-trip test file: `test_builder.py` (Hail-backed providers) or `test_<dataset>.py` (anndata providers).
- Mark with `@pytest.mark.hail` if Hail is required. Use the `hail_session` fixture.
- Fixtures: `tests/testdata/raw/<dataset>/`. Snapshots: `tests/snapshots/`.
- Assert against snapshots with `hvantk.tests._snapshot_utils`. The `--regenerate-snapshots` flag rewrites snapshots in place.

## 9. Validation contract

Every per-resource `SKILL.md` MUST declare these paths, which MUST match the `tests:` block in `plugin.yaml`. All paths are resolved relative to the plugin folder, and normally live inside it:

- `fixture` — input file or directory used by the round-trip test. Normally plugin-local (`tests/testdata/raw/<dataset>/`). A fixture that is genuinely shared with cross-cutting tests — `hvantk/tests/test_plugin_conformance.py`, or an integration test in another package — instead lives in the repo-level tree and is referenced in place as `../../tests/testdata/raw/<dataset>`, rather than being duplicated per consumer. Say which form applies in the plugin's `SKILL.md`, since the two are not interchangeable. `dbnsfp`, `ensembl_gene`, `gevir` and `gnomad_metrics` use the shared form.
- `schema_snapshot` — `tests/snapshots/schema.json`
- `row_snapshot` — `tests/snapshots/sample_rows.json`. Keys used to select snapshot rows must be unique-in-table — `_snapshot_utils.collect_sample_rows` does not deduplicate, so a duplicated key yields non-deterministic snapshots. For builders that legitimately produce multi-row keys (e.g., GWAS Catalog), maintain `tests/snapshots/sample_keys.json` listing the singleton-key subset to sample.
- `drift_fingerprint` — `tests/drift_fingerprint.json` (the expected fingerprint; see § 12).
- `command` — the pytest invocation (append `-m hail` only when the test is Hail-marked).

## 10. Hard guardrails

- NEVER invent Hail field names. Read the schema from a real run.
- NEVER invent VCF/TSV column names. Read the file header first.
- NEVER assume catalog content. Read the plugin's `catalog/datasets.json` (or run `hvantk catalog show <accession>`).
- NEVER paste code from a builder into a skill. Reference the file path.
- When uncertain, READ existing code (cite which file).

## 11. Out of scope for any skill

- Hail context init. Tests use `hail_session`; runtime uses `init_hail()`.
- Cross-resource utilities. Genome/locus helpers (`contig_recoding`) live in `hvantk/core/utils/genome.py`; QTL helpers (`parse_gtex_variant_id`, `strip_ensembl_version`) in `hvantk/core/utils/qtl_helpers.py`. Gene-ID mapping is owned by `HGNCGeneCatalogStreamer` in `hvantk/skills/hgnc/streamers.py` (the old `GeneMapper` / `gene_mapper.py` / `gene_aliases.py` are retired).
- "How to use the product" — analytical guidance is downstream.

## 12. Drift probe contract

Each dataset declares a `drift_probe.module` + `function` in `plugin.yaml`. The probe is a zero-arg callable returning a dict with this exact shape:

```python
{
    "probe_version": int,        # bump on probe-logic change
    "source_version": str | None, # upstream version string (Last-Modified, release tag, etc.)
    "headers":  {"<file>": [str, ...]},   # column or section headers from the live source
    "checksums": {"<file>": str},          # sha256 over the bytes used to derive `headers`
    "fetched_at": str,            # ISO-8601 UTC timestamp
}
```

The expected fingerprint lives at `hvantk/skills/<provider>/[<dataset>/]tests/drift_fingerprint.json`. `hvantk drift <provider:dataset>` compares the live probe output against this file. Update the fingerprint when an intentional upstream change has been validated; do not silently regenerate it in the same PR as a behavioural change.

**Sharing one `drift_fingerprint` across datasets is meaningful, not a shortcut.** Datasets that point at the same baseline *and* the same probe are declaring that they have **one** drift signal between them, not one each — which is the honest declaration when the probe is provider-scoped. `ucsc-cellbrowser` is the worked example: `default`, `adult-ctx` and `dev-ctx` are distinct *schema* variants (their obs cell-type column is `celltype`, `Class` and `Type_v2` respectively, so each earns its own snapshot), but `fetch_fingerprint()` takes no arguments and fingerprints the provider-wide catalog, so all three always report identically.

`hvantk drift` probes such a group once and fans the result out — every dataset still gets its own report entry — and the drift workflow opens a single PR per signal. Without that, one upstream event produced three identical PRs whose branches all wrote the same file, so merging any one made the others conflict.

Two datasets may **not** share a baseline while declaring *different* probes: each would overwrite the other's file, and whichever regenerated last would define "clean" for both. `hvantk plugins validate` rejects that.

### Automated drift workflow

A scheduled GitHub Actions workflow (`.github/workflows/drift.yml`) runs `hvantk drift --all --json` daily at 06:00 UTC. For each plugin reporting `status: drifted`, the workflow:

1. Branches `drift/<provider>-<dataset>` from the base branch (`env.BASE_BRANCH`, defaulting to `dev`) — or `drift/<provider>` when several of that provider's datasets share one drift signal, so the group gets a single PR.
2. Regenerates `drift_fingerprint.json` via `hvantk drift --regenerate <provider:dataset>`.
3. Opens (or updates) a draft PR via `gh pr create` / `gh pr edit`, with the structured diff embedded in the body and the regenerated fingerprint already committed. If the plugin's `plugin.yaml` declares `maintainers:` whose entries look like GitHub handles, those handles are `cc`'d in the PR body.

The branch is bot-owned and uses `--force-with-lease`, so if the same dataset drifts again before the previous PR is merged the same branch is updated in place rather than spawning a new PR. Datasets reporting `status: probe_failed` are logged to the job summary but never trigger a PR -- those are infrastructure failures, not data drift.

An agent or human reviews the PR to decide whether the change is a compatible upstream update (just merge the snapshot bump), a breaking schema change (also update `builder.py`), or a spurious probe difference (fix the probe).

The workflow's `workflow_dispatch` trigger accepts a `dry_run` input that runs the helper in `--dry-run` mode, so an operator can validate the workflow plumbing without making real commits.

## 13. Lifecycle stages

Manifests MAY declare `lifecycle.download` and `lifecycle.parse` callables. `hvantk reprocess <provider:dataset>` chains:

1. `lifecycle.download` (if declared) — fetch raw inputs into a working dir.
2. `lifecycle.parse` (if declared) — normalise raw inputs into the builder's expected layout.
3. `builder` — produce the Hail Table or AnnData artifact (always required).
4. `drift_probe` — run a post-build drift check against the committed fingerprint (warning, not failure, unless `--strict` is passed).

Both lifecycle stages are optional; a download-only provider (e.g., a static URL) may omit `parse`, and a vendor-supplied tarball may omit `download`. When present, each is `(module, function)` resolved lazily and surfaced on the `DatasetSpec` (as `download_fn` / `parse_fn`) for the `reprocess` runner.
