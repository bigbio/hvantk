---
name: hvantk:resource-alphagenome
description: AlphaGenome per-variant deep-learning effect predictions (expression, chromatin, other molecular phenotypes), fetched live from a credentialed API and built into a Hail Table.
status: provisional
backend: hail
domain: genomics
---

# alphagenome

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

Provisional. Covers the `alphagenome:predictions` builder, which drives the AlphaGenome deep-learning API (Google DeepMind) for a given set of input variants and produces a Hail Table of per-variant effect predictions keyed by `(locus, alleles)`.

Out of scope for this skill (per `_conventions` § 11):
- Provisioning API credentials. Resolved by `load_config` in `hvantk/skills/alphagenome/pipelines.py` (config `api.key` > `ALPHAGENOME_API_KEY` env var > error).
- Downloading/acquiring the raw AlphaGenome SDK. There is no `lifecycle.download` in `plugin.yaml` — this is a credentialed live prediction service, not a static file.
- Per-modality Hail Table assembly of the raw prediction structure. `AlphaGenomePipeline._assemble_outputs` (`hvantk/skills/alphagenome/pipelines.py`) currently merges checkpointed batches into a consolidated `predictions.json`; deferred until the SDK response structure is validated (see module docstring in `pipelines.py`).

## 2. Source identity

AlphaGenome is a deep learning model from Google DeepMind that predicts the functional effects of genetic variants on gene expression, chromatin accessibility, and other molecular phenotypes at nucleotide resolution, served through a credentialed prediction API (`hvantk/skills/alphagenome/pipelines.py` module docstring).

There is no static upstream file: predictions are generated on demand per variant, per interval, against a live model. `source.catalog_ref: alphagenome` is declared in `plugin.yaml`, but this skill does not restate catalog content (per `_conventions` § 2) — query it via `hvantk catalog show alphagenome`.

The AlphaGenome Python SDK (the client this builder imports — `from alphagenome.data import genome`, `from alphagenome.models import dna_client`, `pipelines.py` `_import_alphagenome`) is published openly on PyPI and is what the drift probe fingerprints (see § 3).

## 3. Backend choice + reasoning

`hail`, per `plugin.yaml` (`backend: hail`). The builder produces a Hail Table keyed by `(locus, alleles)` (`_conventions` § 3 variant-domain convention), consistent with every other variant-level annotation source in the toolkit, even though the underlying data originates from per-variant API calls rather than a bulk file import.

## 4. Raw format & gotchas

- `parsed_input` (the builder's first positional argument) must be either a Hail Table path ending in `.ht`, or a TSV with `chrom`/`pos`/`ref`/`alt` columns — both handled in `_run_alphagenome_pipeline` (`hvantk/skills/alphagenome/builder.py`). A `.ht` input must already contain `locus` and `alleles` fields, or the builder raises `ValueError`. A TSV missing any of `chrom`/`pos`/`ref`/`alt` also raises `ValueError` naming the missing columns.
- `config_path` is a **required** `**params` kwarg — no default. It must point at an AlphaGenome YAML config with two required top-level sections, `api` and `ontology` (`pipelines.py` `_REQUIRED_SECTIONS`). The fixture shape is `hvantk/skills/alphagenome/tests/testdata/alphagenome_config.yaml`: `api.key`/`api.max_retries`/`api.retry_backoff`/`api.request_timeout`, `ontology.terms` (UBERON ontology terms), `ontology.output_types` (e.g. `RNA_SEQ`, `CHROMATIN` — validated at stream time against `alphagenome.models.dna_client.OutputType`, `pipelines.py` `stream()`), and an optional `intervals` block (`default_size`, `adaptive`, `adaptive_max_size`, `density_window`), defaulted from `ALPHAGENOME_DEFAULT_INTERVAL_SIZE` (1,048,576 = 1 Mbp) and `ALPHAGENOME_DEFAULT_DENSITY_WINDOW` (50,000 = 50 kb) in `hvantk/skills/alphagenome/shared/constants.py` if omitted.
- `no_resume` (bool, default `False`) is the only other supported `**params` key; `output_path`/`overwrite` are explicitly stripped from `params` before being forwarded to the pipeline (`build_alphagenome_predictions`, `builder.py`) since output handling belongs to the orchestrator, not the builder.
- Adaptive interval grouping (`compute_intervals` / `_compute_adaptive_intervals`, `pipelines.py`): nearby variants on the same chromosome within `density_window` are batched into one API call to reduce request count; variants beyond `adaptive_max_size` apart are split into sub-groups. Disable via `intervals.adaptive: false` in the config to get one interval per variant instead.
- `RateLimitedCaller` (`pipelines.py`) retries transient errors (429/5xx/timeout/connection) with exponential backoff and jitter, and enters a 60-second cooldown after 3 consecutive rate-limit hits (`_COOLDOWN_THRESHOLD`, `_COOLDOWN_SECONDS`). Failed variants (after retries exhausted) are recorded via `CheckpointManager.record_failed_variant` rather than aborting the whole run.
- `CheckpointManager` (`pipelines.py`) writes per-batch JSON files under `<output_dir>/_checkpoints/` and a `state.json` tracking completed intervals, enabling resumption; `no_resume=True` clears this state before starting.
- Loading a `.ht` input caps at `_MAX_VARIANTS_COLLECT = 100_000` rows (`AlphaGenomePipeline._load_variants_from_hail_table`) — it calls `.collect()`, materializing all rows on the driver, and raises `ValueError` above that cap.
- Currently only single-nucleotide/simple ref/alt pairs are exercised; the builder does no allele normalization of its own.

## 5. Output contract

Hail Table keyed by `(locus, alleles)`, per `_conventions` § 3 variant convention. Builder returns an `AnnotationTable` (`hvantk.core.models.AnnotationTable`) via `AnnotationTable.from_hail(ht, provenance=ctx.provenance(schema_id="alphagenome-v1"))` (`build_alphagenome_predictions`, `builder.py`).

Row-level schema is not fixed: prediction output fields depend on the AlphaGenome config (`ontology.output_types`) and the model/SDK version in use at build time, so this skill does not enumerate columns (see § 9 for why no schema snapshot exists). The current builder emits a minimal locus/alleles-keyed table from the input variants and checkpoints the raw prediction JSON alongside it (`_run_alphagenome_pipeline`, `builder.py`); full per-modality field assembly into the Hail Table is deferred (see § 1).

## 6. hvantk integration points

- Plugin manifest: `hvantk/skills/alphagenome/plugin.yaml` (`api_version: 2`, dataset `predictions`, `artifact_type: AnnotationTable`, `schema_id: alphagenome-v1`). The loader resolves it via `get_registry().get_dataset("alphagenome:predictions")`.
- Builder: `build_alphagenome_predictions` in `hvantk/skills/alphagenome/builder.py`. Signature: `(parsed_input, ctx, **params) -> AnnotationTable`.
- Pipeline: `AlphaGenomePipeline` in `hvantk/skills/alphagenome/pipelines.py`, driven internally by `_run_alphagenome_pipeline` (`builder.py`) via `pipeline.setup()` → iterate `pipeline.stream()` → `pipeline.teardown()`.
- Constants: `ALPHAGENOME_DEFAULT_INTERVAL_SIZE`, `ALPHAGENOME_DEFAULT_DENSITY_WINDOW`, `ALPHAGENOME_DEFAULT_RETRY_BACKOFF`, `ALPHAGENOME_DEFAULT_MAX_RETRIES`, `ALPHAGENOME_DEFAULT_REQUEST_TIMEOUT` in `hvantk/skills/alphagenome/shared/constants.py`.
- Drift probe: `fetch_fingerprint` in `hvantk/skills/alphagenome/drift_probe.py`, wired via `plugin.yaml`'s `drift_probe:` block; drives `hvantk drift alphagenome:predictions`.
- No CLI `lifecycle.download`/`lifecycle.parse` blocks and no `cli:` block are declared in `plugin.yaml` — there is no built-in downloader.
- Existing tests: `hvantk/skills/alphagenome/tests/test_alphagenome.py` (registration-only + skipped round-trip), `test_pipelines.py` (unit tests for `load_config`, interval computation, etc. — mocked, no live API), `test_drift_probe.py` (offline via `requests_mock`).

## 7. Workflow steps

When invoked to build or refresh an AlphaGenome predictions table:

1. **Confirm API credentials are available.** Either `api.key` in the config YAML or the `ALPHAGENOME_API_KEY` environment variable must resolve (`load_config`, `pipelines.py`); this workflow does not provision credentials (§ 1).
2. **Confirm `parsed_input` shape.** A `.ht` with `locus`/`alleles`, or a TSV with `chrom`/`pos`/`ref`/`alt` columns (§ 4). Row count matters if using a `.ht`: `_MAX_VARIANTS_COLLECT = 100_000`.
3. **Confirm the config YAML** has `api` and `ontology` sections and valid `ontology.output_types` (validated against `alphagenome.models.dna_client.OutputType` at stream time — an invalid value raises before any API calls are made).
4. **Build** via `hvantk reprocess alphagenome:predictions --raw-dir <dir> --output <out.ht> --plugin-arg config_path=<path/to/config.yaml> [--plugin-arg no_resume=true]`.
5. **Expect API cost and latency.** Each interval issues live API calls with retry/backoff; large variant sets should rely on checkpoint resumption (`CheckpointManager`) rather than restarting from scratch — re-run with the same `output_path` and `no_resume=False` (default) to resume.
6. **Sanity-check the output.** Confirm the table is keyed by `(locus, alleles)` and row count matches the input variant count minus any `failed_variants` recorded in `<output_dir>/_checkpoints/state.json`.

## 8. Update playbook

Triggered when the AlphaGenome SDK publishes a new release (per the drift probe, § 3) or when the model/API behavior is suspected to have changed server-side.

1. Run `hvantk drift alphagenome:predictions` (or wait for the fortnightly automated drift workflow) to check whether the PyPI SDK release stream has moved since the committed `tests/drift_fingerprint.json`.
2. On a schema-risk signal (SDK version bump), re-read the drift probe's module docstring (`drift_probe.py`): it explicitly documents that this probe detects a new SDK release but **cannot** detect a server-side model update shipped without an SDK release — that gap is a property of the service, not a probe bug.
3. If the new SDK release changes the `alphagenome.data.genome` / `alphagenome.models.dna_client` response shape, update `_serialize_value`/`_serialize_prediction` (`pipelines.py`) accordingly and re-validate with a manual smoke test (no fixture exists to automate this — § 9).
4. Regenerate the fingerprint with `hvantk drift --regenerate alphagenome:predictions` once the change is validated; do not regenerate in the same PR as a behavioral change (`_conventions` § 12).

## 9. Validation contract

Declared in `plugin.yaml`'s `tests:` block (all paths plugin-relative under `hvantk/skills/alphagenome/`):

- `fixture`: `tests/testdata/raw/alphagenome`
- `schema_snapshot`: `tests/snapshots/schema.json`
- `row_snapshot`: `tests/snapshots/sample_rows.json`
- `drift_fingerprint`: `tests/drift_fingerprint.json`
- `test_command`: `pytest hvantk/skills/alphagenome/tests`

**Snapshot status:** `alphagenome:predictions` is on the `KNOWN_INCOMPLETE` ledger in `hvantk/tests/test_plugin_contract_artifacts.py`, missing `fixture`, `schema_snapshot`, and `row_snapshot`. Per that file's comment block, AlphaGenome is cause (1) of three: "a credentialed live prediction API" with "no static upstream artifact" that "cannot be snapshotted" — there is no raw file to commit as a fixture, and the row-level schema depends on the live model/config at build time (§ 5), so no fixed schema or row snapshot can exist either. `tests/test_alphagenome.py::test_alphagenome_predictions_round_trip` is `@pytest.mark.skip`'d for this reason ("No fixture available for alphagenome (requires AlphaGenome API access); manual smoke-test only"). Only `drift_fingerprint` is populated (from a live probe run against PyPI) — see § 3 and § 8. `tests/testdata/alphagenome_config.yaml` is a checked-in config fixture used by `test_pipelines.py`'s `TestLoadConfig` unit tests, distinct from the `fixture` path in the `tests:` block (which remains unpopulated).
