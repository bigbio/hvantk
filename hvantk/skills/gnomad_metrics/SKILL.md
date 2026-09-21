---
name: hvantk:resource-gnomad-metrics
description: gnomAD per-gene constraint metrics (pLI, oe_lof/LOEUF, mis_z) — Hail Table keyed by gene_id.
status: provisional
backend: hail
domain: genomics
---

# gnomAD constraint metrics

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

Provisional. Covers the `gnomad-metrics:metrics` builder, which imports the gnomAD per-gene constraint table (pLI, oe_lof/LOEUF, mis_z, and related intolerance statistics) into a Hail Table. Two upstream releases are supported: v2.1.1 `by_gene` (GRCh37, keyed by `gene_id`, the hvantk default) and v4.0 `constraint_metrics` (GRCh38, per-transcript, no `gene_id` column).

Out of scope for this skill:
- The full gnomAD frequency/genotype VCF release (hundreds of GB) — documentation-only per `docs_site/guide/data-sources.md`'s size threshold; not built by this plugin.
- Downstream consumers that join this table against gene-level artifacts.

## 2. Source identity

gnomAD (Genome Aggregation Database, Broad Institute) publishes per-gene loss-of-function constraint tables in the public `gcp-public-data--gnomad` GCS bucket, served over plain HTTPS with no auth (`hvantk/skills/gnomad_metrics/shared/constants.py`, `GNOMAD_RELEASE_BASE_URL`).

> Catalog note: `hvantk/skills/gnomad_metrics/catalog/datasets.json` (accession `gnomAD_v4.1`) describes the **full-genome/exome sites VCF release** (`gnomad.genomes.v4.1.sites.vcf.bgz`, ~150 GB; `gnomad.exomes.v4.1.sites.vcf.bgz`, ~80 GB) — a different, much larger artifact than the small constraint table this plugin actually builds. Do not treat that catalog entry as describing `gnomad-metrics:metrics`'s source files; the authoritative URLs/versions for what this plugin builds live in `shared/constants.py` (`GNOMAD_CONSTRAINT_TABLES`), not the catalog.

Two releases, both hosted under `GNOMAD_RELEASE_BASE_URL`:
- `v2.1.1` (path segment `2.1.1`, no `v` prefix) — `by_gene` (default) and `by_transcript`, flat column names, one row per gene/transcript, includes `gene_id`.
- `v4.0` (path segment `v4.0`) — `constraint_metrics`, one row per transcript, **dotted** column names (`lof.pLI`, `lof.oe_ci.upper` = LOEUF, `mis.z_score`, …), **no `gene_id` column**. gnomAD did not re-release constraint for v4.1, so v4.0 is the newest constraint table.

## 3. Backend choice + reasoning

`hail`, per `plugin.yaml`. Gene-level annotation table keyed by `gene_id` — matches `_conventions` § 3's gene-domain keying convention, and lets consumers join it directly against other Hail-based gene tables.

## 4. Raw format & gotchas

- File: tab-separated, single header row. Confirmed against the shared fixture `hvantk/tests/testdata/raw/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.chr20.tsv.bgz` — header includes `gene`, `transcript`, `gene_id`, `pLI`, `oe_lof`, `constraint_flag`, `brain_expression`, `chromosome`, etc.
- Missing-value handling: `hl.import_table(..., impute=True, ...)` (`builder.py`) — no explicit `missing=` argument, so Hail's default sentinel `"NA"` applies. The fixture's `brain_expression` column carries literal `NA` for missing genes and is imported as `null`, not the string `"NA"`.
- Empty string is a *different* case from `NA`: `constraint_flag` holds an empty field (two adjacent tabs) rather than `NA` for unflagged genes, and is imported as `""` (an empty string), not `null` — do not conflate the two when filtering.
- No column renaming: the builder passes the upstream header straight through (`hl.import_table` with no `.rename`), unlike `hgnc`'s renamed-field pattern.
- Key column: default `key="gene_id"` (v2.1.1 `by_gene`/`by_transcript`). v4.0 `constraint_metrics` has no `gene_id` column, so callers must pass `--plugin-arg key=transcript` (or another present column) — the builder does not infer this automatically.
- **The committed fixture is v2.1.1, whose column names are flat** (78 columns, none containing a dot). The dotted-name advice below therefore applies to a v4.0 build only and is NOT exercised by the round-trip test.
- v4.0's column names are dotted — `lof.pLI`, `lof.oe_ci.upper` (= LOEUF), `mis.z_score` — and v4.0 ships **no `gene_id` column** at all; both facts are documented in `shared/constants.py` (see its module docstring). A dot is not valid in a Python identifier, so once such a column is imported these fields must be subscripted (`ht['lof.pLI']`): attribute access `ht.lof.pLI` reads as field `lof` then `.pLI` and fails. That consequence is Hail's, not something `constants.py` states.
- `parsed_input` may be a file path or a `raw_dir` (`builder.py`): when it is a directory, the builder globs `*.bgz` then `*.tsv` inside it. Zero matches raises `FileNotFoundError`; more than one match raises `ValueError` and refuses to silently pick one — a `raw_dir` must hold exactly one version's constraint file (mirrors the `gevir` builder's fail-loud pattern, per PR #222 review referenced in `builder.py`).
- Import options: `impute=True`, `min_partitions=100` — contrast with `hgnc`, which imports with `impute=False` and keeps every field as a string.
- Optional `fields` param (`list[str]`) selects a column subset post-import via `ht.select(*fields)`.

## 5. Output contract

Hail Table keyed by `gene_id` (v2.1.1) or the caller-supplied `key` (v4.0, typically `transcript`). Builder returns an `AnnotationTable` (`hvantk.core.models.AnnotationTable`) stamped with `schema_id="gnomad-metrics-v1"`. Schema is the source of truth — see `hvantk/skills/gnomad_metrics/tests/snapshots/schema.json`.

Row fields (from the v2.1.1 `by_gene` schema snapshot) fall into: identifiers (`gene`, `gene_id`, `transcript`, `transcript_type`, `transcript_level`, `gene_type`), constraint scores (`pLI`, `pNull`, `pRec`, `oe_lof`, `oe_lof_lower`, `oe_lof_upper` = LOEUF, `oe_mis`, `oe_syn`, `lof_z`, `mis_z`, `syn_z`, `mu_lof`, `mu_mis`, `mu_syn`), observed/expected counts (`obs_lof`, `exp_lof`, `obs_mis`, `exp_mis`, `obs_syn`, `exp_syn`, `obs_het_lof`, `obs_hom_lof`, `exp_hom_lof`, `possible_lof`/`mis`/`syn`), population allele-frequency fields (`classic_caf`, `classic_caf_<pop>`, `p`, `p_<pop>` for `afr`/`amr`/`asj`/`eas`/`fin`/`nfe`/`oth`/`sas`), gene structure (`gene_length`, `cds_length`, `num_coding_exons`, `chromosome`, `start_position`, `end_position`), ExAC cross-references (`exac_pLI`, `exac_obs_lof`, `exac_exp_lof`, `exac_oe_lof`), rank/bin fields (`oe_lof_upper_bin`, `oe_lof_upper_bin_6`, `oe_lof_upper_rank`, `oe_mis_upper_bin`), and misc (`constraint_flag`, `brain_expression`, `n_sites`, `no_lofs`, `defined`, `max_af`). All fields are scalar (`str`/`int32`/`float64`) — no arrays or structs.

## 6. hvantk integration points

- Plugin manifest: `hvantk/skills/gnomad_metrics/plugin.yaml`, dataset key `gnomad-metrics:metrics`.
- Builder: `build_gnomad_metrics_metrics` in `hvantk/skills/gnomad_metrics/builder.py`, signature `(parsed_input, ctx, **params) -> AnnotationTable`; recognized `params`: `fields`, `key`.
- Downloader: `download_gnomad_metrics` / `download_dataset` (lifecycle entry point) / `download_cmd` in `hvantk/skills/gnomad_metrics/cli.py`; CLI command `gnomad-metrics-download` (wired via `plugin.yaml`'s `cli:` block); `--version {v2.1.1,v4.0}` and `--table <name>` select the object.
- Drift probe: `fetch_fingerprint` in `hvantk/skills/gnomad_metrics/drift_probe.py` — HEADs every object in `GNOMAD_CONSTRAINT_TABLES` (all three: v2.1.1 `by_gene`, `by_transcript`, v4.0 `constraint_metrics`) over one shared `requests.Session`, and compares the MD5 `ETag`, `Content-Length`, and `x-goog-generation` per object under `headers` (no `checksums`, since no body is fetched). `Last-Modified` is recorded under `informational` only.
- Constants: `GNOMAD_RELEASE_BASE_URL`, `GNOMAD_CONSTRAINT_TABLES`, `DEFAULT_VERSION`, `DEFAULT_TABLE`, `resolve_table`, `constraint_url`, `constraint_filename` in `hvantk/skills/gnomad_metrics/shared/constants.py`.
- Tests: `hvantk/skills/gnomad_metrics/tests/test_gnomad_metrics.py` (registration + round-trip via `run_builder_for_spec`), `test_builder.py`, `test_download.py`, `test_drift_probe.py`.

## 7. Workflow steps

1. **Confirm the raw file is present.** If absent, run `hvantk download gnomad-metrics --output <path> [--version v2.1.1|v4.0] [--table <name>]`.
2. **If pointing `--raw-dir` at a directory**, ensure it holds exactly one `*.bgz`/`*.tsv` constraint file — the builder raises `ValueError` on more than one candidate rather than guessing.
3. **Build the table**: `hvantk reprocess gnomad-metrics:metrics --raw-dir <dir> --output <out>.ht [--plugin-arg fields='["gene_id","pLI","oe_lof"]'] [--plugin-arg key=transcript]` (the `key` override is required for v4.0).
4. **Sanity-check the output.** Confirm the key matches what was requested (`gene_id` by default), and that `constraint_flag` values are empty strings rather than nulls for unflagged genes (§ 4).
5. **Run the round-trip test** (§ 9) after any builder or fixture change.

## 8. Update playbook

Triggered when a new gnomAD constraint release ships, or the drift probe reports a schema-tier change.

1. `hvantk drift gnomad-metrics:metrics` compares live ETag/Content-Length/`x-goog-generation` per object against `tests/drift_fingerprint.json`. Because the compared surface lives entirely under `headers`, every drift event on this plugin is tiered `"schema"` (never `"routine"`) by `classify_risk` — see `_conventions` § 12 and the reasoning documented in `drift_probe.py`.
2. If a genuinely new release adds/replaces an object path, update `GNOMAD_CONSTRAINT_TABLES` (and `DEFAULT_VERSION`/`DEFAULT_TABLE` if the default should move) in `shared/constants.py` — the drift probe's `OBJECT_PATHS` is derived from that constant, so no separate probe edit is needed.
3. Regenerate the fingerprint with `hvantk drift --regenerate gnomad-metrics:metrics` after validating the change is a genuine upstream republish, not a probe artifact.
4. If the header of the upstream TSV changed shape, regenerate `tests/snapshots/{schema,sample_rows}.json` with `pytest hvantk/skills/gnomad_metrics/tests -m hail --regenerate-snapshots` and review the diff before committing.

## 9. Validation contract

Declared in `plugin.yaml`'s `tests:` block for dataset `gnomad-metrics:metrics`. Per `_conventions` § 9, this plugin uses the **shared repo-level fixture form** (not plugin-local) — the fixture lives under `hvantk/tests/testdata/raw/gnomad/`, shared with `hvantk/tests/test_plugin_conformance.py`, and is referenced in place rather than duplicated:

- `fixture`: `../../tests/testdata/raw/gnomad` (resolved relative to the plugin folder → `hvantk/tests/testdata/raw/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.chr20.tsv.bgz`)
- `schema_snapshot`: `tests/snapshots/schema.json`
- `row_snapshot`: `tests/snapshots/sample_rows.json`
- `drift_fingerprint`: `tests/drift_fingerprint.json`
- `command`: `pytest hvantk/skills/gnomad_metrics/tests -m hail`
