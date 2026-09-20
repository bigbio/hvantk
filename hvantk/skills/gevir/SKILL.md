---
name: hvantk:resource-gevir
description: GeVIR gene-level intolerance-to-variation percentile ranks (GeVIR/VIRLoF), keyed by gene_id.
status: provisional
backend: hail
domain: genomics
---

# gevir

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

Provisional. Covers the GeVIR (Gene Variation Intolerance Rank) metrics Hail Table builder (`gevir:metrics`), which imports the per-gene GeVIR/VIRLoF percentile-rank TSV derived from Abramovs, Brass & Tassabehji 2020 (PMID 31873297, *Nature Genetics* 52(1):35-39, DOI 10.1038/s41588-019-0560-2), Supplementary Table 2.

Out of scope for this skill (per `_conventions` § 11):
- Downloading or deriving the raw TSV — no `lifecycle.download` is declared (see § 4); upstream files are materialized externally.
- Cross-resource gene-ID mapping — owned by `HGNCGeneCatalogStreamer` (`hvantk/skills/hgnc/streamers.py`).
- Downstream consumers. Those reference the built table by path.

## 2. Source identity

GeVIR ranks 19,361 protein-coding genes by intolerance to functional variation, derived from the density and spatial distribution of protein-coding variants observed across ~138,632 gnomAD exome and genome sequences (gnomAD v2, GRCh37; the metric itself is build-agnostic). It is a gene-level metric — **not** a variant-level pathogenicity score.

Distributed as Supplementary Table 2 of the paper, served as the article's MOESM3 object from Springer's static-content CDN (`.xlsx`, 10,270,511 bytes per the committed drift fingerprint) — the only public one of the paper's six MOESM slots (the rest answer 403). The authors' code repository (`github.com/gevirank/gevir`) ships analysis code only; its `tables/` directory holds a placeholder file and is not a usable source.

Catalog entry: `hvantk/skills/gevir/catalog/datasets.json`, accession `GeVIR_v1.0` — `license: "Academic use"`, `genome_build: "GRCh37"`, `files: [{"path": "gevir_metrics_pmid31873297.tsv.bgz", ...}]`. That file is the **derived** TSV the builder consumes, distinct from the `.xlsx` the drift probe watches (§ 4).

## 3. Backend choice + reasoning

`hail`. Gene-keyed table (key `gene_id`, per `_conventions` § 3), imported with `hl.import_table(..., key="gene_id")` — consistent with the rest of the gene-metric ecosystem this table is joined against.

## 4. Raw format & gotchas

All facts below are from `hvantk/skills/gevir/builder.py` and `hvantk/skills/gevir/drift_probe.py`.

- **The probe target and the build input are different files.** The drift probe HEADs the upstream `.xlsx` workbook (`GEVIR_SUPPLEMENTARY_URL`); the builder reads a bgzipped TSV (`gevir_metrics_pmid31873297.tsv.bgz`, ~2,000,000 bytes per the catalog entry) manually derived from sheet `table_2` of that workbook. No downloader or parse step performs that extraction — the catalog `files` entry describes the derived TSV precisely because fetching the upstream URL does not yield it without an extract-and-convert step.
- `_resolve_gevir_path()`: `hvantk reprocess` hands the builder the raw *directory* for datasets with no `lifecycle.parse` (`raw_dir` is forwarded verbatim as `parsed_input`), and `hl.import_table` cannot read a directory. A directory input is resolved by globbing it for the single non-hidden file it contains, raising `ValueError` if zero or more than one candidate is found. A path already pointing at a file is returned unchanged, which is what keeps direct `build(<file>)` calls in tests/snapshots working.
- Import: `hl.import_table(paths=..., impute=True, min_partitions=100, key="gene_id")` — type inference is **enabled** here, unlike `hgnc`/`cosmic_cgc`, which import as all-string.
- No field renaming or value transforms beyond an optional `fields` selection — the builder trusts the TSV header names as-is.

## 5. Output contract

`AnnotationTable` (`schema_id="gevir-metrics-v1"`), keyed by `gene_id` (Ensembl gene ID). Schema per the committed snapshot (`hvantk/skills/gevir/tests/snapshots/schema.json`):

| field | type | description |
| --- | --- | --- |
| `gene_id` | str (key) | Ensembl gene identifier (`ENSG...`). |
| `gnomad_gene_name` | str | HGNC gene symbol as used in gnomAD. |
| `canonical_transcript` | str | Ensembl canonical transcript (`ENST...`) used to compute the metric. |
| `gevir_pct` | float64 | GeVIR percentile rank (0-100). Lower = more intolerant/constrained. |
| `virlof_pct` | float64 | VIRLoF percentile rank (0-100), combining GeVIR with the gnomAD LOEUF constraint metric; lower = stronger overall intolerance. |

Live release row count: 19,361 genes (§ 2). The committed row snapshot (`tests/snapshots/sample_rows.json`) samples 6 real genes from the fixture build (`ENSG00000000003`, `ENSG00000000005`, `ENSG00000000419`, `ENSG00000000457`, `ENSG00000000460`, `ENSG00000000938`).

## 6. hvantk integration points

- Manifest: `hvantk/skills/gevir/plugin.yaml` — dataset `gevir:metrics`, `artifact_type: AnnotationTable`, `schema_id: gevir-metrics-v1`, `catalog: catalog/datasets.json`. No `lifecycle:` or `cli:` block declared.
- Builder: `build_gevir_metrics` in `hvantk/skills/gevir/builder.py`. Signature `(parsed_input, ctx, **params) -> AnnotationTable`. Recognised `params`: `fields` (optional list to select).
- Drift probe: `fetch_fingerprint` in `hvantk/skills/gevir/drift_probe.py` (`PROBE_VERSION = 2`). Issues a single HEAD against `GEVIR_SUPPLEMENTARY_URL` with `Accept-Encoding: identity`, comparing `Content-Length` and a normalized `ETag` (Springer serves it as an MD5 content digest) under `headers`; `checksums` stays empty because the probe never fetches a body and so has no schema signal to offer. `Last-Modified` is demoted to `informational` (excluded from drift comparison), following the `hgnc` precedent where every regeneration moved only the timestamp. Fails closed if either validator is missing/empty, or if the response is not identity-encoded (a compressing proxy would make `Content-Length` describe the compressed body).
- No `streamers.py` exists for this plugin (confirmed absent from the plugin directory) and no `cli.py`.
- Catalog: `hvantk/skills/gevir/catalog/datasets.json` (`GeVIR_v1.0`).

## 7. Workflow steps

1. Confirm the raw derived TSV is present (`gevir_metrics_pmid31873297.tsv.bgz`). There is no downloader, so it must be sourced externally by extracting sheet `table_2` from the MOESM3 workbook and BGZF-compressing it (see § 8).
2. Build: `hvantk reprocess gevir:metrics --raw-dir <dir-containing-the-tsv> --output <out.ht>`, optionally `--plugin-arg fields=<comma-list>`. A directory input works because `_resolve_gevir_path` globs it for the single file it contains.
3. Sanity-check the output: table keyed by `gene_id`, ~19,361 rows for a full build, `gevir_pct`/`virlof_pct` within `[0, 100]`.
4. Run the snapshot round-trip test (§ 9); regenerate with `--regenerate-snapshots` after an intentional change and review the diff before committing.

## 8. Update playbook

Triggered by an erratum or updated ESM object on the article (rare — the supplementary data is publication-fixed) or by a need to re-derive the local TSV.

1. `hvantk drift gevir:metrics` compares the live `Content-Length`/`ETag` against `tests/drift_fingerprint.json` (currently `content_length: "10270511"`, `etag: "6423adf134a669acc357f619d1162009"`). A drift here means only that the **upstream `.xlsx`** changed — it says nothing about whether the derived TSV was re-extracted, which is a separate, manual step (§ 4).
2. If sheet `table_2`'s layout changes, the derivation (extract + BGZF-compress; not scripted in this repo) must be re-run to regenerate the TSV, and the fixture plus snapshots regenerated from it.
3. Regenerate snapshots: `pytest hvantk/skills/gevir/tests/test_builder.py -m hail --regenerate-snapshots`, review the diff, commit.
4. A real downloader for this source (extract-and-convert, not a thin URL fetch) is a recommended follow-up per the project's downloader decision framework, but is not yet implemented.

## 9. Validation contract

Declared in `plugin.yaml`'s `tests:` block (`hvantk/skills/gevir/plugin.yaml`), paths resolved relative to the plugin folder:

- `fixture`: `../../tests/testdata/raw/gevir` — the **shared repo-level form**, not a plugin-local fixture. Per `_conventions` § 9, `gevir` is one of four datasets (alongside `dbnsfp`, `ensembl_gene`, `gnomad_metrics`) that use this form because the fixture is shared with cross-cutting tests rather than duplicated. It resolves to `hvantk/tests/testdata/raw/gevir/gevir_metrics_pmid31873297.tsv.bgz`.
- `schema_snapshot`: `tests/snapshots/schema.json` (committed — see § 5)
- `row_snapshot`: `tests/snapshots/sample_rows.json` (committed — see § 5)
- `drift_fingerprint`: `tests/drift_fingerprint.json` (committed — see § 8)
- `test_command`: `pytest hvantk/skills/gevir/tests -m hail`

All four artifacts are shipped. `gevir` is **not** on the `KNOWN_INCOMPLETE` ledger in `hvantk/tests/test_plugin_contract_artifacts.py`; that file's header comment records that `gevir` "shipped its drift fingerprint as part of the gevir plugin-review work, so it left this list."

Tests: `hvantk/skills/gevir/tests/test_builder.py::test_gevir_snapshot_round_trip` (hail-marked snapshot round trip against the fixture) plus three unmarked unit tests for `_resolve_gevir_path` (file-path passthrough, single-file-directory resolution, `ValueError` on an ambiguous directory). `hvantk/skills/gevir/tests/test_drift_probe.py` runs 9 offline (`requests_mock`) tests: fingerprint shape, `checksums` staying empty, `Last-Modified` demotion to `informational`, weak/strong ETag normalization to the same value, fail-closed on a missing content signal, fail-closed on an ETag-only response, fail-closed on an empty ETag, the identity-encoding request header, and fail-closed on a compressed (`Content-Encoding: gzip`) response.
