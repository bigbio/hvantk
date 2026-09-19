---
name: hvantk:resource-cosmic-cgc
description: COSMIC Cancer Gene Census (CGC) gene-level cancer catalogue, keyed by gene_symbol or hgnc_id.
status: provisional
backend: hail
domain: genomics
---

# cosmic-cgc

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

Provisional. Covers the COSMIC Cancer Gene Census (CGC) submissions Hail Table builder (`cosmic-cgc:submissions`), which imports a licence-gated COSMIC CGC TSV export into a gene-level annotation table.

Out of scope for this skill (per `_conventions` § 11):
- Downloading the raw file. No `lifecycle.download`/`cli.py` exists for this plugin (confirmed: the plugin folder has no `cli.py`) — acquisition is manual.
- Cross-resource gene-symbol → HGNC resolution logic. Owned by `HGNCGeneCatalogStreamer` (`hvantk/skills/hgnc/streamers.py`); this builder only *consumes* a `gene_catalog` object handed to it.
- Downstream consumers. Those reference the built table by path.

## 2. Source identity

COSMIC Cancer Gene Census, maintained by the Wellcome Sanger Institute (`https://cancer.sanger.ac.uk/census`). A curated catalogue of genes with causal roles in human cancer: tier classification (Tier 1 — strong experimental evidence; Tier 2 — literature-supported), somatic/germline mutation context, tumour types, and role in cancer (oncogene/TSG/fusion).

`plugin.yaml` declares `source.catalog_ref: cosmic-cgc`, but — unlike `gevir` — this plugin ships **no** `catalog/datasets.json` (no `catalog:` key in `plugin.yaml`, and no `catalog/` folder in the plugin directory). There is no catalog entry to reference for URLs/version/license; nothing here should be inferred beyond what the code shows.

The Census *data* page (`cancer.sanger.ac.uk/census`) answers 302 to `/cosmic/login` — it is login- and licence-gated, so acquisition is manual. The release-notes page (`https://cancer.sanger.ac.uk/cosmic/release_notes`, no trailing slash) does not require login and is the only public, probeable artifact (see § 3).

## 3. Backend choice + reasoning

`hail`. Gene-keyed table, consistent with `_conventions` § 3. COSMIC CGC is consumed through `CosmicCGCGeneDiseaseTableStreamer` (§ 6), which subclasses the Hail-based `GeneDiseaseTableStreamer`, so a Hail Table avoids materializing a separate format for that consumer.

## 4. Raw format & gotchas

All facts below are from `hvantk/skills/cosmic_cgc/builder.py` and `hvantk/skills/cosmic_cgc/shared/constants.py`.

- Input is a TSV (`delimiter="\t"`); compression is auto-resolved via `resolve_compression()` (`hvantk/core/utils/file_utils.py`), which detects BGZF/gzip/uncompressed and sets `force_bgz` accordingly. Imported with `hl.import_table(impute=False, min_partitions=10, ...)` — fields stay strings, no type inference.
- Field rename map `COSMIC_CGC_FIELDS` (`shared/constants.py`), applied via `build_rename_map()` (`hvantk/core/utils/table_utils.py`, case/separator-insensitive matching): `"Gene Symbol"→gene_symbol`, `"Name"→gene_name`, `"Entrez GeneId"→entrez_id`, `"Genome Location"→genome_location`, `"Tier"→classification`, `"Hallmark"→hallmark`, `"Chr Band"→chr_band`, `"Somatic"→somatic`, `"Germline"→germline`, `"Tumour Types(Somatic)"→tumour_types_somatic`, `"Tumour Types(Germline)"→tumour_types_germline`, `"Cancer Syndrome"→cancer_syndrome`, `"Tissue Type"→tissue_type`, `"Molecular Genetics"→molecular_genetics`, `"Role in Cancer"→role_in_cancer`, `"Mutation Types"→mutation_types`, `"Translocation Partner"→translocation_partner`, `"Other Germline Mut"→other_germline_mut`, `"Other Syndrome"→other_syndrome`, `"Synonyms"→synonyms`.
- Tier normalization: raw digit-only `classification` values (`"1"`, `"2"`) are rewritten to `"Tier 1"` / `"Tier 2"`; already-prefixed values pass through unchanged.
- `classification_level` (int) is added via `annotate_classification_level()` against `COSMIC_CGC_CLASSIFICATION_LEVELS = ["Tier 1", "Tier 2"]` — ordinal rank, most-confident first; anything unmatched (including missing) sorts last.
- `somatic`, `germline`, `hallmark` are coerced to boolean via `str_to_bool()`: recognises `yes`/`y`/`true`/`1` case-insensitively; missing or anything else → `False`.
- Multi-value fields `tumour_types_somatic`, `tumour_types_germline`, `role_in_cancer`, `mutation_types` are comma-split, stripped, and emptied entries dropped; missing/empty input becomes `[]`, not a missing array.
- Keying: if a `gene_catalog` (a `GeneCatalogStreamer`) build param is supplied, symbols are resolved to `hgnc_id` via `map_ids()`, the `HGNC:` prefix is stripped, unresolved rows are dropped, and the table is keyed by `hgnc_id`. Otherwise the table is keyed by `gene_symbol` (a warning is logged).
- `COSMIC_TISSUE_TYPES` (`E/L/M/O → Epithelial/Lymphoid/Mesenchymal/Other`) exists in `shared/constants.py` but the builder never applies it — `tissue_type` stays as the raw single-letter code in the built table; the decoding map is unused by the builder (only relevant if a consumer wants to decode it).
- No downloader is implemented; upstream files must be externally materialized because COSMIC gates the CGC export behind account login.

## 5. Output contract

`AnnotationTable` (`schema_id="cosmic-cgc-v1"`), one row per gene, keyed by `hgnc_id` (if a `gene_catalog` was supplied) or `gene_symbol` (default). No schema snapshot is committed (§ 9), so the field list below is derived from `COSMIC_CGC_FIELDS` plus the builder's added fields, not from a live schema dump: `gene_symbol`, `gene_name`, `entrez_id`, `genome_location`, `classification` (`"Tier 1"`/`"Tier 2"`), `classification_level` (int), `hallmark` (bool), `chr_band`, `somatic` (bool), `germline` (bool), `tumour_types_somatic` (`array<str>`), `tumour_types_germline` (`array<str>`), `cancer_syndrome`, `tissue_type`, `molecular_genetics`, `role_in_cancer` (`array<str>`), `mutation_types` (`array<str>`), `translocation_partner`, `other_germline_mut`, `other_syndrome`, `synonyms`, plus `hgnc_id` when resolved.

## 6. hvantk integration points

- Manifest: `hvantk/skills/cosmic_cgc/plugin.yaml` — dataset `cosmic-cgc:submissions`, `artifact_type: AnnotationTable`, `schema_id: cosmic-cgc-v1`. No `lifecycle:` or `cli:` block declared.
- Builder: `build_cosmic_cgc_submissions` in `hvantk/skills/cosmic_cgc/builder.py`. Signature `(parsed_input, ctx, **params) -> AnnotationTable`. Recognised `params`: `mutation_context` (`"both"` default, `"somatic"`, or `"germline"`), `min_classification` (one of `COSMIC_CGC_CLASSIFICATION_LEVELS`), `gene_catalog` (a `GeneCatalogStreamer`), `fields`.
- Constants: `COSMIC_CGC_FIELDS`, `COSMIC_CGC_CLASSIFICATION_LEVELS`, `COSMIC_TISSUE_TYPES`, `COSMIC_MUTATION_CONTEXTS` in `hvantk/skills/cosmic_cgc/shared/constants.py`.
- Drift probe: `fetch_fingerprint` in `hvantk/skills/cosmic_cgc/drift_probe.py` (`PROBE_VERSION = 2`). Probes the public release-notes page only — the gated Census data is never fetched — and anchors on `id="v<N>"` HTML attributes, never on prose (a prior prose-matching approach picked up unrelated `Actionability`-product version mentions). Fails closed if the response redirected (the host 302s the trailing-slash variant to a login form) or if no anchors are found.
- Streamer: `CosmicCGCGeneDiseaseTableStreamer` in `hvantk/skills/cosmic_cgc/streamers.py`, subclassing `GeneDiseaseTableStreamer` (`hvantk/core/streamers/gene_disease_table.py`). Adds `filter_by_mutation_context`, `get_geneset_per_tumour_type`/`_role`/`_tissue`, `mutation_context_summary`, `role_summary`, `classification_summary`, `describe`. Because CGC is one row per gene (not per gene-disease assertion), the MONDO/disease-dependent base methods (`get_genes_by_mondo_id`, `categorize_by_ontology`, `get_genes_by_disease`, `get_genes_by_moi`, `get_geneset_per_disease`, `categorize_by_ontology_summary`) are overridden to raise `NotImplementedError`.

## 7. Workflow steps

1. Obtain the raw COSMIC CGC export manually — requires a COSMIC account/login at `cancer.sanger.ac.uk/census`. No in-repo downloader exists.
2. Build: `hvantk reprocess cosmic-cgc:submissions --raw-dir <dir> --output <out.ht>`, optionally with `--plugin-arg mutation_context=somatic|germline|both`, `--plugin-arg min_classification=<level>`, `--plugin-arg fields=<comma-list>`.
3. To resolve `hgnc_id`, pass a `gene_catalog` (an `HGNCGeneCatalogStreamer`) build param via the Python API — `--plugin-arg` values are strings, so a streamer object cannot be passed that way; call `build_cosmic_cgc_submissions(...)` directly for that path.
4. Sanity-check the output: `classification` values are `"Tier 1"`/`"Tier 2"`, `somatic`/`germline`/`hallmark` are true booleans (not strings), and the multi-value fields are arrays.
5. Run the registered-in-loader test (§ 9). No round-trip/snapshot assertion can run — there is no fixture.

## 8. Update playbook

Triggered when COSMIC publishes a new release (the `v101`…`v104`… series) or the CGC export's column layout changes.

1. `hvantk drift cosmic-cgc:submissions` compares the live release index against `tests/drift_fingerprint.json`. Regenerate via `hvantk drift --regenerate cosmic-cgc:submissions` once a genuine new release is confirmed.
2. A moved fingerprint means only that a new COSMIC release exists (`source_version` bumped, e.g. `"v104"` → `"v105"`) — it does **not** confirm the Census gene-level contents changed; per `drift_probe.py`, that is undetectable without an authenticated fetch.
3. If COSMIC's CGC export column headers change, update `COSMIC_CGC_FIELDS` (and `COSMIC_CGC_CLASSIFICATION_LEVELS` / `COSMIC_MUTATION_CONTEXTS` if the tier or mutation-context vocabulary itself changes) in `shared/constants.py`.
4. There is no fixture or snapshot to regenerate (§ 9); verify the builder manually against a freshly (manually) downloaded export.

## 9. Validation contract

Declared in `plugin.yaml`'s `tests:` block (`hvantk/skills/cosmic_cgc/plugin.yaml`), paths plugin-relative:

- `fixture`: `tests/testdata/raw/cosmic-cgc`
- `schema_snapshot`: `tests/snapshots/schema.json`
- `row_snapshot`: `tests/snapshots/sample_rows.json`
- `drift_fingerprint`: `tests/drift_fingerprint.json`
- `test_command`: `pytest hvantk/skills/cosmic_cgc/tests`

**`fixture`, `schema_snapshot`, and `row_snapshot` do not exist and cannot be created.** `cosmic-cgc:submissions` is on the `KNOWN_INCOMPLETE` ledger in `hvantk/tests/test_plugin_contract_artifacts.py` with exactly those three fields missing, because COSMIC's licence forbids redistributing rows — this is cause 2 of the three the ledger's header comment documents ("licence forbids redistributing rows (cosmic-cgc)"), distinct from "not yet seeded." This is a permanent gap, not a to-do.

What does exist and run: `hvantk/skills/cosmic_cgc/tests/test_cosmic_cgc.py::test_cosmic_cgc_submissions_registered` (loader-only — confirms the manifest resolves `cosmic-cgc:submissions` to `AnnotationTable` / `cosmic-cgc-v1`); `test_cosmic_cgc_submissions_round_trip` is explicitly `@pytest.mark.skip`'d, reason: `"No fixture available for cosmic-cgc (COSMIC requires account login); manual smoke-test only"`. `hvantk/skills/cosmic_cgc/tests/test_drift_probe.py` runs 7 offline (`requests_mock`) tests of the release-index probe: anchor-only matching excludes unrelated `Actionability`-product prose versions, editorial prose edits don't move the signal, forward-looking prose can't bump `source_version`, numeric sort, new-release detection, and fail-closed on a login redirect or missing anchors. `tests/drift_fingerprint.json` **is** committed (`source_version: "v104"`, `headers.cosmic-release-index: ["v101", "v102", "v103", "v104"]`) — the release-index signal is a real, running drift check even though no row/schema snapshot can ever exist for this plugin.
