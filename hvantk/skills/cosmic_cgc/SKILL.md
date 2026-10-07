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

Header, delimiter, compression, and value conventions below were confirmed against a
licensed `Cosmic_CancerGeneCensus_v103_GRCh38` export (763 rows) on 2026-10-06 — counts
and conventions only, never rows (the COSMIC licence forbids redistributing those; see
§ 9). Everything else is from `hvantk/skills/cosmic_cgc/builder.py` and
`hvantk/skills/cosmic_cgc/shared/constants.py`.

- **Current (v103+) header**, 21 tab-separated columns, this exact order: `GENE_SYMBOL
  NAME COSMIC_GENE_ID CHROMOSOME GENOME_START GENOME_STOP CHR_BAND SOMATIC GERMLINE
  TUMOUR_TYPES_SOMATIC TUMOUR_TYPES_GERMLINE CANCER_SYNDROME TISSUE_TYPE
  MOLECULAR_GENETICS ROLE_IN_CANCER MUTATION_TYPES TRANSLOCATION_PARTNER
  OTHER_GERMLINE_MUT OTHER_SYNDROME TIER SYNONYMS`. No quoting (confirmed: zero
  double-quote characters across a 200-line sample).
- **Legacy exports** instead used Title-Case headers (`"Gene Symbol"`, `"Chr Band"`, …)
  plus a combined `"Genome Location"` string, a separate `"Entrez GeneId"`, and a
  `"Hallmark"` column — none of which the current export carries at all (confirmed
  absent, not just unmapped) — and had no `COSMIC_GENE_ID`, `CHROMOSOME`,
  `GENOME_START` or `GENOME_STOP`. Each header generation has its own rename map in
  `shared/constants.py`: `COSMIC_CGC_CURRENT_FIELDS` (the 21 columns above) and
  `COSMIC_CGC_LEGACY_FIELDS` (20 columns). `build_rename_map()`
  (`hvantk/core/utils/table_utils.py`) matches case- and separator-insensitively, so
  `Chr Band` and `CHR_BAND` both match either spelling.
- **Header generation.** The builder renames with both maps (`COSMIC_CGC_FIELDS`), then
  takes the generation the header is closest to (fewest missing columns, then fewest
  columns that generation does not define; `COSMIC_CGC_HEADER_GENERATIONS`). If that
  generation lacks any column of its map, the build raises `ValueError` naming the
  generation, the missing columns and the columns it does not define, so a partly
  renamed header (e.g. `GENOME_START` → `GENOME_START_POS`, or a legacy `Genome
  Location` renamed) never gets a schema ID whose fields it lacks. A header that has
  every column plus extra ones builds and keeps the extras, with a warning that names
  them; an extra column that the other generation's map knows is renamed (and
  `Hallmark` is cast to bool).
- Compression: the real `.tsv.gz` is **standard gzip, not BGZF** — confirmed from its
  magic bytes (`FLG=0x08`/FNAME set; no BGZF `"BC"` extra-field subfield).
  `resolve_compression()` (`hvantk/core/utils/file_utils.py`) correctly classifies it as
  plain gzip (`force_bgz=False`), but `hl.import_table()` refuses to read ANY non-BGZF
  `.gz` path at all without `force=True` — regardless of `min_partitions` (confirmed
  empirically with it unset, 1, and 10) — raising `HailException` for exactly this
  input. The builder now passes `force=True` whenever
  `force_bgz` is `False` — confirmed harmless for uncompressed input (partition count is
  unaffected either way) and required for the real export's actual compression.
- Missing-value convention: the empty string. No literal `NA`/`N/A`/`-`/`null` sentinel
  appears anywhere as a field value in the licensed export.
- `TIER`: bare digit `"1"` or `"2"` (590 / 173 of 763 rows respectively).
- `SOMATIC`/`GERMLINE`: `"y"`/`"n"` — **never** `"yes"`, and never empty in the licensed
  export (719/44 and 111/652 of 763 rows respectively). `str_to_bool()`'s broader
  truthy set (`yes`/`y`/`true`/`1`, case-insensitive) still applies correctly since
  `"y"` is a member.
- `CHROMOSOME`: bare `1`-`22` or `X` (no `chr` prefix); no `Y`/`MT` values observed.
- `GENOME_START`/`GENOME_STOP`: numeric or empty (6/763 rows are empty). Both are
  cast to `int32` via `hl.parse_int32()`, which is missing-tolerant — empty or
  non-numeric input becomes missing rather than raising (Hail 0.2 API:
  `functions/string.html#hail.expr.functions.parse_int32`). The builder counts the
  non-empty values that fail to parse (e.g. `1000000.0`) and logs a warning with the
  count for each column that has any. A current header must carry both columns, so the
  cast and the count always run for it; legacy exports carry neither.
- `COSMIC_GENE_ID`: shape `COSG` + 5-6 digits; always populated.
- Tier normalization: raw digit-only `classification` values (`"1"`, `"2"`) are rewritten to `"Tier 1"` / `"Tier 2"`; already-prefixed values pass through unchanged.
- `classification_level` (int) is added via `annotate_classification_level()` against `COSMIC_CGC_CLASSIFICATION_LEVELS = ["Tier 1", "Tier 2"]` — ordinal rank, most-confident first; anything unmatched (including missing) sorts last.
- `somatic`, `germline`, `hallmark` are coerced to boolean via `str_to_bool()`: recognises `yes`/`y`/`true`/`1` case-insensitively; missing or anything else → `False`. `hallmark` is only annotated when the input actually carries a `Hallmark` column (legacy exports) — the current header has none, and the builder does not fabricate it.
- Multi-value fields `tumour_types_somatic`, `tumour_types_germline`, `role_in_cancer`, `mutation_types` are comma-separated, usually `", "` but not always — a handful of rows in the licensed export have a bare `,` with no following space, which is why the builder splits on `,` and strips each token rather than splitting on the literal `", "`. Tokens are stripped and emptied entries dropped; missing/empty input becomes `[]`, not a missing array. `ROLE_IN_CANCER` vocabulary: `oncogene`, `TSG`, `fusion`. `MUTATION_TYPES` vocabulary includes short codes (`Mis`, `N`, `F`, `D`, `A`, `O`, `S`, `T`) and a couple of multi-word tokens (e.g. one with an embedded period and space) — real COSMIC data, not a parsing artifact; the comma-split-then-strip approach handles both without special-casing.
- `TISSUE_TYPE` vocabulary: `E`/`L`/`M`/`O`, decoded by `COSMIC_TISSUE_TYPES` in `shared/constants.py` (`Epithelial`/`Lymphoid`/`Mesenchymal`/`Other`), though the builder never applies that decode — `tissue_type` stays as the raw code. **Known gap** (not fixed by this update, flagged for whoever picks it up next): `TISSUE_TYPE` is itself sometimes multi-valued in the licensed export — comma-separated, same convention as the four fields above, roughly 18% of rows carry 2+ codes — but the builder leaves it as a raw scalar string rather than splitting it into an array. `CosmicCGCGeneDiseaseTableStreamer.get_geneset_per_tissue()` (§ 6) groups by it as a scalar, so splitting it here would need a matching streamer change; out of scope for this builder-focused update.
- `SYNONYMS` is comma-separated but **not** split into an array by the builder — it stays a raw string. Its separator spacing is the opposite of the four fields above: a bare `,` with no following space is the common case, `", "` the exception.
- Keying: if a `gene_catalog` (a `GeneCatalogStreamer`) build param is supplied, symbols are resolved to `hgnc_id` via `map_ids()`, the `HGNC:` prefix is stripped, unresolved rows are dropped, and the table is keyed by `hgnc_id`. The builder logs how many symbols resolved, as a warning naming the dropped count (and up to ten examples) when any did not; if the table has symbols but the catalog resolves none of them, the build raises `ValueError` instead of writing an empty table. Otherwise the table is keyed by `gene_symbol` (a warning is logged).
- No downloader is implemented; upstream files must be externally materialized because COSMIC gates the CGC export behind account login.

## 5. Output contract

`AnnotationTable` (`schema_id="cosmic-cgc-v2"` for v103+ headers, or `cosmic-cgc-legacy-v1` for legacy headers; a header missing a column of its generation raises, § 4), one row per gene, keyed by `hgnc_id` (if a `gene_catalog` was supplied) or `gene_symbol` (default). Field list below is the committed `tests/snapshots/schema.json` (§ 9), built from the current (v103+) header via the synthetic fixture with a `gene_catalog` supplied (hence `hgnc_id` present and keyed):

`gene_symbol` (str), `gene_name` (str), `cosmic_gene_id` (str), `chromosome` (str), `genome_start` (int32), `genome_stop` (int32), `chr_band` (str), `somatic` (bool), `germline` (bool), `tumour_types_somatic` (`array<str>`), `tumour_types_germline` (`array<str>`), `cancer_syndrome` (str), `tissue_type` (str), `molecular_genetics` (str), `role_in_cancer` (`array<str>`), `mutation_types` (`array<str>`), `translocation_partner` (str), `other_germline_mut` (str), `other_syndrome` (str), `classification` (str, `"Tier 1"`/`"Tier 2"`), `synonyms` (str), `classification_level` (int32), `hgnc_id` (str, only when a `gene_catalog` was supplied).

`hallmark` (bool), `entrez_id` (str), and `genome_location` (str) come from `COSMIC_CGC_LEGACY_FIELDS` and are **legacy-export-only**: a legacy header must carry all three (§ 4), and they appear in the built table, in place of `cosmic_gene_id`, `chromosome`, `genome_start` and `genome_stop`, for a `cosmic-cgc-legacy-v1` build. The current v103+ header carries none of the three, so none of them appear in the schema above, and the builder does not fabricate them.

## 6. hvantk integration points

- Manifest: `hvantk/skills/cosmic_cgc/plugin.yaml` — dataset `cosmic-cgc:submissions`, `artifact_type: AnnotationTable`, `schema_id: cosmic-cgc-v2` (current header) and `schema_ids: [cosmic-cgc-legacy-v1]` (legacy header). No `lifecycle:` or `cli:` block declared.
- Builder: `build_cosmic_cgc_submissions` in `hvantk/skills/cosmic_cgc/builder.py`. Signature `(parsed_input, ctx, **params) -> AnnotationTable`. Recognised `params`: `mutation_context` (`"both"` default, `"somatic"`, or `"germline"`), `min_classification` (one of `COSMIC_CGC_CLASSIFICATION_LEVELS`), `gene_catalog` (a `GeneCatalogStreamer`), `fields`.
- Constants: `COSMIC_CGC_CURRENT_FIELDS`, `COSMIC_CGC_LEGACY_FIELDS`, `COSMIC_CGC_HEADER_GENERATIONS` (schema ID → rename map), `COSMIC_CGC_FIELDS` (both maps merged), `COSMIC_CGC_CLASSIFICATION_LEVELS`, `COSMIC_TISSUE_TYPES`, `COSMIC_MUTATION_CONTEXTS` in `hvantk/skills/cosmic_cgc/shared/constants.py`.
- Drift probe: `fetch_fingerprint` in `hvantk/skills/cosmic_cgc/drift_probe.py` (`PROBE_VERSION = 2`). Probes the public release-notes page only — the gated Census data is never fetched — and anchors on `id="v<N>"` HTML attributes, never on prose (a prior prose-matching approach picked up unrelated `Actionability`-product version mentions). Fails closed if the response redirected (the host 302s the trailing-slash variant to a login form) or if no anchors are found.
- Streamer: `CosmicCGCGeneDiseaseTableStreamer` in `hvantk/skills/cosmic_cgc/streamers.py`, subclassing `GeneDiseaseTableStreamer` (`hvantk/core/streamers/gene_disease_table.py`). Adds `filter_by_mutation_context`, `get_geneset_per_tumour_type`/`_role`/`_tissue`, `mutation_context_summary`, `role_summary`, `classification_summary`, `describe`. Because CGC is one row per gene (not per gene-disease assertion), the MONDO/disease-dependent base methods (`get_genes_by_mondo_id`, `categorize_by_ontology`, `get_genes_by_disease`, `get_genes_by_moi`, `get_geneset_per_disease`, `categorize_by_ontology_summary`) are overridden to raise `NotImplementedError`.
- Tests: `pytest hvantk/skills/cosmic_cgc/tests -m hail` — artifact paths in § 9.

## 7. Workflow steps

1. Obtain the raw COSMIC CGC export manually — requires a COSMIC account/login at `cancer.sanger.ac.uk/census`. No in-repo downloader exists.
2. Build: `hvantk reprocess cosmic-cgc:submissions --raw-dir <dir> --intermediate <dir>/<file>.tsv.gz --skip-parse --skip-download --output <out.ht>`, optionally with `--plugin-arg mutation_context=somatic|germline|both`, `--plugin-arg min_classification=<level>`, `--plugin-arg fields=<comma-list>`, `--plugin-arg hgnc_ht=<hgnc-lookup.ht>` (see step 3).
   `--intermediate <file> --skip-parse` is required: the dataset declares no
   `lifecycle.parse`, so without it `hvantk reprocess` hands the `--raw-dir` directory
   to the builder, which needs the file itself and fails with `IsADirectoryError`.
   `gwas-catalog:associations` documents the same pattern
   (`hvantk/skills/gwas_catalog/SKILL.md` § 6). On the committed synthetic fixture this
   builds a 5-row table keyed by `gene_symbol`.
3. To resolve `hgnc_id`, add `--plugin-arg hgnc_ht=<path-to-hgnc-lookup.ht>` (or
   `hgnc_path=<path>`): `hvantk reprocess` builds an `HGNCGeneCatalogStreamer` from it
   and passes it to the builder as `gene_catalog`, as for `pqtl`. From Python, pass any
   `GeneCatalogStreamer` as `gene_catalog` to `build_cosmic_cgc_submissions(...)`.
4. Sanity-check the output: `classification` values are `"Tier 1"`/`"Tier 2"`, `somatic`/`germline`/`hallmark` are true booleans (not strings), and the multi-value fields are arrays.
5. Run the test suite (§ 9): a loader-registration test, the offline drift-probe tests, and the Hail builder tests, among them a round trip that builds from the committed synthetic fixture and asserts schema and sample rows against committed snapshots.

## 8. Update playbook

Triggered when COSMIC publishes a new release (the `v101`…`v104`… series) or the CGC export's column layout changes.

1. `hvantk drift cosmic-cgc:submissions` compares the live release index against `tests/drift_fingerprint.json`. Regenerate via `hvantk drift --regenerate cosmic-cgc:submissions` once a genuine new release is confirmed.
2. A moved fingerprint means only that a new COSMIC release exists (`source_version` bumped, e.g. `"v104"` → `"v105"`) — it does **not** confirm the Census gene-level contents changed; per `drift_probe.py`, that is undetectable without an authenticated fetch.
3. If COSMIC's CGC export column headers change, the export is a new header generation, which needs three things: its own column map in `shared/constants.py`; an entry for that map in `COSMIC_CGC_HEADER_GENERATIONS` under a new schema ID, which is the detection update (the builder picks generations from that registry); and the new schema ID in `plugin.yaml`'s `schema_ids` (`run_builder_for_spec` rejects a stamp the manifest does not accept). Without them, a header that renames or drops a column raises `ValueError`, and one that only adds columns builds under the old schema ID with a warning naming the extra columns (§ 4). Also update `COSMIC_CGC_CLASSIFICATION_LEVELS` / `COSMIC_MUTATION_CONTEXTS` if the tier or mutation-context vocabulary itself changes.
4. The fixture is synthetic (§ 9) and does not track the real upstream export, so a header or value-convention change is not caught automatically. Verify the builder manually against a freshly (manually) downloaded real export first, then update the synthetic fixture (`tests/testdata/raw/cosmic-cgc/`) and its README to match the new header/conventions, and regenerate snapshots: `pytest hvantk/skills/cosmic_cgc/tests/test_builder.py -m hail --regenerate-snapshots`, inspect the diff by hand, then rerun without the flag to confirm.

## 9. Validation contract

Declared in `plugin.yaml`'s `tests:` block (`hvantk/skills/cosmic_cgc/plugin.yaml`), paths plugin-relative:

- `fixture`: `tests/testdata/raw/cosmic-cgc`
- `schema_snapshot`: `tests/snapshots/schema.json`
- `row_snapshot`: `tests/snapshots/sample_rows.json`
- `drift_fingerprint`: `tests/drift_fingerprint.json`
- `command`: `pytest hvantk/skills/cosmic_cgc/tests -m hail`

`-m hail` selects the five `@pytest.mark.hail`-marked tests in `test_builder.py`, 5 of
the 13 tests collected under `hvantk/skills/cosmic_cgc/tests`: the round trip
(`test_cosmic_cgc_submissions_round_trip`, below), a full legacy header
(`test_legacy_header_fields_are_renamed`), the non-integer coordinate count
(`test_non_integer_coordinates_are_counted`), a current header with `GENOME_START`
renamed raising (`test_partly_renamed_header_raises`), and a gene catalog that resolves
no symbol raising (`test_gene_catalog_resolving_no_symbol_raises`). The other 8 (the
`test_cosmic_cgc.py` registration check plus 7 offline `test_drift_probe.py` cases)
carry no `hail` marker and run under the default (non-`-m hail`) pytest selection
instead.

The COSMIC CGC licence forbids redistributing rows, so `fixture` is a **synthetic,
format-faithful** TSV (`tests/testdata/raw/cosmic-cgc/cosmic-cgc-synthetic.tsv.gz`, with
a README beside it) rather than a truncated real export: the real v103 header,
delimiter, compression, and value conventions (§ 4), but five fabricated rows
(`SYNTHA`-`SYNTHE`) with invented gene symbols and obviously fake free text. A licence
that forbids redistributing *real* rows limits what a fixture may contain, not whether
one exists.

`hvantk/skills/cosmic_cgc/tests/test_builder.py::test_cosmic_cgc_submissions_round_trip`
(`@pytest.mark.hail`) builds from that fixture with a stub `gene_catalog` (maps
`SYNTHA`-`SYNTHD` to invented HGNC ids; `SYNTHE` is deliberately left unmapped) and
asserts: the stamped `schema_id` is `cosmic-cgc-v2`; the built schema matches
`schema_snapshot`; three sampled rows (keyed by `hgnc_id`) match `row_snapshot`; the
unresolved `SYNTHE` row is dropped from the table; `genome_start`/`genome_stop` are
`int32`; and none of the raw `COSMIC_GENE_ID`/`CHROMOSOME`/`GENOME_START`/`GENOME_STOP`
column names leak into the built schema.
`hvantk/skills/cosmic_cgc/tests/test_cosmic_cgc.py::test_cosmic_cgc_submissions_registered`
is loader-only (confirms the manifest resolves `cosmic-cgc:submissions` to
`AnnotationTable`, `schema_id` `cosmic-cgc-v2` and `schema_ids`
`("cosmic-cgc-legacy-v1",)`). `hvantk/skills/cosmic_cgc/tests/test_drift_probe.py`
runs 7 offline (`requests_mock`) tests of the release-index probe: anchor-only matching
excludes unrelated `Actionability`-product prose versions, editorial prose edits don't
move the signal, forward-looking prose can't bump `source_version`, numeric sort,
new-release detection, and fail-closed on a login redirect or missing anchors.
`tests/drift_fingerprint.json` **is** committed (`source_version: "v104"`,
`headers.cosmic-release-index: ["v101", "v102", "v103", "v104"]`) — the release-index
signal and the fixture/snapshot round-trip are both real, running checks.
