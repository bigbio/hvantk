---
name: hvantk:resource-pqtl
description: Protein quantitative trait loci (pQTL) association metrics from Fang et al. (2025) GTEx TMT proteomics, built into a Hail Table keyed by (locus, alleles, gene_id).
status: provisional
backend: hail
domain: proteomics
---

# pqtl

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

Provisional. Covers the `pqtl:metrics` builder, which converts Fang et al. (2025) cis-pQTL allpairs summary statistics into a Hail Table keyed by `(locus, alleles, gene_id)`, for use as the protein-level arm of the QTL cascade join against eQTL data.

Out of scope for this skill (per `_conventions` § 11):
- Acquiring the raw supplementary data. There is no `lifecycle.download`/`lifecycle.parse` block in `plugin.yaml` and no downloader is implemented — the file must be obtained manually as preprint supplementary material.
- The HGNC lookup table the gene-symbol → Ensembl-ID mapping depends on. See `hvantk/skills/hgnc/SKILL.md`.
- QTL cascade join logic itself (coloc, cascade traversal). Lives in `hvantk/algorithms/qtlcascade/`.

## 2. Source identity

Fang et al. (2025), "Regulation of protein abundance in normal human tissues", medRxiv, doi:`10.1101/2025.01.10.25320181` — 10,841 proteins across >700 GTEx samples in five tissues, quantified by TMT mass spectrometry; the cis-pQTL allpairs files are what the builder consumes (`hvantk/skills/pqtl/drift_probe.py` module docstring). `source.catalog_ref: pqtl` is declared in `plugin.yaml`; this skill does not restate catalog content (per `_conventions` § 2) — query it via `hvantk catalog show pqtl`.

The five profiled tissues are `Colon`, `Heart`, `Liver`, `Lung`, `Thyroid` (`FANG_TISSUES`, `hvantk/skills/pqtl/shared/constants.py`), `FANG_TISSUE_EQTL_MAPPING` in the same file records how these correspond to GTEx eQTL sub-tissue names (e.g. `Colon` → `Colon_Transverse`, `Colon_Sigmoid`; `Heart` → `Heart_Left_Ventricle`, `Heart_Atrial_Appendage`), but **nothing applies it** — it has zero call sites outside its own definition. The builder writes the raw Fang tissue name through verbatim, and `build_cascade` (`hvantk/algorithms/qtlcascade/cascade.py`) filters both the eQTL and pQTL tables by exact `==` against one `tissue` string. So `tissue="Colon"` will **not** match GTEx labels like `Colon_Transverse`; a caller wanting that alignment must apply the mapping itself.

There is no direct data URL: the summary statistics are supplementary material to a preprint, and acquisition is manual (`drift_probe.py` docstring). Because this is publication-only material, the drift probe treats the publication's own version metadata as the upstream signal (see § 3).

Only `source="gtex_fang"` is currently implemented; `PQTL_SOURCES = ("gtex_fang",)` in `shared/constants.py` is the full supported set, and `build_pqtl_metrics` raises `NotImplementedError` for any other value even if additional names were ever added to `PQTL_SOURCES`.

## 3. Backend choice + reasoning

`hail`, per `plugin.yaml` (`backend: hail`). The builder produces a Hail Table keyed by `(locus, alleles, gene_id)` — the eQTL/pQTL domain convention (`_conventions` § 3) — so it can be joined directly against eQTL Hail Tables in the QTL cascade pipeline without an intermediate conversion step.

## 4. Raw format & gotchas

- Raw format: Fang allpairs files, space-delimited gzip, columns `gene_name SNP CHR BP A1 NMISS BETA STAT P` (`_import_gtex_fang` docstring, `hvantk/skills/pqtl/builder.py`). `scan_tissue_files` (`hvantk/core/utils/qtl_helpers.py`) looks for `.txt.gz`/`.tsv.gz` files; tissue name is inferred from the filename prefix before the first dot (e.g. `Colon.v8.allpairs.txt.gz` → `Colon`).
- Imported with `hl.import_table(delimiter=" ", force=True, types={"BETA": hl.tfloat64, "STAT": hl.tfloat64, "P": hl.tfloat64})` (`_import_gtex_fang`). Rows where `STAT == 0.0` are dropped before further processing, since SE cannot be derived from them.
- Variant IDs are GTEx-format (`chr1_1000050_C_T_b38`), parsed by the **shared** `parse_gtex_variant_id` in `hvantk/core/utils/qtl_helpers.py` (also used by the eQTL builder) — contig prefix normalized to match `reference_genome` (default `GRCh38`), build suffix discarded.
- SE is **not** present in the Fang allpairs files. It is derived in the builder as `se = |BETA / STAT|` (`build_pqtl_metrics`, `builder.py`) — this is why `STAT == 0` rows must be filtered first (division by zero).
- Gene-symbol → Ensembl-ID mapping is **required by default**. `build_pqtl_metrics` raises `ValueError` if `gene_catalog` is `None` and `no_gene_map` is not `True`, because the cascade join keys on `(locus, alleles, gene_id)` with Ensembl IDs on the eQTL side — raw gene symbols would produce zero matches. The mapping is performed via `GeneCatalogStreamer.map_ids(list(symbols), source_type="gene_symbol", target_type="ensembl_gene_id")` (`hvantk/core/streamers/gene_catalog.py`), collecting the full distinct symbol set with `hl.agg.collect_as_set` and broadcasting the result back as `hl.literal(ensembl_mapping)`. Symbols with no mapping fall back to the raw symbol itself (`hl.or_else(mapping_literal.get(ht.gene_symbol), ht.gene_symbol)`).
- Passing `no_gene_map=True` opts out of mapping and keys `gene_id` as the raw gene symbol — the resulting table will **not** join with eQTL tables in cascade analysis (explicit `logger.warning` in `build_pqtl_metrics`).
- Optional `p_threshold` filters rows to `p_value <= p_threshold` after gene mapping. Optional `fields` selects a column subset via `ht.select(*fields)` at the end.
- Every row gets `source=source` (the string passed in, e.g. `"gtex_fang"`) and `is_cis=True` annotated unconditionally (`build_pqtl_metrics`) — there is no trans-pQTL support yet.
- Final key is `("locus", "alleles", "gene_id")`, set via `ht.key_by(...)` after the `gene_id` annotation.

## 5. Output contract

Hail Table keyed by `(locus, alleles, gene_id)`. Builder returns an `AnnotationTable` (`hvantk.core.models.AnnotationTable`) via `AnnotationTable.from_hail(ht, provenance=ctx.provenance(schema_id="pqtl-v1"))` (`build_pqtl_metrics`, `builder.py`).

Fields confirmed present from the builder's transform chain (§ 4): `locus`, `alleles`, `gene_id` (key fields); `gene_symbol` (retained pre-mapping value); `beta`, `se`, `p_value`, `tissue`, `source`, `is_cis` (value fields). `variant_id` and `stat` are explicitly dropped (`ht.drop("stat", "variant_id")`) before the key is set. This skill does not claim a fixed field list beyond what `builder.py` shows, since no schema snapshot has ever been generated (see § 9) — treat the above as builder-verified, not as the full committed schema.

## 6. hvantk integration points

- Plugin manifest: `hvantk/skills/pqtl/plugin.yaml` (`api_version: 2`, dataset `metrics`, `artifact_type: AnnotationTable`, `schema_id: pqtl-v1`). The loader resolves it via `get_registry().get_dataset("pqtl:metrics")`.
- Builder: `build_pqtl_metrics` in `hvantk/skills/pqtl/builder.py`. Signature: `(parsed_input, ctx, *, reference_genome="GRCh38", source="gtex_fang", tissue=None, gene_catalog=None, no_gene_map=False, p_threshold=None, fields=None) -> AnnotationTable`.
- Shared helpers: `parse_gtex_variant_id`, `scan_tissue_files` in `hvantk/core/utils/qtl_helpers.py` (shared with the eQTL builder — do not duplicate this logic in a pqtl-local module).
- Constants: `PQTL_SOURCES`, `FANG_TISSUES`, `FANG_TISSUE_EQTL_MAPPING` in `hvantk/skills/pqtl/shared/constants.py`.
- Gene mapping dependency: `GeneCatalogStreamer` (base class, `hvantk/core/streamers/gene_catalog.py`), `map_ids(...)` method — the `gene_catalog` param is normally an `HGNCGeneCatalogStreamer` instance built from `hvantk reprocess hgnc:lookup` output (see `hvantk/skills/hgnc/SKILL.md`). The `hvantk reprocess` CLI resolves `--plugin-arg hgnc_ht=<path>` (or `hgnc_path`) into `gene_catalog=HGNCGeneCatalogStreamer.from_path(<path>)` in `hvantk/tools/plugins/reprocess_cli.py`.
- Drift probe: `fetch_fingerprint` in `hvantk/skills/pqtl/drift_probe.py`, wired via `plugin.yaml`'s `drift_probe:` block; drives `hvantk drift pqtl:metrics`.
- No CLI `lifecycle:`/`cli:` blocks are declared in `plugin.yaml` — there is no built-in downloader.
- Existing tests: `hvantk/skills/pqtl/tests/test_pqtl.py` (registration-only + skipped round-trip), `test_drift_probe.py` (probe behavior, offline).

## 7. Workflow steps

When invoked to build or refresh the pQTL metrics table:

1. **Confirm the raw Fang allpairs file(s) are present** as `.txt.gz`/`.tsv.gz` under `--raw-dir` (§ 4). Acquisition is manual — this workflow does not download them (§ 1).
2. **Build the HGNC lookup table first**, unless intentionally opting out with `no_gene_map=true`: `hvantk reprocess hgnc:lookup --raw-dir <dir> --output <hgnc.ht>`.
3. **Build** via `hvantk reprocess pqtl:metrics --raw-dir <dir-with-fang-allpairs> --output <out.ht> --plugin-arg hgnc_ht=<hgnc-lookup.ht>`. The CLI translates `hgnc_ht` (or `hgnc_path`) into the builder's `gene_catalog` kwarg by constructing `HGNCGeneCatalogStreamer.from_path(...)` (`hvantk/tools/plugins/reprocess_cli.py`) — do not pass `--plugin-arg gene_catalog=...` directly. To opt out: add `--plugin-arg no_gene_map=true` instead.
4. **Restrict to one tissue** with `--plugin-arg tissue=<TissueName>` (one of `FANG_TISSUES`) if a cascade run only needs a matched subset.
5. **Sanity-check the output.** Confirm the table is keyed by `(locus, alleles, gene_id)`, that `gene_id` values look like Ensembl IDs (`ENSG...`) rather than raw symbols (unless `no_gene_map=true` was intentional), and that `se` has no `inf`/`NaN` values (would indicate a `STAT == 0` row slipped through).
6. **Do not** pass a `source` other than `"gtex_fang"` — `build_pqtl_metrics` raises `NotImplementedError` for anything else.

## 8. Update playbook

Triggered when the drift probe reports a new preprint version or journal publication for doi:`10.1101/2025.01.10.25320181` (§ 2, § 3).

1. Run `hvantk drift pqtl:metrics` (or wait for the fortnightly automated drift workflow) to check the medRxiv API's `version`/`date`/`published` fields against the committed `tests/drift_fingerprint.json`.
2. Re-read the drift probe's module docstring (`drift_probe.py`) before acting: it detects a new preprint version or journal publication, but explicitly **cannot** detect an in-place replacement of the supplementary file under an unchanged version — that gap is a property of how the source is published, not a probe bug.
3. If `published` moves from `"NA"` to a journal DOI, or `version` increments, manually re-download the supplementary allpairs files and diff the column header against `_import_gtex_fang`'s expected columns (`gene_name SNP CHR BP A1 NMISS BETA STAT P`) before rebuilding.
4. Regenerate the fingerprint with `hvantk drift --regenerate pqtl:metrics` once the change is validated; do not regenerate in the same PR as a behavioral change (`_conventions` § 12).

## 9. Validation contract

Declared in `plugin.yaml`'s `tests:` block (all paths plugin-relative under `hvantk/skills/pqtl/`):

- `fixture`: `tests/testdata/raw/pqtl`
- `schema_snapshot`: `tests/snapshots/schema.json`
- `row_snapshot`: `tests/snapshots/sample_rows.json`
- `drift_fingerprint`: `tests/drift_fingerprint.json`
- `command`: `pytest hvantk/skills/pqtl/tests`

**Snapshot status:** `pqtl:metrics` is on the `KNOWN_INCOMPLETE` ledger in `hvantk/tests/test_plugin_contract_artifacts.py`, missing `fixture`, `schema_snapshot`, and `row_snapshot`. Per that file's comment block, pQTL is cause (1) of the two the ledger now records: it is "publication-only supplementary data" with "no static upstream artifact" that "cannot be snapshotted" — the Fang et al. allpairs files are preprint supplementary material with no stable, redistributable download URL to derive a committable fixture from, so no fixture and no schema/row snapshot can be generated. `tests/test_pqtl.py::test_pqtl_metrics_round_trip` is `@pytest.mark.skip`'d for this reason ("No fixture available for pqtl; manual smoke-test only"). Only `drift_fingerprint` is populated (from a live probe run against the medRxiv API) — see § 3 and § 8.
