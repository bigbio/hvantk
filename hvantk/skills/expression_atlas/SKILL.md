---
name: hvantk:resource-expression-atlas
description: Build an AnnData object from an EBI Expression Atlas baseline bulk-RNA-seq experiment (TPM expression matrix + SDRF sample metadata).
status: provisional
backend: anndata
domain: transcriptomics
---

# EBI Expression Atlas resource skill

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes every convention there.

## 1. Status & scope

- **Status:** provisional. Builder, downloader, dataset class and downloader tests are in place, and the round-trip contract is **seeded**: `tests/testdata/raw/expression-atlas/` holds a truncated fixture, `tests/snapshots/` holds `schema.json` + `sample_rows.json`, and `tests/test_builder.py` asserts against both (§ 9). The fixture reproduces the **real** upstream header, including the third `GeneID` transcript column that used to break the builder (#342, fixed).
- **In scope:** any single Expression Atlas baseline bulk-RNA-seq experiment with a gene-centric TPM matrix (genes x samples) and a paired SDRF metadata file, converted to an AnnData object keyed `samples x genes`.
- **Out of scope:** scRNA-seq cell-level matrices (use the UCSC Cell Browser plugin); differential expression matrices; cross-accession multi-experiment merging; gene-symbol / accession normalization (downstream).

## 2. Source identity

- **Provider:** EBI Expression Atlas (<https://www.ebi.ac.uk/gxa>).
- **Catalog entry:** the per-accession dataset list lives in `hvantk/skills/expression_atlas/catalog/datasets.json` (filterable via `data_source == "Expression_Atlas"`). Browse with `hvantk catalog list --data-source Expression_Atlas` or `hvantk catalog show E-GTEX-8`. A top-level `expression-atlas` provider-level catalog entry remains a TODO; when added, point the plugin manifest's `source.catalog_ref` at it and remove this note.
- Per-accession download URLs derive from `EXPRESSION_ATLAS_BASE_URL` in `hvantk/skills/expression_atlas/shared/constants.py` and the FTP prefix `/pub/databases/microarray/data/atlas/experiments/<accession>/`.

## 3. Backend choice + reasoning

**`backend: anndata`, `domain: transcriptomics`.** Per `_conventions` § 3 "Expression (dense, Hail-friendly) → MatrixTable" vs "Expression (sparse, single-cell) → anndata h5ad" — Expression Atlas baseline bulk-RNA-seq sits between the two: matrices are dense floats but typically small enough (tens of thousands of genes x hundreds of samples) that the scanpy / AnnData ecosystem is the natural fit. Downstream consumers (cross-tissue baseline expression queries, tissue-specificity scoring via `tspex`) operate on AnnData. Hail offers no advantage at this scale.

## 4. Raw format & gotchas

The exact raw-format gotchas live where the parser does:

- Expression TSV parser + AnnData assembly: `hvantk/skills/expression_atlas/shared/expression_atlas.py` (`create_anndata_from_expression_atlas`).
- SDRF long → wide reshape: same file (`_import_sdrf`, `_reshape_sdrf_long_to_wide_format`, `convert_sdrf_to_dataframe`).
- High-level builder narrative: docstring on `build_expression_atlas` in `hvantk/skills/expression_atlas/builder.py`.

Stable notes:

- Expression matrix is gene-centric (rows = genes, columns = samples), transposed by the builder to AnnData's samples-as-obs convention.
- **There are THREE leading metadata columns, not two.** The real header is `Gene ID`, `Gene Name`, `GeneID`, then one float column per sample. `Gene ID` is the gene; **`GeneID` (no space) is the per-row transcript id** — the two differ by a single space, so any matching that is loose about whitespace will drop the wrong one. Treating the third as a sample is what made every build from an unmodified download die with `ValueError: could not convert string to float: 'ENSMUST...'` (#342).
- `create_anndata_from_expression_atlas` therefore sets aside, as annotations, any leading column that holds values of which **none** are numeric, and preserves them in `var` rather than dropping them. Name them in `extra_annotation_columns` to skip the inference. The rule keys on content, not position: an all-missing sample column coerces to NaN (numeric) and correctly stays a sample.
- Because a `*-transcripts-tpms.tsv` export is transcript-level, **one gene id spans several rows**, so `var_names` are not unique unless rows are de-duplicated (the committed fixture keeps the first row per gene). `var["GeneID"]` is what disambiguates them.
- SDRF is tab-separated, no header, and condensed-long: each row is `(accession, unused, sample_id, column_type, column_name, column_value)`. The `unused` column is dropped; `column_type` is either `characteristic` or `factor`; duplicates on `(sample_id, column_name)` keep the **last** value (see `_reshape_sdrf_long_to_wide_format`).
- Column names from SDRF are normalized: spaces → underscores, parentheses stripped (`organism_part_(group)` → `organism_part_group`). Downstream `obs` column names follow this rule.
- TODO: enumerate remaining per-accession quirks (mixed-type factor columns, missing SDRF rows, multi-pipeline TPM variants) as they are encountered.

## 5. Output contract

- **Object:** `anndata.AnnData` optionally saved to `<output_path>.h5ad`.
- **Shape:** `obs = samples`, `var = genes`. `X` is `float32` (genes-x-samples after transpose).
- **`obs`:** indexed by `sample_id`; columns are SDRF characteristics / factors after the long → wide reshape (e.g. `organism`, `tissue`, `cell_type`, ...).
- **`var`:** indexed by `gene_id`. Includes a `Gene Name` column when the source TSV had one.
- **`uns["column_summary"]`:** per-`obs`-column summary annotated by `annotate_column_summary_ad`.
- **Provenance:** stamped on the returned `ExpressionMatrix` via `ctx.provenance(schema_id="expression-atlas-dataset-v1")`; persisted by the platform as a sidecar `.provenance.json`.

## 6. hvantk integration points

- **Builder:** `build_expression_atlas` (signature `(parsed_input, ctx, **params) -> ExpressionMatrix`) in `hvantk/skills/expression_atlas/builder.py`.
- **SDRF / matrix helpers:** `hvantk/skills/expression_atlas/shared/expression_atlas.py`.
- **Dataset / collection classes:** `ExpressionAtlasDataset`, `ExpressionAtlasDatasetCollection` in `hvantk/skills/expression_atlas/shared/datasets.py`.
- **Downloader CLI:** `download_experiments` in `hvantk/skills/expression_atlas/cli.py` (registered as `hvantk expression-atlas-download` and also re-bound under `hvantk download expression-atlas`).
- **Lifecycle entry point:** `download_dataset` in `hvantk/skills/expression_atlas/cli.py`.
- **Build CLI:** `hvantk reprocess expression-atlas:dataset --raw-dir <dir> --output <path>.h5ad` (delegates to the plugin builder; pass builder kwargs via `--plugin-arg key=value`).
- **Plugin manifest:** `hvantk/skills/expression_atlas/plugin.yaml` (drives loader registration; compound dataset key `expression-atlas:dataset`).
- **Tests:** `hvantk/skills/expression_atlas/tests/` (downloader unit + drift-probe sanity present; builder round-trip TODO).

Read the existing files at these paths as ground truth for shape. This skill does not restate code.

## 7. Workflow steps

When invoked to build or update a single Expression Atlas experiment:

1. **Resolve raw paths.** Either download via
   `hvantk expression-atlas-download --accession <E-XXXX-N> --download_path /tmp/atlas`,
   or via the recipe system: `hvantk reprocess expression-atlas:dataset` (lifecycle download → builder).
2. **Build:** `hvantk reprocess expression-atlas:dataset --raw-dir <dir> --output <out>.h5ad`.
   - The builder parses the SDRF, transposes the expression matrix, attaches per-sample metadata into `obs`, annotates provenance, and writes `.h5ad`.
3. **Validate:** run `pytest hvantk/skills/expression_atlas/tests`. That covers the seeded round-trip (`test_builder.py`, schema + sample rows + `n_obs`/`n_vars`) alongside the offline downloader and drift-probe tests. Regenerate snapshots only for an intentional change: `pytest hvantk/skills/expression_atlas/tests/test_builder.py --regenerate-snapshots`. Building from a real, untruncated download works (#342 fixed); verified against the full 320-column `E-MTAB-6798` file.

## 8. Update playbook

TODO. This section will be fleshed out once per-accession drift detection lands (see § 2 catalog note and the drift-probe placeholder in `hvantk/skills/expression_atlas/drift_probe.py`). Expected shape:

1. For each tracked accession in `hvantk/skills/expression_atlas/catalog/datasets.json` (filter `data_source == "Expression_Atlas"`), re-run the per-accession HEAD probe; flag accessions whose `Last-Modified` or `Content-Length` changed.
2. Re-download flagged accessions, rebuild via `hvantk reprocess expression-atlas:dataset`, and diff the new AnnData against the committed snapshots in `tests/snapshots/` (§ 9).
3. If the SDRF column set changed, document the new factor in § 4.

## 9. Validation contract

Per `_conventions` § 9:

- **fixture:** `hvantk/skills/expression_atlas/tests/testdata/raw/expression-atlas/` — seeded. `E-MTAB-6798-transcripts-tpms.tsv` (20 genes x 4 samples, ~0.8 KB) + `E-MTAB-6798.condensed-sdrf.tsv` (the same 4 sample IDs, ~3.3 KB), derived by truncation from the real upstream files under `hvantk/tests/testdata/raw/expression_atlas/`. Recipe recorded at the top of `tests/test_builder.py`.
- **schema_snapshot:** `hvantk/skills/expression_atlas/tests/snapshots/schema.json` — seeded (4 obs x 20 vars).
- **row_snapshot:** `hvantk/skills/expression_atlas/tests/snapshots/sample_rows.json` — seeded.
- **test_command:** `pytest hvantk/skills/expression_atlas/tests`.

The plugin manifest already declares these paths so the loader contract holds. `tests/test_builder.py` now exercises `build_expression_atlas` end-to-end against the committed fixture and asserts both snapshots plus `n_obs`/`n_vars`; regenerate via `pytest hvantk/skills/expression_atlas/tests/test_builder.py --regenerate-snapshots`. No Hail is required — the artifact is AnnData-backed.

> **The fixture is production-shaped.** It keeps the real three-column header
> (`Gene ID`, `Gene Name`, `GeneID`), so the round-trip test exercises the layout an
> actual download has — not a sanitised one. It previously dropped the `GeneID`
> transcript column to work around #342; that bug is fixed and the workaround is gone.
> `tests/test_builder.py` additionally pins the two behaviours directly: the transcript
> column must reach `var` and not become a third sample, and an all-missing sample
> column must stay a sample.
