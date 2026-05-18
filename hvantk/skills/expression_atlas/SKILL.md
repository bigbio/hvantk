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

- **Status:** provisional. Builder, downloader, dataset class, and downloader tests are in place under the plugin folder; round-trip builder snapshots are NOT yet seeded (no fixture, no `schema.json`, no `sample_rows.json`) -- writing the first round-trip test is follow-up work tracked alongside this skill.
- **In scope:** any single Expression Atlas baseline bulk-RNA-seq experiment with a gene-centric TPM matrix (genes x samples) and a paired SDRF metadata file, converted to an AnnData object keyed `samples x genes`.
- **Out of scope:** scRNA-seq cell-level matrices (use the UCSC Cell Browser plugin); differential expression matrices; cross-accession multi-experiment merging; gene-symbol / accession normalization (downstream).

## 2. Source identity

- **Provider:** EBI Expression Atlas (<https://www.ebi.ac.uk/gxa>).
- **Catalog entry:** TODO. `hvantk/resources/catalog.yaml` does not currently expose a top-level `expression-atlas` entry; the per-accession dataset list lives in `hvantk/resources/registry/transcriptomics/datasets.json` (filterable via `data_source == "Expression_Atlas"`). When a dedicated catalog entry lands, point the plugin manifest's `source.catalog_ref` at it and remove this TODO.
- Per-accession download URLs derive from `EXPRESSION_ATLAS_BASE_URL` in `hvantk/core/constants.py` and the FTP prefix `/pub/databases/microarray/data/atlas/experiments/<accession>/`.

## 3. Backend choice + reasoning

**`backend: anndata`, `domain: transcriptomics`.** Per `_conventions` § 3 "Expression (dense, Hail-friendly) → MatrixTable" vs "Expression (sparse, single-cell) → anndata h5ad" — Expression Atlas baseline bulk-RNA-seq sits between the two: matrices are dense floats but typically small enough (tens of thousands of genes x hundreds of samples) that the scanpy / AnnData ecosystem is the natural fit. Downstream consumers (cross-tissue baseline expression queries, tissue-specificity scoring via `tspex`) operate on AnnData. Hail offers no advantage at this scale.

## 4. Raw format & gotchas

The exact raw-format gotchas live where the parser does:

- Expression TSV parser + AnnData assembly: `hvantk/skills/expression_atlas/shared/expression_atlas.py` (`create_anndata_from_expression_atlas`).
- SDRF long → wide reshape: same file (`_import_sdrf`, `_reshape_sdrf_long_to_wide_format`, `convert_sdrf_to_dataframe`).
- High-level builder narrative: docstring on `build_expression_atlas_ad` in `hvantk/skills/expression_atlas/builder.py`.

Stable notes:

- Expression matrix is gene-centric (rows = genes, columns = samples). The first two columns are `Gene ID` and `Gene Name`; everything after that is one float column per sample. The builder transposes to AnnData's samples-as-obs convention.
- SDRF is tab-separated, no header, and condensed-long: each row is `(accession, unused, sample_id, column_type, column_name, column_value)`. The `unused` column is dropped; `column_type` is either `characteristic` or `factor`; duplicates on `(sample_id, column_name)` keep the **last** value (see `_reshape_sdrf_long_to_wide_format`).
- Column names from SDRF are normalized: spaces → underscores, parentheses stripped (`organism_part_(group)` → `organism_part_group`). Downstream `obs` column names follow this rule.
- TODO: enumerate per-accession quirks (e.g. mixed-type factor columns, missing SDRF rows, multi-pipeline TPM variants) once the round-trip fixture is seeded. Add specific cases here as they are encountered.

## 5. Output contract

- **Object:** `anndata.AnnData` optionally saved to `<output_path>.h5ad`.
- **Shape:** `obs = samples`, `var = genes`. `X` is `float32` (genes-x-samples after transpose).
- **`obs`:** indexed by `sample_id`; columns are SDRF characteristics / factors after the long → wide reshape (e.g. `organism`, `tissue`, `cell_type`, ...).
- **`var`:** indexed by `gene_id`. Includes a `Gene Name` column when the source TSV had one.
- **`uns["hvantk_metadata"]`:** provenance dict from `hvantk.core.anndata_utils.build_anndata_metadata("ExpressionAtlas", expression_matrix_path)`.
- **`uns["column_summary"]`:** per-`obs`-column summary annotated by `annotate_column_summary_ad`.

## 6. hvantk integration points

- **Builder:** `build_expression_atlas_ad` in `hvantk/skills/expression_atlas/builder.py`.
- **SDRF / matrix helpers:** `hvantk/skills/expression_atlas/shared/expression_atlas.py`.
- **Dataset / collection classes:** `ExpressionAtlasDataset`, `ExpressionAtlasDatasetCollection` in `hvantk/skills/expression_atlas/shared/datasets.py`.
- **Downloader CLI:** `download_experiments` in `hvantk/skills/expression_atlas/cli.py` (registered as `hvantk expression-atlas-download` and also re-bound under `hvantk download expression-atlas`).
- **Lifecycle entry point:** `download_dataset` in `hvantk/skills/expression_atlas/cli.py`.
- **Build CLI:** `hvantk mkmatrix expression-atlas` in `hvantk/tools/build/make_matrix_cli.py` (delegates to the plugin builder).
- **Plugin manifest:** `hvantk/skills/expression_atlas/plugin.yaml` (drives loader registration; compound dataset key `expression-atlas:dataset`).
- **Tests:** `hvantk/skills/expression_atlas/tests/` (downloader unit + drift-probe sanity present; builder round-trip TODO).

Read the existing files at these paths as ground truth for shape. This skill does not restate code.

## 7. Workflow steps

When invoked to build or update a single Expression Atlas experiment:

1. **Resolve raw paths.** Either download via
   `hvantk expression-atlas-download --accession <E-XXXX-N> --download_path /tmp/atlas`,
   or via the recipe system: `hvantk reprocess expression-atlas:dataset` (lifecycle download → builder).
2. **Build:** `hvantk mkmatrix expression-atlas -e <expression>.tsv -s <accession>.sdrf.tsv -o <out>.h5ad`.
   - The builder parses the SDRF, transposes the expression matrix, attaches per-sample metadata into `obs`, annotates provenance, and writes `.h5ad`.
3. **Validate:** TODO — once the round-trip fixture is seeded, run `pytest hvantk/skills/expression_atlas/tests`. Until then, the offline downloader unit tests + drift-probe placeholder test are what guard this plugin.

## 8. Update playbook

TODO. This section will be fleshed out once per-accession drift detection lands (see § 2 catalog note and the drift-probe placeholder in `hvantk/skills/expression_atlas/drift_probe.py`). Expected shape:

1. For each tracked accession in `registry/transcriptomics/datasets.json` (filter `data_source == "Expression_Atlas"`), re-run the per-accession HEAD probe; flag accessions whose `Last-Modified` or `Content-Length` changed.
2. Re-download flagged accessions, rebuild via `hvantk mkmatrix expression-atlas`, and diff the new AnnData against the snapshotted shape / `obs` columns.
3. If the SDRF column set changed, document the new factor in § 4.

## 9. Validation contract

Per `_conventions` § 9:

- **fixture:** `hvantk/skills/expression_atlas/tests/testdata/raw/expression-atlas/` (directory present, not yet seeded). The first round-trip test will populate this with a small TSV + condensed-SDRF pair.
- **schema_snapshot:** `hvantk/skills/expression_atlas/tests/snapshots/schema.json` (TODO — created on first `--regenerate-snapshots` run).
- **row_snapshot:** `hvantk/skills/expression_atlas/tests/snapshots/sample_rows.json` (TODO — same).
- **test_command:** `pytest hvantk/skills/expression_atlas/tests`.

The plugin manifest already declares these paths so the loader contract holds. The downloader unit tests + drift-probe placeholder test pass today; the builder round-trip is the gap to close in a follow-up PR.

> **Snapshot status:** schema.json and sample_rows.json have NOT yet been seeded
> for this plugin. On first round-trip run in a hail-enabled environment, use
> `pytest hvantk/skills/expression_atlas/tests/test_builder.py --regenerate-snapshots`
> to bootstrap them, then commit. Until seeded, the round-trip test cannot verify
> output against a fixed schema.
