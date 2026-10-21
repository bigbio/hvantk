---
name: hvantk:resource-expression-atlas
description: Build an AnnData object from an EBI Expression Atlas baseline bulk-RNA-seq experiment (TPM expression matrix + SDRF sample metadata).
status: provisional
backend: anndata
domain: transcriptomics
---

# EBI Expression Atlas resource skill

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

- **Status:** provisional. Builder, downloader, dataset class and downloader tests are in place, and the round-trip contract is **seeded**: `tests/testdata/raw/expression-atlas/` holds a truncated fixture, `tests/snapshots/` holds `schema.json` + `sample_rows.json`, and `tests/test_builder.py` asserts against both (§ 9). The fixture reproduces the **real** upstream header, including the third `GeneID` transcript column that used to break the builder (#342, fixed).
- **In scope:** any single Expression Atlas baseline bulk-RNA-seq experiment with a gene-centric TPM matrix (genes x samples) and a paired SDRF metadata file, converted to an AnnData object keyed `samples x genes` (`samples x transcripts` for a transcript-level export — see § 4/§ 5).
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
- `create_anndata_from_expression_atlas` therefore sets aside, as annotations, leading columns that hold values of which **none** are numeric, and preserves them in `var` rather than dropping them. Name them in `extra_annotation_columns` to skip the inference.
- The inference is bounded on three sides, because misfiring is worse than the original crash — a real sample reclassified as an annotation disappears from the matrix and the build *succeeds*:
  1. a column named in the SDRF (`metadata_df.index`) is a sample whatever it contains, so a bad sentinel (`-`, `n/a `, a comma decimal) reaches the float cast and raises;
  2. only the **leading run** of columns is considered — everything from the first sample onward is data, so a non-numeric column mid-matrix is malformed, not metadata;
  3. conversely, a **numeric** leading column the SDRF does *not* declare raises rather than being absorbed as a sample. A numeric annotation column (`Entrez`, gene length, an analytics `p-value`) is indistinguishable by content from a sample the SDRF omits, and the old code guessed "sample" — putting NCBI gene ids into the matrix as expression values, detectable only as an all-NaN `obs` row.
  Within those bounds content decides, since position alone is not enough (the transcript column is not always third). Note an all-missing column reads as all-NaN, which `to_numeric` also reports as **non**-numeric; it survives as a sample only because of the `values.notna().any()` guard, not because NaN counts as numeric. That guard is load-bearing.
- Because a `*-transcripts-tpms.tsv` export is transcript-level, **one gene id spans several rows** — `Gene ID` repeats once per transcript of that gene. `var_names` therefore cannot be the gene id: anndata accepts duplicate `var_names` and only warns, so a naive build round-trips an `.h5ad` whose index is non-unique and silently returns multiple columns from `adata[:, gene_id]` downstream (#349). `create_anndata_from_expression_atlas` detects the transcript id column (`transcript_id_column`, default `"GeneID"`) among the inferred annotation columns and, when present, indexes `var` by it instead — unique per row by construction — demoting `gene_id_column` to a regular `var` column (`"Gene ID"`) so the gene id is kept, just no longer as the index. A gene-level export (no transcript column detected) is unaffected and still indexes `var` by the gene id.
- **`hvantk reprocess` scopes the expression matrix to the SDRF's accession, and picks the TPM file whose columns are the SDRF's runs.** `_resolve_expression_atlas_inputs` (`builder.py`) reads the accession out of the sole `*.condensed-sdrf.tsv`'s own filename (`<accession>.condensed-sdrf.tsv`) and only considers TPM files named for that same accession -- SDRF and TPM candidates used to be counted independently, so a `--raw-dir` holding one accession's SDRF and a *different* accession's TPM file (a reused directory, two downloads merged by hand) silently passed and combined mismatched expression + metadata. `--plugin-arg accession=...` with no `config_path` filter fetches *every* file in the experiment's FTP directory (`ftp.nlst()`, `hvantk/skills/expression_atlas/cli.py`), so `<accession>-tpms.tsv` (gene-level) and `<accession>-transcripts-tpms.tsv` can both be present, and the name does not say which can be built: E-MTAB-6798's gene-level file has one column per assay group (`g1` … `g94`, each cell a string of quartiles), its transcript-level file one column per run. So the resolver reads each candidate's header line and keeps the files that have at least one SDRF run id among their columns. Exactly one is used; of several, the gene-level one, since a per-run gene file needs no aggregation. It raises when the SDRF's accession matches no TPM file, when no candidate has a run column, or when the choice is still ambiguous; the file actually chosen is logged at INFO.
- **`level="gene"` collapses a transcript-level export back to one row per gene.** `hvantk expression markers`/`summarize` and the PTM constraint expression loader read `var_names` as the gene id with no override; against a transcript-keyed `var` they would silently emit per-transcript results labelled as per-gene (#383). This only matters when the resolved input is transcript-level (see above) — pass `level=gene` (e.g. `--plugin-arg level=gene` on `hvantk reprocess`) to sum each sample's transcript-level values onto their shared `Gene ID` and re-key `var` by the gene id (index name `gene_id`, matching a gene-level build; `Gene Name` carried along). Summation is correct for the abundance units Expression Atlas ships (TPM, FPKM, counts) — all linear, so a gene's total is the sum of its transcripts' — using `min_count=1` so a group of all-missing inputs sums to `NaN`, not `0.0`; a transcript with no `Gene ID` raises rather than being dropped by `groupby`. Against an already gene-level export, `level=gene` is a harmless no-op. The default, `level=transcript`, is unchanged. Any other value raises. Those same two consumers refuse a transcript-keyed AnnData outright (`hvantk.core.models.anndata_utils.require_gene_keyed`, checked via `var.index.name == "transcript_id"`) rather than silently mislabelling transcripts as genes, and the error names `--plugin-arg level=gene` as the fix. `hvantk annotate prepare`'s matrix source type refuses one too (defense in depth), but `level=gene` does not make an Expression Atlas build usable there: it resolves `var_names` as HGNC gene symbols, not Ensembl gene ids, so it needs a symbol-keyed reduction regardless of `level`.
- SDRF is tab-separated, no header, and condensed-long: each row is `(accession, unused, sample_id, column_type, column_name, column_value)`. The `unused` column is dropped; `column_type` is either `characteristic` or `factor`; duplicates on `(sample_id, column_name)` keep the **last** value (see `_reshape_sdrf_long_to_wide_format`).
- Column names from SDRF are normalized: spaces → underscores, parentheses stripped (`organism_part_(group)` → `organism_part_group`). Downstream `obs` column names follow this rule.
- TODO: enumerate remaining per-accession quirks (mixed-type factor columns, missing SDRF rows, multi-pipeline TPM variants) as they are encountered.

## 5. Output contract

- **Object:** `anndata.AnnData` optionally saved to `<output_path>.h5ad`.
- **Shape:** `obs = samples`, `var = genes` (`= transcripts` for a transcript-level
  export — see below). `X` is `float32`, **samples x var** — the source TSV is
  genes/transcripts x samples and the builder transposes it (`df[sample_cols].values.T`).
- **`obs`:** indexed by `sample_id`; columns are SDRF characteristics / factors after the long → wide reshape (e.g. `organism`, `tissue`, `cell_type`, ...).
- **`var`:** indexed by `gene_id` for a gene-level export (which file `hvantk reprocess` builds from: § 4). For a **transcript-level** export (the real Expression Atlas per-transcript shape — see § 4) built with the default `level=transcript`, indexed by `transcript_id` instead, with `"Gene ID"` carried as a regular column so the gene id is not lost. Built with `level=gene`, a transcript-level export is summed onto the gene id and indexed by `gene_id` like a gene-level build (see § 4). Includes a `Gene Name` column when the source TSV had one.
- **`uns["column_summary"]`:** per-`obs`-column summary annotated by `annotate_column_summary_ad`.
- **Provenance:** stamped on the returned `ExpressionMatrix` via `ctx.provenance(schema_id="expression-atlas-dataset-v1")`; persisted by the platform as a sidecar `.provenance.json`.

## 6. hvantk integration points

- **Builder:** `build_expression_atlas` (signature `(parsed_input, ctx, **params) -> ExpressionMatrix`) in `hvantk/skills/expression_atlas/builder.py`.
- **SDRF / matrix helpers:** `hvantk/skills/expression_atlas/shared/expression_atlas.py`.
- **Dataset / collection classes:** `ExpressionAtlasDataset`, `ExpressionAtlasDatasetCollection` in `hvantk/skills/expression_atlas/shared/datasets.py`.
- **Downloader CLI:** `download_experiments` in `hvantk/skills/expression_atlas/cli.py` (declared in the manifest's `cli:` block as `expression-atlas-download`; the loader strips the `-download` suffix and binds it under the `download` group, so the invocation is `hvantk download expression-atlas`).
- **Lifecycle entry point:** `download_dataset` in `hvantk/skills/expression_atlas/cli.py`.
- **Build CLI:** `hvantk reprocess expression-atlas:dataset --raw-dir <dir> --output <path>.h5ad` (delegates to the plugin builder; pass builder kwargs via `--plugin-arg key=value`).
- **Plugin manifest:** `hvantk/skills/expression_atlas/plugin.yaml` (drives loader registration; compound dataset key `expression-atlas:dataset`).
- **Drift probe:** `fetch_fingerprint` in `hvantk/skills/expression_atlas/drift_probe.py`, compared against `tests/drift_fingerprint.json` by `hvantk drift expression-atlas:dataset` (see § 12 of `_conventions`).
- **Tests:** `hvantk/skills/expression_atlas/tests/` (downloader unit, drift-probe sanity, and
  `test_builder.py` — builder round-trip against a fixture derived from the real 320-column
  export, plus the column-inference regression guards).

Read the existing files at these paths as ground truth for shape. This skill does not restate code.

## 7. Workflow steps

When invoked to build or update a single Expression Atlas experiment:

1. **Resolve raw paths.** Either download via
   `hvantk download expression-atlas --accession <E-XXXX-N> --download-path data/expression-atlas`,
   or let `hvantk reprocess expression-atlas:dataset` run the lifecycle download stage itself (download → builder).
2. **Build:** `hvantk reprocess expression-atlas:dataset --raw-dir <dir> --output <out>.h5ad --skip-download` (the raw dir already holds step 1's download; without `--skip-download`, pass `--plugin-arg accession=<E-XXXX-N>` too, since `hvantk reprocess` forwards it to the download stage).
   - The builder parses the SDRF, transposes the expression matrix, attaches per-sample metadata into `obs`, annotates provenance, and writes `.h5ad`.
   - It builds from the TPM file whose columns are the SDRF's runs (the **gene-level** one when both files qualify). For E-MTAB-6798 that is the **transcript-level** file, since its gene-level file has a column per assay group (§4).
   - A transcript-level build needs `--plugin-arg level=gene` to collapse it to one row per gene before handing the output to `hvantk expression markers`/`summarize` or the PTM constraint loader, neither of which accepts a transcript-keyed matrix; a gene-level build already works without it. `hvantk annotate prepare`'s matrix source type needs more than `level=gene` either way (§4): it resolves `var_names` as HGNC gene symbols, not Ensembl gene ids.
3. **Validate:** run `pytest hvantk/skills/expression_atlas/tests`. That covers the seeded round-trip (`test_builder.py`, schema + sample rows + `n_obs`/`n_vars`) alongside the offline downloader and drift-probe tests. Regenerate snapshots only for an intentional change: `pytest hvantk/skills/expression_atlas/tests/test_builder.py --regenerate-snapshots`. Building from a real, untruncated download works (#342 fixed); verified against the full 320-column `E-MTAB-6798` file.

## 8. Update playbook

`hvantk drift expression-atlas:dataset` runs a real probe (`drift_probe.py`, `PROBE_VERSION = 1`), but it is atlas-wide, not per-accession: it fetches the single global experiment index (`https://www.ebi.ac.uk/gxa/json/experiments`), hashes the sorted `(experimentAccession, lastUpdate)` list over every experiment EBI hosts, and reports drift when *any* experiment anywhere in the atlas changes (`source_version` is just the experiment count, e.g. `"4562 experiments"`). A clean result means nothing moved atlas-wide; a drifted one does not say which accession moved, so updating one tracked accession is still a manual procedure:

1. For each tracked accession in `hvantk/skills/expression_atlas/catalog/datasets.json` (filter `data_source == "Expression_Atlas"`), re-download and compare by hand: `hvantk download expression-atlas --accession <E-XXXX-N> --download-path <dir>`.
2. Rebuild via `hvantk reprocess expression-atlas:dataset --raw-dir <dir> --output <out>.h5ad`, and diff the new AnnData against the committed snapshots in `tests/snapshots/` (§ 9).
3. If the SDRF column set changed, document the new factor in § 4.
4. Once the atlas-wide move is accounted for, accept it with `hvantk drift --regenerate expression-atlas:dataset --accepted-in "PR #N"` (`_conventions` § 12).

## 9. Validation contract

Per `_conventions` § 9:

- **fixture:** `hvantk/skills/expression_atlas/tests/testdata/raw/expression-atlas/` — seeded. `E-MTAB-6798-transcripts-tpms.tsv` (21 transcript rows / 20 distinct genes x 4 samples, ~1.2 KB — one gene carries two transcripts, so the fixture exercises the repeated-gene-id shape a real transcripts-tpms export has, #349) + `E-MTAB-6798.condensed-sdrf.tsv` (the same 4 sample IDs, ~3.3 KB), derived by truncation from the real upstream E-MTAB-6798 files (recipe at the top of `tests/test_builder.py`).
- **schema_snapshot:** `hvantk/skills/expression_atlas/tests/snapshots/schema.json` — seeded (4 obs x 21 vars).
- **row_snapshot:** `hvantk/skills/expression_atlas/tests/snapshots/sample_rows.json` — seeded.
- **drift_fingerprint:** `hvantk/skills/expression_atlas/tests/drift_fingerprint.json` — the expected fingerprint compared by `hvantk drift` (see § 12).
- **command:** `pytest hvantk/skills/expression_atlas/tests`.

The plugin manifest already declares these paths so the loader contract holds. `tests/test_builder.py` now exercises `build_expression_atlas` end-to-end against the committed fixture and asserts both snapshots plus `n_obs`/`n_vars`; regenerate via `pytest hvantk/skills/expression_atlas/tests/test_builder.py --regenerate-snapshots`. No Hail is required — the artifact is AnnData-backed.

The fixture is real E-MTAB-6798 data (CC BY 4.0); see its `THIRD_PARTY_DATA.md` entry.

> **The fixture is production-shaped.** It keeps the real three-column header
> (`Gene ID`, `Gene Name`, `GeneID`), so the round-trip test exercises the layout an
> actual download has — not a sanitised one. It previously dropped the `GeneID`
> transcript column to work around #342; that bug is fixed and the workaround is gone.
> It also previously kept only one row per gene, which meant the snapshot never
> exercised a repeated gene id; it now carries a real second transcript row for one
> gene (#349), so `var` is indexed by `transcript_id` and `"Gene ID"` repeats across
> two rows exactly as a real transcripts-tpms export would produce.
> `tests/test_builder.py` additionally pins these behaviours directly: the transcript
> column must key `var` and not become a third sample, a detected transcript column
> must leave `var_names` unique with the gene id demoted to a column, an
> all-missing sample column must stay a sample, and (#383) `level=gene` against the
> same two-transcript gene sums its transcripts' values per sample and re-keys `var`
> by the gene id (an all-missing sample or gene sums to `NaN`, not `0.0`; a
> transcript with no gene id raises). `hvantk/tests/test_reprocess_expression_atlas_cli.py`
> -- not under this skill's own `tests/`, since it imports `hvantk.tools.plugins.reprocess_cli`
> and `skills/` must not import `tools/` -- drives the real `reprocess` command against
> this fixture: with and without `level=gene`; with a fabricated per-run gene-level
> file alongside the fixture's transcript-level one, which
> `_resolve_expression_atlas_inputs` prefers; with an assay-group gene-level file
> (`g1`, `g2` columns, the shape E-MTAB-6798 ships), which it passes over for the
> transcript-level one; and two mismatched/ambiguous raw-dir shapes that must raise
> rather than silently combine or guess.
