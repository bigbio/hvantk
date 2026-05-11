---
name: hvantk:resource-msigdb
description: Build a Hail Table from an MSigDB GMT gene-set file (e.g., C2 Canonical Pathways) for enrichment / burden / overlap analyses.
status: provisional
backend: hail
domain: mapping
---

# MSigDB (Molecular Signatures Database)

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes every convention there.

## 1. Status & scope

- **Status:** provisional. The builder, fixture, snapshots, and round-trip test exist; this skill is the design contract and update reference.
- **In scope:** any GMT file from MSigDB → single Hail Table keyed by `set_name`. Verified against the C2 Canonical Pathways human gene-symbols collection (`c2.cp.v2026.1.Hs.symbols.gmt`, 4,115 sets, 1.6 MB).
- **Out of scope:** downloader (MSigDB requires login + license click-through, so per `_conventions` § 11 acquisition is manual); other MSigDB collections (H, C1, C3-C8) — the builder is collection-agnostic but each collection should be tracked in the catalog separately if onboarded; gene-symbol normalization / alias resolution (use `hvantk/data/gene_mapper.py` downstream); cross-format variants beyond GMT (GMX, XML).

## 2. Source identity

- **Provider:** Broad Institute MSigDB. **Variant pinned by this skill's catalog entry:** C2 / CP / human / gene symbols / v2026.1.
- **Catalog entry:** present. `hvantk/resources/registry/genomics/datasets.json` contains `MSigDB_C2_CP_v2026.1.Hs.symbols`, and `hvantk/resources/catalog.yaml` records `datasets.genomics.count: 10`. URLs / cadence / license / citation live in the catalog and registry — not here.

Stable note (not in catalog): MSigDB ships per-collection GMT files. The GMT format is the same across collections, so this builder works for any MSigDB GMT, but each onboarded collection needs its own catalog entry to record license / version / file path.

## 3. Backend choice + reasoning

**`backend: hail`, `domain: mapping`.** A C2 CP GMT is ~4k rows × variable-width gene columns. Per `_conventions` § 3 "Lookup / mapping" allows a Hail Table or pandas DataFrame. Hail wins here because downstream consumers (enrichment / burden / overlap, e.g., `hvantk/enrichex/`) join against Hail Tables keyed on gene symbols — producing a Hail Table avoids re-materialization at every join site, mirroring the HGNC decision. Key by `set_name` (string, unique-in-file).

> Catalog placement note: this skill's catalog entry lives in `registry/genomics/datasets.json`, not a dedicated `mapping/` registry directory (no such directory exists today). The decision is consistent with `GWAS_Catalog_v1.0_*` (also gene-symbol / variant-adjacent annotation curation) and avoids invasive changes to `catalog.yaml registry_structure` + `hvantk/resources/unified_registry.py`. Revisit if a `mapping/` domain is later introduced.

## 4. Raw format & gotchas

GMT is **tab-separated with variable-width rows**:

- Column 1: `set_name` (string, unique within a single GMT — e.g., `KEGG_APOPTOSIS`).
- Column 2: description text. For MSigDB-issued GMTs this is **always** a gsea-msigdb URL like `https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/KEGG_APOPTOSIS`. The spec allows arbitrary text, so the builder preserves it as-is in `source_url` without parsing. Confirmed against the local 4,115-row source: 100% of rows had a `https://www.gsea-msigdb.org/` URL in column 2.
- Columns 3..N: gene members. For human gene-symbol GMTs (`.Hs.symbols.gmt`), these are HGNC-approved gene symbols. For other variants (`.entrez.gmt`, `.Mm.symbols.gmt`) the contents differ; the builder is symbol-agnostic — it carries strings.

**Variable-width-row gotcha (this is the key deviation):** `hl.import_table` rejects rows with inconsistent column counts. Instead, the builder uses **`hl.import_lines`** (one row per line, single `text: str` field) and splits on `\t` inside the transform. Gene members are sliced as `parts[2:]` into an `array<str>`. Set lengths in the source range from 5 (MSigDB minimum) to 1,497 genes; Hail arrays handle this without issue.

**Empty-line gotcha:** `hl.import_lines` does yield rows for blank lines (with `text == ""`). The transform filters `set_name == ""` after split as a defensive guard.

**No type coercions needed** — every field is string or `array<str>` by design.

## 5. Output contract

- **Object:** `hl.Table` checkpointed to `output_path` (a `.ht` directory).
- **Key:** `[set_name]` (string, unique-in-table).
- **Globals:** `hvantk_metadata` set by `_create_table_base`.
- **Fields:**
  - `set_name: str` — gene-set identifier (e.g., `KEGG_APOPTOSIS`).
  - `source_url: str` — GMT column 2, verbatim. For MSigDB-issued GMTs this is a `https://www.gsea-msigdb.org/...` URL.
  - `genes: array<str>` — gene members. Order is preserved from the source file.
- **Reference genome:** N/A. Gene-set membership is genome-independent; the catalog entry records `GRCh38` only because the schema requires it.

Per `_conventions` § 9, set names are unique-in-table for a single GMT, so **no `sample_keys.json` is maintained** — the snapshot test reads keys directly from a small inline list. The post-#101 conventions §9 "unique key, no workaround" rule applies cleanly here.

## 6. hvantk integration points

- **Builder:** `create_msigdb_tb` in `hvantk/tables/table_builders.py`, via `_create_table_base()`. Signature per `_conventions` § 5: `(input_path, output_path, overwrite=False, export_tsv=False)`.
- **Registry:** `TABLE_BUILDERS["msigdb"] = create_table_adapter("hvantk.tables.table_builders", "create_msigdb_tb")` in `hvantk/tables/registry.py`.
- **CLI:** `mktable_msigdb` in `hvantk/commands/make_table_cli.py` (command name `msigdb`), decorated with `@_raw_input_opt` / `@_output_ht_opt` / `@_overwrite_opt` / `@_export_tsv_opt`. No `--ref-genome` flag (gene-set membership is genome-independent).
- **Catalog wiring:** see § 2. **Downloader:** out of scope (manual acquisition).

## 7. Workflow steps

1. **Resolve raw path.** Caller passes the unzipped `.gmt` path. Acquire from <https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp> (login required).
2. **Import.** `hl.import_lines(input_path, min_partitions=4)` — yields one row per line with `text: str`.
3. **Transform** (`transform_func` passed to `_create_table_base`):
   - `parts = ht.text.split("\t")`.
   - `set_name = parts[0]`, `source_url = parts[1]`, `genes = parts[2:]` (Hail array slice).
   - `ht.filter(ht.set_name != "")` — defensive against blank lines.
   - `ht.key_by("set_name")`.
4. **Checkpoint + globals + optional TSV.** Handled by `_create_table_base` via `overwrite` and `export_tsv` kwargs.

## 8. Update playbook

MSigDB releases ~annually (versioned `v<year>.<n>`, e.g., `v2026.1`, `v2025.1`). Per release:

1. Acquire the new GMT (manual download). Update the file path / version in the corresponding `registry/genomics/datasets.json` entry; bump `accession` (`MSigDB_C2_CP_v2026.1.Hs.symbols` → `MSigDB_C2_CP_v2027.1.Hs.symbols`).
2. Re-run the round-trip test (§ 9). If it passes, no builder change.
3. The GMT format has been stable for ~15 years; column 1 / column 2 / variable-tail shape has not changed. If MSigDB ever changes the description column (column 2) away from a URL, the `source_url` field name becomes misleading — rename to `description` and update this skill.
4. To onboard a different collection (e.g., C5 GO, H Hallmark): add a new catalog entry with the new accession; the same `create_msigdb_tb` builder works without modification. Add a parallel fixture and snapshot directory if the new collection has structural quirks (e.g., GMTs with embedded null bytes).

## 9. Validation contract

Per `_conventions` § 9:

- **fixture:** `hvantk/tests/testdata/raw/msigdb/c2.cp-sample.gmt`. 20 gene sets, ~24 KB, sliced from the v2026.1 C2 CP source via `local/planning/skills-tier3.5-msigdb-fixture-slicer.py` (gitignored — see `local/planning/`). Exercises the short edge (size 5: BIOCARTA, SA), medium sets (60-330 genes), a long set (`REACTOME_CELL_CYCLE`, 688 genes), and the extra-long tail (`REACTOME_POST_TRANSLATIONAL_PROTEIN_MODIFICATION`, 1,497 genes). All 20 fixture rows have a `https://www.gsea-msigdb.org/` URL in column 2 (matches the live-file invariant).
- **schema_snapshot:** `hvantk/tests/snapshots/msigdb/schema.json`.
- **row_snapshot:** `hvantk/tests/snapshots/msigdb/sample_rows.json`. `set_name` keys are unique-in-table, so no `sample_keys.json` is maintained per `_conventions` § 9 (post-#101). The round-trip test inlines the small key list.
- **test_command:** `pytest hvantk/tests/test_msigdb_builder.py -m hail`.

Round-trip test asserts: builder idempotent with `overwrite=True`; checkpointed schema matches `schema.json`; deterministic sorted row slice matches `sample_rows.json`. Regenerate via `--regenerate-snapshots` when the schema changes (rare — see § 8).
