---
name: hvantk:resource-gwas-catalog
description: Build a Hail Table from the EBI GWAS Catalog v1.0 full-associations TSV for variant-annotation joins.
status: provisional
backend: hail
domain: variants
---

# GWAS Catalog (EBI) — v1.0 full associations

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes every convention there.

## 1. Status & scope

- **Status:** provisional. The builder, fixture, snapshots, and round-trip test now exist; this skill remains the design contract and update reference for future agents.
- **In scope:** raw GWAS Catalog v1.0 full-associations TSV → single Hail Table, checkpointed, keyed for variant joins.
- **Out of scope:** downloader; v1.0.2 schema (EFO URIs — possible follow-on); trait-aware enrichment / burden / PS-ROC; cross-source joins.

## 2. Source identity

- **Provider:** EBI GWAS Catalog. **Schema pinned by this skill:** v1.0 (34 columns, no `MAPPED_TRAIT_URI`).
- **Catalog entry:** present. `hvantk/resources/registry/genomics/datasets.json` contains `GWAS_Catalog_v1.0_e115_r2026-04-27`, and `hvantk/resources/catalog.yaml` records `datasets.genomics.count: 9`. URLs / cadence / license / citation live in the catalog and registry — not here.

## 3. Backend choice + reasoning

**`backend: hail`, `domain: variants`.** The catalog is a flat ~10⁶-row TSV consumed by joining onto variant tables (ClinVar, dbNSFP, gnomAD) keyed by `(locus, alleles)`. Hail gives `hl.parse_locus` over `CHR_ID` + `CHR_POS`, cheap key-joins into the existing variant ecosystem, and free checkpointing + globals via `_create_table_base`. A pandas builder would force callers to materialize ~1M rows per join. Domain is `variants` (not `genes`) — every row anchors to a SNP, not `MAPPED_GENE`.

Key by `(locus, alleles)` with sentinel ALT — `alleles = [<risk_allele>, "N"]`. The risk allele is real; the ALT is a synthetic `"N"` sentinel. Callers joining against ClinVar / gnomAD must normalize on REF/ALT separately; this skill does not do that work.

### Judgment call #1 — Keying

**Decision: A — `(locus, alleles)` with sentinel ALT.** Matches the standard hvantk variants convention; downstream normalization for cross-source joins is left to callers.

## 4. Raw format & gotchas

Single 556 MB tab-separated UTF-8 TSV inside a zip (`gwas-catalog-download-associations-v1.0-full.tsv`), ~1.1M rows, **34 columns** with whitespace-and-slash-laden headers (`DISEASE/TRAIT`, `STRONGEST SNP-RISK ALLELE`, `OR or BETA`, `95% CI (TEXT)`, `PLATFORM [SNPS PASSING QC]`). Use `hl.import_table(..., delimiter='\t', quote=None, missing='', impute=False)` — column names with spaces, mixed types, and sparse missing values make Hail's imputation brittle. Rename to snake_case in the transform.

Type coercions in transform (all string at import):

- `chr_pos` → `int32` (after multi-chrom cleanup — see #3).
- `p_value`, `pvalue_mlog`, `or_or_beta`, `risk_allele_frequency`, `upstream_gene_distance`, `downstream_gene_distance` → `float64`. `p_value` arrives in scientific notation (`3E-12`); `hl.float` parses directly.
- `intergenic`, `cnv` → `bool` from `0/1` or `Y/N` (confirm against fixture).
- `pubmedid`, `snp_id_current`, `merged` → `int32`; `missing=''` so blanks → `NA`.

### Gotchas (maintainer's empirical pre-inspection)

| Gotcha | Frequency | Where it bites |
|---|---|---|
| `STRONGEST SNP-RISK ALLELE = "rsXXX-?"` | ~30% of rows | Judgment call #2 |
| Multi-chrom / haplotype / interaction: `CHR_ID = "6;7"`, `"2;2"`, or `"1 x 10"` | small | `CHR_POS` int cast + locus parse fail — judgment call #3 |
| Same rsID × many trait rows (rs704 × 591 in 100k) | high; normal shape | Drives judgment call #1 |
| `P-VALUE` scientific notation | universal | `hl.float` handles it; benign |
| `OR or BETA` single column, mutually exclusive | universal | Carry as-is; callers disambiguate |
| No REF/ALT in raw file | universal | Drives judgment call #1 |

### Judgment call #2 — Missing-allele rows (`rsXXX-?`)

**Decision: A — Drop.** Filter rows where `STRONGEST SNP-RISK ALLELE` ends in `-?`. Accept ~30% data loss; every retained row carries a real risk allele.

### Judgment call #3 — Multi-chromosome / haplotype rows

**Decision: A — Drop.** Filter rows where `CHR_ID` is NOT a canonical single contig — i.e., does not match `^(chr)?(\d+|X|Y|MT?)$`. This drops `;`-separated multi-chromosome rows (`"6;7"`, `"2;2"`), interaction-pair rows (`"1 x 10"` with `x`-separator), and any other malformed contig shape encountered in the wild. Haplotype / interaction associations are silently lost; one locus per row keeps the `(locus, alleles)` key clean. Implement as `ht.filter(ht.chr_id.matches("^(chr)?(\\d+|X|Y|MT?)$"))`.

### Judgment call #4 — Trait ontology

**Decision: A — N/A for v1.0.** v1.0 has no `MAPPED_TRAIT` / `MAPPED_TRAIT_URI` — carry `DISEASE/TRAIT` as plain string. No EFO work in this skill.

## 5. Output contract

- **Object:** `hl.Table` checkpointed to `output_path` (a `.ht` directory).
- **Key:** `[locus, alleles]` where `alleles = [<risk_allele>, "N"]` (sentinel ALT).
- **Globals:** `hvantk_metadata` set by `_create_table_base`.
- **Fields:** snake_case 1:1 renames of the surviving raw 34 columns. Do not drop raw columns; let callers `select()`. No judgment-call flag columns (`is_haplotype`, `has_risk_allele`) — both disqualifying conditions are filtered upstream.
- **Reference genome:** `GRCh38`. No liftover.

## 6. hvantk integration points

- **Builder:** `create_gwas_catalog_tb` in `hvantk/skills/gwas_catalog/builder.py`, via `_create_table_base()`. Signature per conventions §5. A re-export shim in `hvantk/tables/table_builders.py` keeps the old import path working.
- **Registry:** registered via the plugin manifest at `hvantk/skills/gwas_catalog/plugin.yaml` (`gwas-catalog:associations`). Plugin discovery wires it into `TABLE_BUILDERS` at import time.
- **CLI:** `mktable_gwas_catalog` in `hvantk/tools/build/make_table_cli.py` (command name `gwas-catalog`), decorated with `@_raw_input_opt` / `@_output_ht_opt`.
- **Catalog wiring:** see §2. **Downloader:** out of scope.

## 7. Workflow steps

1. **Resolve raw path.** Caller passes the unzipped TSV path; the builder does not unzip.
2. **Import.** `hl.import_table(input_path, delimiter='\t', quote=None, missing='', impute=False)`.
3. **Transform** (`transform_func` passed to `_create_table_base`):
   - Rename the 34 columns to snake_case 1:1 (do not drop).
   - Filter rows where `STRONGEST SNP-RISK ALLELE` ends in `-?` (judgment call #2).
   - Filter rows where `CHR_ID` is not canonical (no match against `^(chr)?(\d+|X|Y|MT?)$`) — drops `;`-separated and `x`-separated malformed shapes (judgment call #3).
   - Cast numerics (`chr_pos`, `p_value`, `pvalue_mlog`, `or_or_beta`, `risk_allele_frequency`, `upstream_gene_distance`, `downstream_gene_distance`).
   - Recode `chr_id` to GRCh38 contig form: the catalog ships bare contigs (`"7"`, `"12"`), but Hail's `GRCh38` reference expects `"chr7"` etc. Use `contig = hl.if_else(chr_id.startswith("chr"), chr_id, "chr" + chr_id)`, then `locus = hl.parse_locus(contig + ':' + hl.str(chr_pos), reference_genome=reference_genome)`.
   - Extract `risk_allele` from `STRONGEST SNP-RISK ALLELE` (split on `-`, take suffix).
   - Construct `alleles = [risk_allele, "N"]` and `key_by(locus, alleles)`.
4. **Checkpoint + globals + optional TSV.** Handled by `_create_table_base` via `overwrite` and `export_tsv` kwargs.

## 8. Update playbook

Releases ~quarterly (tags like `e116_r2026-08-xx`). Per release:

1. Update the GWAS Catalog entry in `registry/genomics/datasets.json` (`last_updated`, file path). Bump accession only on schema change.
2. Re-run the round-trip test (§9). If it passes, no builder change.
3. On schema change: bump accession suffix, update the rename map in the builder, regenerate snapshots with `pytest --regenerate-snapshots`, revisit §3–§4 judgment calls.
4. If EBI publishes v1.0.2 alongside v1.0: revisit judgment call #4 (option B).

## 9. Validation contract

Per conventions §9:

- **fixture:** `hvantk/skills/gwas_catalog/tests/testdata/raw/gwas-catalog/gwas-catalog-sample.tsv`. Agent #2 slices from the 556 MB local file (path from orchestrator). Target ≤ 100 KB, ~50–100 rows. Fixture must exercise **surviving rows only** — no `STRONGEST SNP-RISK ALLELE` ending in `-?`, no `CHR_ID` containing `;`. Must include: multi-trait-per-variant rows (same `(locus, alleles)` × multiple traits), scientific-notation p-values, and at least one row with `OR or BETA` populated plus one without.
- **schema_snapshot:** `hvantk/skills/gwas_catalog/tests/snapshots/schema.json`.
- **row_snapshot:** `hvantk/skills/gwas_catalog/tests/snapshots/sample_rows.json`.
- **sample_keys:** `hvantk/skills/gwas_catalog/tests/snapshots/sample_keys.json`. Lists the `(locus, alleles)` keys used for `row_snapshot` assertions. These keys MUST be unique-in-table — `_snapshot_utils.collect_sample_rows` does not deduplicate, so a duplicated key produces non-deterministic snapshots. Multi-trait-per-variant rows in this catalog routinely share keys; the snapshot subset must use singleton-key rows.
- **test_command:** `pytest hvantk/skills/gwas_catalog/tests -m hail`.

Round-trip test asserts: builder idempotent with `overwrite=True`; checkpointed schema matches `schema.json`; deterministic sorted row slice matches `sample_rows.json`. Regenerate via `--regenerate-snapshots` when a judgment call resolves or the schema changes.
