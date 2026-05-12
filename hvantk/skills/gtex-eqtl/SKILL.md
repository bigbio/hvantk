---
name: hvantk:resource-gtex-eqtl
description: Build a Hail Table from GTEx cis-eQTL summary statistics (per-tissue significant variant-gene pairs) for qtlcascade and variant annotation.
status: provisional
backend: hail
domain: qtl
---

# GTEx eQTL (cis-QTL, significant variant-gene pairs)

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes every convention there.

## 1. Status & scope

- **Status:** provisional. The builder, helper functions, CLI, and registry adapter pre-exist this skill; this skill is the design contract and update reference.
- **In scope:** GTEx v11 cis-eQTL **significant variant-gene pairs** (`*.v11.eQTLs.signif_pairs.parquet`) — one parquet file per tissue. Verified against the v11 Liver release sliced to an 800-row fixture.
- **Anchored variant:** v11 / parquet. The same builder dispatches to v8 (TSV) and eqtlgen (TSV, different schema); see § 3. Triple-keyed output is identical across sources.
- **Out of scope:** downloader (GTEx distribution is via the GTEx portal / Google Cloud — manual acquisition per `_conventions` § 11); GTEx eGenes / sGenes / apaGenes / SuSiE fine-mapping files (different schemas, would need separate skills); GTEx allpairs files (>1 GB per tissue — documentation-only per the >1 GB downloader rule); cis-sQTL / cis-apaQTL / smallRNA cis-eQTL (the signif_pairs schema is shared but the `phenotype_id` semantics differ — see § 4 Gap 2); cross-tissue meta-analysis or trans-eQTLs.

## 2. Source identity

- **Provider:** GTEx Consortium. **Variant pinned by this skill:** v11 / cis-eQTL / signif_pairs / per-tissue parquet.
- **Catalog entry:** present. `hvantk/resources/registry/genomics/datasets.json` contains `GTEx_v11_eQTL_signif_pairs`. URLs / cadence / license / citation live in the registry — not here.
- **Source schema documentation:** `local/data/qtl_data/README_eQTL_v11.txt` (ships with the GTEx v11 release). The skill cites this README as the authoritative format reference; do NOT restate column definitions here.

Stable note (not in catalog): GTEx ships per-tissue parquet files named `<Tissue>.v11.eQTLs.signif_pairs.parquet`. The builder scans a directory and infers tissue from the filename prefix before the first dot (e.g., `Liver.v11.eQTLs.signif_pairs.parquet` → `Liver`). Single-file input is also accepted.

## 3. Backend choice + reasoning

**`backend: hail`, `domain: qtl`.** Per `_conventions` § 3, variant-level annotations key on `(locus, alleles)` — but eQTL summary stats have an additional gene axis: one variant can be an eQTL for multiple genes. The builder uses the **triple key** `(locus, alleles, gene_id)` to prevent information loss and enable correct cascade joins in `hvantk/qtlcascade/`. This is the first skill to exercise the triple-key shape; `_conventions` § 3 declares it but no prior skill anchored it.

**Multi-source dispatch under one builder.** `create_eqtl_tb` accepts `source ∈ {"gtex_v11", "gtex_v8", "eqtlgen"}` (see `hvantk/qtlcascade/constants.py:EQTL_SOURCES`). Each source uses a different import helper:

| Source | Import helper | Raw format | Anchored by skill? |
|---|---|---|---|
| `gtex_v11` (default) | `_import_eqtl_gtex_parquet` | Per-tissue parquet via `spark.read.parquet` → `hl.Table.from_spark` | **Yes** |
| `gtex_v8` | `_import_eqtl_gtex_tsv` | Per-tissue gzipped TSV via `hl.import_table` | No — alternative |
| `eqtlgen` | `_import_eqtl_eqtlgen` | Single TSV with different column names | No — separate provider |

The post-import `transform_func` is shared across all three: variant ID parsing, gene-version stripping, p-value filtering, triple-keying. Output schema is identical regardless of dispatch source.

**Spark/Parquet bridge note.** v11 ingestion uses `pyspark.sql.SparkSession.builder.getOrCreate()` → `spark.read.parquet(...)` → `hl.Table.from_spark(sdf)`. This is the first hvantk skill to document this path. The Spark session is reachable from within Hail (Hail itself runs on Spark); `_import_eqtl_gtex_parquet` calls `Path(fp).resolve()` because Spark requires absolute paths for parquet (per the builder docstring's "prototype lesson #4" annotation).

## 4. Raw format & gotchas

The v11 `*.signif_pairs.parquet` actually contains 12 columns. Cite `local/data/qtl_data/README_eQTL_v11.txt` for the canonical schema, but be aware of two README/file drifts (observed when slicing the Liver fixture):

**Drift 1: column name `start_distance`, not `tss_distance`.** The README documents column 4 as `tss_distance`; the actual parquet column is named `start_distance`. The builder does NOT read this column (only `phenotype_id`, `variant_id`, `slope`, `slope_se`, `pval_nominal`, and optionally `maf` — see Gap 1 below), so the output Hail Table is unaffected. Documentation drift, harmless.

**Drift 2: no `gene_id` column in signif_pairs.parquet.** The README schema for signif_pairs lists both `gene_id` and `phenotype_id` as separate columns. The actual file has only `phenotype_id`. For eQTL files this is fine (`phenotype_id ≡ gene_id` for eQTLs per the README), but **sQTL / apaQTL files use `phenotype_id` to encode intron coordinates or transcript-cluster IDs**, NOT bare gene IDs. The builder reads `phenotype_id` as the gene-ID source unconditionally — applying it to sQTL/apaQTL signif_pairs without modification would produce wrong gene keys. The skill scope (§ 1) is cis-eQTL only; sQTL/apaQTL would need a separate skill or a builder branch.

**Gap 1 (real builder bug, TODO follow-up): `maf` vs `af`.** `_import_eqtl_gtex_parquet` (in `hvantk/tables/table_builders.py:1563-1567`) probes for a column named `maf` and falls back to `hl.missing(hl.tfloat64)` if absent. v11 signif_pairs has `af` (ALT allele frequency, in-sample), **not `maf`** (minor allele frequency). Output `maf` field is therefore always null for v11 input. Real bug — but the round-trip snapshot will fix this null-`maf` shape as the recorded contract for the v11 source. Follow-up PR should:
- Read `af` from v11 parquet, optionally derive `maf` as `min(af, 1-af)` if downstream cares.
- Verify v8 path: v8 TSV signif_variant_gene_pairs may have `maf` directly — `_import_eqtl_gtex_tsv` should retain the existing probe.
- Cite: `hvantk/tables/table_builders.py:1563-1567` for the probe; `hvantk/tables/table_builders.py:1597-1620` for the v8 select.

**Gene ID versioning.** `gene_id` (read from `phenotype_id`) is a versioned GENCODE/Ensembl ID like `ENSG00000268903.1`. The shared `transform_func` strips the version suffix via `_strip_ensembl_version(ht.gene_id_raw)` before keying, producing `ENSG00000268903`. This is required for cross-table joins (e.g., HGNC, Ensembl gene metrics) that key on unversioned Ensembl IDs.

**Variant ID format.** `variant_id` is `{chr}_{pos}_{ref}_{alt}_b38` (GRCh38). The shared `_parse_gtex_variant_id` helper parses this into `locus` (a Hail `locus<GRCh38>`) and `alleles` (a length-2 `array<str>`). The trailing `_b38` suffix is consumed.

**Empty-tissue scan gotcha.** `_import_eqtl_gtex_parquet` raises `FileNotFoundError` if no parquet files match the requested tissue, but matching is case-sensitive on the filename prefix. Spelling errors silently produce a "no files" error; defensive callers should validate the tissue name against the known GTEx tissue list.

## 5. Output contract

- **Object:** `hl.Table` checkpointed to `output_path` (a `.ht` directory).
- **Key:** `[locus, alleles, gene_id]` (triple key).
- **Globals:** `hvantk_metadata` set by `_create_table_base`.
- **Fields (post-transform):**
  - `locus: locus<GRCh38>` — parsed from `variant_id`.
  - `alleles: array<str>` (length 2) — `[ref, alt]` from `variant_id`.
  - `gene_id: str` — unversioned Ensembl ID (e.g., `ENSG00000268903`).
  - `beta: float64` — regression slope from `slope`.
  - `se: float64` — SE of slope from `slope_se`.
  - `p_value: float64` — nominal p-value from `pval_nominal`.
  - `maf: float64` — **always null for v11** (Gap 1). Populated for v8 and eqtlgen.
  - `tissue: str` — inferred from filename prefix (e.g., `Liver`).
  - `gene_symbol: str` — always null for v11/v8 (`hl.missing(hl.tstr)` in the parquet importer); populated for eqtlgen.
  - `source: str` — value of the `source` parameter (e.g., `gtex_v11`).
  - `is_cis: bool` — always `True` (this skill is cis-eQTL scoped).
- **Reference genome:** GRCh38. The `_parse_gtex_variant_id` helper constructs `hl.Locus(..., reference_genome=reference_genome)`.

**Snapshot key uniqueness.** For the v11 Liver signif_pairs source, `(locus, alleles, gene_id)` is empirically unique-in-table — one row per (variant, gene) pair. Test inlines the small key list (per `_conventions` § 9 post-#101 rule); no `sample_keys.json` maintained.

## 6. hvantk integration points

- **Builder:** `create_eqtl_tb` in `hvantk/tables/table_builders.py:1644`, dispatched via `import_func` to one of three source-specific helpers. Uses `_create_table_base()` (so `overwrite` / `export_tsv` / `fields` flow through).
- **Source-specific helpers** (private):
  - `_import_eqtl_gtex_parquet` in `hvantk/tables/table_builders.py:1540` (v11 / Spark / Parquet).
  - `_import_eqtl_gtex_tsv` in `hvantk/tables/table_builders.py:1578` (v8 / TSV).
  - `_import_eqtl_eqtlgen` (eqtlgen TSV, separate schema).
- **Source constants:** `EQTL_SOURCES = ("gtex_v11", "gtex_v8", "eqtlgen")` in `hvantk/qtlcascade/constants.py:77`.
- **Registry:** `TABLE_BUILDERS["eqtl"] = create_table_adapter("hvantk.tables.table_builders", "create_eqtl_tb")` in `hvantk/tables/registry.py`.
- **CLI:** `mktable_eqtl` in `hvantk/commands/make_table_cli.py:604` (command name `eqtl`), with `--source` flag (Click choice, default `gtex_v11`), `--tissue`, `--p-threshold`, plus the standard input/output/overwrite/export options.
- **Downstream consumer:** `hvantk/qtlcascade/` — the eQTL Hail Table is one half of the eQTL ⊕ pQTL cascade join.

## 7. Workflow steps

1. **Resolve raw path.** Caller passes either a single parquet file (`Liver.v11.eQTLs.signif_pairs.parquet`) or a directory containing per-tissue parquet files. Acquire from the GTEx portal (https://gtexportal.org) — manual download, no skill-side acquisition.
2. **Import (v11 path).** `_import_eqtl_gtex_parquet` scans `.parquet` files under `input_path`, optionally filters by `tissue`, reads each via `spark.read.parquet`, converts via `hl.Table.from_spark`, and unions the per-tissue tables. Source columns selected: `phenotype_id`, `variant_id`, `slope`, `slope_se`, `pval_nominal`, optionally `maf` (null in v11 — see Gap 1).
3. **Transform** (shared `transform_func` in `create_eqtl_tb`):
   - `_parse_gtex_variant_id(ht, "variant_id", reference_genome)` — yields `locus`, `alleles`, drops the raw `variant_id` after annotation.
   - `gene_id = _strip_ensembl_version(ht.gene_id_raw)` — strips `.N` version suffix.
   - `ht.drop("gene_id_raw", "variant_id")`.
   - If `p_threshold > 0`: `ht.filter(ht.p_value <= p_threshold)`. **Set `p_threshold=0` to retain all rows** (signif_pairs files are already pre-filtered to significant pairs by GTEx; this is what coloc / qtlcascade typically want).
   - `ht.annotate(source=source, is_cis=True)`.
   - `ht.key_by("locus", "alleles", "gene_id")`.
4. **Checkpoint + globals + optional TSV.** Handled by `_create_table_base` via `overwrite` and `export_tsv` kwargs.

**Default usage (CLI, single tissue):**

```bash
hvantk mktable eqtl \
    --raw-input /path/to/gtex_v11_eqtl/ \
    --output-ht /path/to/gtex_v11_eqtl_Liver.ht \
    --source gtex_v11 \
    --tissue Liver \
    --p-threshold 0 \
    --ref-genome GRCh38
```

## 8. Update playbook

GTEx releases major versions every few years (v8 → v9 → v10 → v11). Per release:

1. **Acquire** the new release's per-tissue signif_pairs files (manual download).
2. **Check format.** Compare the new release's `README_eQTL_*.txt` against the v11 schema documented in this skill. If columns added/removed: extend `_import_eqtl_gtex_parquet` and re-record the snapshot.
3. **Add a new source constant** to `EQTL_SOURCES` (e.g., `"gtex_v12"`); add a `gtex_v12` branch in `create_eqtl_tb.import_func` if the import shape diverges from v11; bump the registry entry's `accession` (`GTEx_v11_eQTL_signif_pairs` → `GTEx_v12_eQTL_signif_pairs`).
4. **Re-run round-trip (§ 9).** If schema changes, regenerate snapshots via `--regenerate-snapshots`.
5. **Resolve Gap 1 (maf / af).** Until the follow-up PR lands, every new v11+ release will continue to produce null `maf`. The follow-up should fix the probe to read `af` and optionally derive `maf`.

**Cross-source compatibility:** v8 TSV → v11 parquet was a breaking format change (text → columnar binary). v11 → future releases may break again. Anchor each release to its own source constant rather than reusing `gtex_v11` after the next migration.

## 9. Validation contract

Per `_conventions` § 9:

- **fixture:** `hvantk/tests/testdata/raw/gtex-eqtl/Liver.v11.eQTLs.signif_pairs.parquet`. 800 rows sliced from the v11 Liver source via `local/planning/skills-gtex-eqtl-fixture-slicer.py` (gitignored). ~36 KB. Preserves the parquet binary format (deterministic via `pq.Table.slice(0, N)`).
- **schema_snapshot:** `hvantk/tests/snapshots/gtex-eqtl/schema.json`.
- **row_snapshot:** `hvantk/tests/snapshots/gtex-eqtl/sample_rows.json`. Triple key `(locus, alleles, gene_id)` is unique-in-table for signif_pairs, so the test inlines a small key list and no `sample_keys.json` is maintained (per `_conventions` § 9 post-#101).
- **test_command:** `pytest hvantk/tests/test_gtex_eqtl_builder.py -m hail`.

Round-trip test asserts: builder idempotent with `overwrite=True`; checkpointed schema matches `schema.json`; deterministic sample-row slice matches `sample_rows.json`. Test uses `source="gtex_v11"`, `tissue="Liver"`, `p_threshold=0` (retain all rows in the fixture). Spark session is initialized on demand by `_import_eqtl_gtex_parquet`; the `hail_session` fixture ensures Hail (and therefore Spark) is up.

Regenerate via `--regenerate-snapshots` when:
- The shared `transform_func` adds / removes / renames fields.
- Gap 1 (maf/af) is fixed — the `maf` field will start carrying real values, breaking the current snapshot's null-`maf` contract.
- A new GTEx release is anchored that changes the input schema (per § 8).
