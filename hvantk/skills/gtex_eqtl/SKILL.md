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

- **Status:** provisional. The builder, helper functions, and CLI pre-exist this skill; this skill is the design contract and update reference.
- **In scope:** GTEx v11 cis-eQTL **significant variant-gene pairs** (`*.v11.eQTLs.signif_pairs.parquet`) — one parquet file per tissue. Verified against the v11 Liver release sliced to an 800-row fixture.
- **Anchored variant:** v11 / parquet. The same builder dispatches to v8 (TSV) and eqtlgen (TSV, different schema); see § 3. Triple-keyed output is identical across sources.
- **Out of scope:** downloader (GTEx distribution is via the GTEx portal / Google Cloud — manual acquisition per `_conventions` § 11); GTEx eGenes / sGenes / apaGenes / SuSiE fine-mapping files (different schemas, would need separate skills); GTEx allpairs files (>1 GB per tissue — documentation-only per the >1 GB downloader rule); cis-sQTL / cis-apaQTL / smallRNA cis-eQTL (the signif_pairs schema is shared but the `phenotype_id` semantics differ — see § 4 Gap 2); cross-tissue meta-analysis or trans-eQTLs.

## 2. Source identity

- **Provider:** GTEx Consortium. **Variant pinned by this skill:** v11 / cis-eQTL / signif_pairs / per-tissue parquet.
- **Catalog entry:** present. `hvantk/resources/registry/genomics/datasets.json` contains `GTEx_v11_eQTL_signif_pairs`. URLs / cadence / license / citation live in the registry — not here.
- **Source schema documentation:** the `README_eQTL_v11.txt` that GTEx ships inside the v11 cis-QTL release archive (download via the GTEx portal). The skill cites this README as the authoritative format reference; do NOT restate column definitions here.

Stable note (not in catalog): GTEx ships per-tissue parquet files named `<Tissue>.v11.eQTLs.signif_pairs.parquet`. The builder scans a directory and infers tissue from the filename prefix before the first dot (e.g., `Liver.v11.eQTLs.signif_pairs.parquet` → `Liver`). Single-file input is also accepted.

## 3. Backend choice + reasoning

**`backend: hail`, `domain: qtl`.** Per `_conventions` § 3, variant-level annotations key on `(locus, alleles)` — but eQTL summary stats have an additional gene axis: one variant can be an eQTL for multiple genes. The builder uses the **triple key** `(locus, alleles, gene_id)` to prevent information loss and enable correct cascade joins in `hvantk/qtlcascade/`. This is the first skill to exercise the triple-key shape; `_conventions` § 3 declares it but no prior skill anchored it.

**Multi-source dispatch under one builder.** `build_eqtl_associations` accepts `source ∈ {"gtex_v11", "gtex_v8", "eqtlgen"}` (see `hvantk/skills/gtex_eqtl/shared/constants.py:EQTL_SOURCES`). Each source uses a different import helper:

| Source | Import helper | Raw format | Anchored by skill? |
|---|---|---|---|
| `gtex_v11` (default) | `_import_gtex_parquet` | Per-tissue parquet via `spark.read.parquet` → `hl.Table.from_spark` | **Yes** |
| `gtex_v8` | `_import_gtex_tsv` | Per-tissue gzipped TSV via `hl.import_table` | No — alternative |
| `eqtlgen` | `_import_eqtlgen` | Single TSV with different column names | No — separate provider |

The post-import transform is shared across all three (inline in `build_eqtl_associations`): variant ID parsing, gene-version stripping, p-value filtering, triple-keying. Output schema is identical regardless of dispatch source.

**Spark/Parquet bridge note.** v11 ingestion uses `pyspark.sql.SparkSession.builder.getOrCreate()` → `spark.read.parquet(...)` → `hl.Table.from_spark(sdf)`. This is the first hvantk skill to document this path. The Spark session is reachable from within Hail (Hail itself runs on Spark); `_import_gtex_parquet` calls `Path(fp).resolve()` because Spark requires absolute paths for parquet (per the builder docstring's "prototype lesson #4" annotation).

## 4. Raw format & gotchas

The v11 `*.signif_pairs.parquet` actually contains 12 columns. Cite the upstream `README_eQTL_v11.txt` (shipped in the GTEx v11 cis-QTL release archive) for the canonical schema, but be aware of two README/file drifts (observed when slicing the Liver fixture):

**Drift 1: column name `start_distance`, not `tss_distance`.** The README documents column 4 as `tss_distance`; the actual parquet column is named `start_distance`. The builder does NOT read this column (only `phenotype_id`, `variant_id`, `slope`, `slope_se`, `pval_nominal`, and `af`/`maf` — see below), so the output Hail Table is unaffected. Documentation drift, harmless.

**Drift 2: no `gene_id` column in signif_pairs.parquet.** The README schema for signif_pairs lists both `gene_id` and `phenotype_id` as separate columns. The actual file has only `phenotype_id`. For eQTL files this is fine (`phenotype_id ≡ gene_id` for eQTLs per the README), but **sQTL / apaQTL files use `phenotype_id` to encode intron coordinates or transcript-cluster IDs**, NOT bare gene IDs. The builder reads `phenotype_id` as the gene-ID source unconditionally — applying it to sQTL/apaQTL signif_pairs without modification would produce wrong gene keys. The skill scope (§ 1) is cis-eQTL only; sQTL/apaQTL would need a separate skill or a builder branch.

**`maf` vs `af`.** `_import_gtex_parquet` probes for an `af` column (ALT allele frequency, in-sample) and falls back to `maf` then `hl.missing(hl.tfloat64)`. v11 signif_pairs has `af`, not `maf`. The v11 select reads `af` directly and **exposes both fields**:
- `af: float64` — read from the v11 parquet's `af` column. Preserves directional info (slope is reported per ALT allele; consumers wanting direction-aware analysis need `af`, not `maf`).
- `maf: float64` — derived as `hl.min(af, 1 - af)`. Half-symmetric: matches the conventional MAF definition for downstream coloc / frequency-filter use.

For v8 (`_import_gtex_tsv`), the source TSV carries `maf` directly (not `af`); both fields are set equal to the TSV's `maf` value as an approximation, since v8 doesn't carry directional info — a real caveat documented here: **v8 inputs lose direction**. Downstream code that depends on directional ALT frequency should use v11.

For eqtlgen, both `af` and `maf` remain `hl.missing` (the source distributes neither).

**Gene ID versioning.** `gene_id` (read from `phenotype_id`) is a versioned GENCODE/Ensembl ID like `ENSG00000268903.1`. The shared transform strips the version suffix via `strip_ensembl_version(ht.gene_id_raw)` before keying, producing `ENSG00000268903`. This is required for cross-table joins (e.g., HGNC, Ensembl gene metrics) that key on unversioned Ensembl IDs.

**Variant ID format.** `variant_id` is `{chr}_{pos}_{ref}_{alt}_b38` (GRCh38). The shared `parse_gtex_variant_id` helper parses this into `locus` (a Hail `locus<GRCh38>`) and `alleles` (a length-2 `array<str>`). The trailing `_b38` suffix is consumed.

**Empty-tissue scan gotcha.** `_import_gtex_parquet` raises `FileNotFoundError` if no parquet files match the requested tissue, but matching is case-sensitive on the filename prefix. Spelling errors silently produce a "no files" error; defensive callers should validate the tissue name against the known GTEx tissue list.

## 5. Output contract

- **Object:** an `AnnotationTable` wrapping a Hail Table (via `AnnotationTable.from_hail`). The `hvantk reprocess` runner checkpoints it to the `--output` `.ht` directory.
- **Key:** `[locus, alleles, gene_id]` (triple key).
- **Provenance:** stamped on the returned `AnnotationTable` via `ctx.provenance(schema_id="gtex-eqtl-eqtls-v1")`; persisted as a sidecar `.provenance.json`.
- **Fields (post-transform):**
  - `locus: locus<GRCh38>` — parsed from `variant_id`.
  - `alleles: array<str>` (length 2) — `[ref, alt]` from `variant_id`.
  - `gene_id: str` — unversioned Ensembl ID (e.g., `ENSG00000268903`).
  - `beta: float64` — regression slope from `slope`.
  - `se: float64` — SE of slope from `slope_se`.
  - `p_value: float64` — nominal p-value from `pval_nominal`.
  - `af: float64` — ALT allele frequency. For v11: read directly from the `af` parquet column. For v8: set equal to the v8 TSV's `maf` (approximation; v8 loses direction). For eqtlgen: null.
  - `maf: float64` — minor allele frequency, `min(af, 1-af)`. For v11: derived from `af`. For v8: read directly from the v8 TSV's `maf` column. For eqtlgen: null.
  - `tissue: str` — inferred from filename prefix (e.g., `Liver`).
  - `gene_symbol: str` — always null for v11/v8 (`hl.missing(hl.tstr)` in the parquet importer); populated for eqtlgen.
  - `source: str` — value of the `source` parameter (e.g., `gtex_v11`).
  - `is_cis: bool` — always `True` (this skill is cis-eQTL scoped).
- **Reference genome:** GRCh38. The `parse_gtex_variant_id` helper constructs `hl.Locus(..., reference_genome=reference_genome)`.

**Snapshot key uniqueness.** For the v11 Liver signif_pairs source, `(locus, alleles, gene_id)` is empirically unique-in-table — one row per (variant, gene) pair. Test inlines the small key list (per `_conventions` § 9 post-#101 rule); no `sample_keys.json` maintained.

## 6. hvantk integration points

- **Builder:** `build_eqtl_associations` in `hvantk/skills/gtex_eqtl/builder.py`, signature `(parsed_input, ctx, *, reference_genome, source, tissue, p_threshold, fields) -> AnnotationTable`. Dispatches on `source` to one of three source-specific import helpers, then applies the shared transform inline (no `_create_table_base`).
- **Source-specific helpers** (private, all in `hvantk/skills/gtex_eqtl/builder.py`):
  - `_import_gtex_parquet` (v11 / Spark / Parquet).
  - `_import_gtex_tsv` (v8 / TSV).
  - `_import_eqtlgen` (eqtlgen TSV, separate schema).
- **Shared QTL helpers:** `parse_gtex_variant_id`, `strip_ensembl_version`, `scan_tissue_files` in `hvantk/core/utils/qtl_helpers.py` (also used by the pQTL builder).
- **Source constants:** `EQTL_SOURCES = ("gtex_v11", "gtex_v8", "eqtlgen")` in `hvantk/skills/gtex_eqtl/shared/constants.py`.
- **Plugin resolution:** declared in `hvantk/skills/gtex_eqtl/plugin.yaml` (dataset `eqtls`, builder `hvantk.skills.gtex_eqtl.builder:build_eqtl_associations`). The plugin loader (`hvantk/core/plugin/loader.py`) auto-resolves `gtex-eqtl:eqtls` via `get_registry().get_dataset(...)`; the build runs through `run_builder_for_spec` (`hvantk/core/plugin/run_builder.py`). No `TABLE_BUILDERS` registry or `create_table_adapter`.
- **CLI:** `hvantk reprocess gtex-eqtl:eqtls --raw-dir <dir> --output <path>.ht`. Source-specific kwargs (`source`, `tissue`, `p_threshold`, `reference_genome`) flow through `--plugin-arg key=value`.
- **Downstream consumer:** `hvantk/qtlcascade/` — the eQTL Hail Table is one half of the eQTL ⊕ pQTL cascade join.

## 7. Workflow steps

1. **Resolve raw path.** Caller passes either a single parquet file (`Liver.v11.eQTLs.signif_pairs.parquet`) or a directory containing per-tissue parquet files. Acquire from the GTEx portal (https://gtexportal.org) — manual download, no skill-side acquisition.
2. **Import (v11 path).** `_import_gtex_parquet` scans `.parquet` files under `parsed_input`, optionally filters by `tissue`, reads each via `spark.read.parquet`, converts via `hl.Table.from_spark`, and unions the per-tissue tables. Source columns selected: `phenotype_id`, `variant_id`, `slope`, `slope_se`, `pval_nominal`, `af` (yields output `af`; `maf` is derived in the same select via `hl.min(af, 1-af)`).
3. **Transform** (inline in `build_eqtl_associations`):
   - `parse_gtex_variant_id(ht, "variant_id", reference_genome)` — yields `locus`, `alleles`.
   - `gene_id = strip_ensembl_version(ht.gene_id_raw)` — strips `.N` version suffix.
   - `ht.drop("gene_id_raw", "variant_id")`.
   - If `p_threshold > 0`: `ht.filter(ht.p_value <= p_threshold)`. **Set `p_threshold=0` to retain all rows** (signif_pairs files are already pre-filtered to significant pairs by GTEx; this is what coloc / qtlcascade typically want).
   - `ht.annotate(source=source, is_cis=True)`.
   - `ht.key_by("locus", "alleles", "gene_id")`.
4. **Wrap + provenance.** `AnnotationTable.from_hail(ht, provenance=ctx.provenance(schema_id="gtex-eqtl-eqtls-v1"))`. The `hvantk reprocess` runner checkpoints the artifact to `--output` and writes the provenance sidecar.

**Default usage (CLI, single tissue):**

```bash
hvantk reprocess gtex-eqtl:eqtls \
    --raw-dir /path/to/gtex_v11_eqtl/ \
    --output /path/to/gtex_v11_eqtl_Liver.ht \
    --skip-download \
    --plugin-arg source=gtex_v11 \
    --plugin-arg tissue=Liver \
    --plugin-arg p_threshold=0 \
    --plugin-arg reference_genome=GRCh38
```

## 8. Update playbook

GTEx releases major versions every few years (v8 → v9 → v10 → v11). Per release:

1. **Acquire** the new release's per-tissue signif_pairs files (manual download).
2. **Check format.** Compare the new release's `README_eQTL_*.txt` against the v11 schema documented in this skill. If columns added/removed: extend `_import_gtex_parquet` and re-record the snapshot.
3. **Add a new source constant** to `EQTL_SOURCES` in `hvantk/skills/gtex_eqtl/shared/constants.py` (e.g., `"gtex_v12"`); add a `gtex_v12` dispatch branch in `build_eqtl_associations` if the import shape diverges from v11; bump the registry entry's `accession` (`GTEx_v11_eQTL_signif_pairs` → `GTEx_v12_eQTL_signif_pairs`).
4. **Re-run round-trip (§ 9).** If schema changes, regenerate snapshots via `--regenerate-snapshots`.
5. **`af` / `maf`.** v11 reads `af` directly and derives `maf = min(af, 1-af)`. If a future release renames `af` → something else, the `af_expr` probe in `_import_gtex_parquet` falls back to `maf` (yielding null direction info, matching v8 behavior). Update the probe order if needed.

**Cross-source compatibility:** v8 TSV → v11 parquet was a breaking format change (text → columnar binary). v11 → future releases may break again. Anchor each release to its own source constant rather than reusing `gtex_v11` after the next migration.

## 9. Validation contract

Per `_conventions` § 9:

- **fixture:** `hvantk/skills/gtex_eqtl/tests/testdata/raw/gtex-eqtl/Liver.v11.eQTLs.signif_pairs.parquet`. 800 rows sliced from the v11 Liver source via `pyarrow.parquet.read_table(...).slice(0, 800)` and `pq.write_table` — deterministic because row-group order is preserved. ~36 KB. Keeps the parquet binary format so the test exercises the real `spark.read.parquet → hl.Table.from_spark` ingestion path.
- **schema_snapshot:** `hvantk/skills/gtex_eqtl/tests/snapshots/schema.json`.
- **row_snapshot:** `hvantk/skills/gtex_eqtl/tests/snapshots/sample_rows.json`. Triple key `(locus, alleles, gene_id)` is unique-in-table for signif_pairs, so the test inlines a small key list and no `sample_keys.json` is maintained (per `_conventions` § 9 post-#101).
- **drift_fingerprint:** `hvantk/skills/gtex_eqtl/tests/drift_fingerprint.json`.
- **test_command:** `pytest hvantk/skills/gtex_eqtl/tests -m hail`.

Round-trip test asserts: checkpointed schema matches `schema.json`; deterministic sample-row slice matches `sample_rows.json`. Test uses `source="gtex_v11"`, `tissue="Liver"`, `p_threshold=0` (retain all rows in the fixture). Spark session is initialized on demand by `_import_gtex_parquet`; the `hail_session` fixture ensures Hail (and therefore Spark) is up.

Regenerate via `--regenerate-snapshots` when:
- The shared transform in `build_eqtl_associations` adds / removes / renames fields.
- A new GTEx release is anchored that changes the input schema (per § 8).
