---
name: hvantk:resource-insider
description: Build a Hail Table from the Interactome Insider genomic BED (protein-protein interface residues projected to GRCh38) for interval-based variant annotation.
status: provisional
backend: hail
domain: protein
---

# INSIDER (Interactome Insider)

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes every convention there.

## 1. Status & scope

- **Status:** provisional. The builder, CLI, registry adapter, and catalog entry pre-exist this skill (the catalog filename is corrected in the same PR — see § 4 Gap 2). This skill is the design contract.
- **Anchored variant:** `Whole_Human_Interactome_Interface_hg38.bed` — the genomic projection product. UCSC-style BED with browser/track metadata; 18.6M data rows across **208,448 named PPI tracks**.
- **In scope:** building an `interval`-keyed Hail Table from the BED file for variant-interval intersection (e.g., "does this variant fall in any predicted interface residue?").
- **Out of scope:**
  - The complementary `H_sapiens_interfacesALL.txt` product (protein-pair keyed; encodes per-protein interface residue arrays with Source = ECLAIR / PDB / I3D). This is a separately-onboardable skill — same resource, different shape and different builder. **Not in this PR.**
  - Downloader. The BED is >1 GB; per `_conventions` § 11 and the downloader strategy (CLAUDE.md), acquisition is manual.
  - Per-PPI variant annotation — see Gap 1 in § 4. The current builder cannot answer "which PPI did this variant intersect?" because track identity is dropped during ingestion.

This skill is the **first interval-keyed skill** in hvantk. Conventions § 3 declares `interval` keying as valid; this skill anchors it.

## 2. Source identity

- **Provider:** Yu lab (Cornell). Wei et al., *Nat Methods* 2017, PMID 29036289.
- **Distribution:** http://interactomeinsider.yulab.org/downloads.html
- **License:** Academic use (per the existing catalog entry).
- **Catalog entry:** `INSIDER_v1.0` in `hvantk/resources/registry/genomics/datasets.json`. **Filename and metadata corrected in the same PR that adds this skill** — the prior entry listed `insider_interaction_sites.tsv` which is not a real INSIDER distribution product (see § 4 Gap 2).

Stable note (not in catalog): INSIDER releases two complementary products from one source. This skill anchors the genomic BED only. The protein-residue TXT is documented in the catalog as a follow-up; see § 8.

## 3. Backend choice + reasoning

**`backend: hail`, `domain: protein`.** Per `_conventions` § 3, protein-level annotations key on `interval` or `protein_id`. The genomic projection BED is naturally interval-keyed: each row is a genomic span (typically 1-3 bp) corresponding to a codon at a predicted interface residue. Hail's `hl.import_bed` produces an `interval<locus<rg>>`-keyed Table directly.

Ingestion path: `hl.import_bed(path, skip_invalid_intervals=True, reference_genome=...)` → `.distinct()`. This skill is the first to anchor the `hl.import_bed` path; conventions § 4 lists it implicitly via `_create_table_base().import_func` accepting any `Callable[[], hl.Table]`.

**Snapshot utility upgrade.** The `_snapshot_utils._jsonable_to_hail_python` helper did not previously handle `hl.tinterval`. This PR adds an interval branch so interval-keyed snapshots round-trip cleanly. Format: `{"start": "<contig>:<pos>", "end": "<contig>:<pos>"}` — symmetric with locus serialization, half-open `[start, end)` per the BED convention. See `hvantk/tests/_snapshot_utils.py:70-82` (the new branch).

## 4. Raw format & gotchas

UCSC-style BED with browser/track metadata. Three structural lines to be aware of:

```
browser hide all
track name=A0A0A0MS80_ppi_P56705 description="Interface Residues for A0A0A0MS80_P56705 (SOURCE: Predicted Interface)" visibility=dense itemRgb="On"
chr11   700235   700235        .   0   +   700235   700235   247,176,91
```

- The `browser` directive (line 1) is ignored by `hl.import_bed`.
- `track name=<P1>_ppi_<P2> description="..."` lines name each PPI. **The track header itself is silently skipped by `hl.import_bed`.** Information loss — see Gap 1.
- BED data rows: `chr, start, end, name, score, strand, thickStart, thickEnd, itemRgb`. The `name` field is `.` for all rows in this file; the `score`, `thick*`, and `itemRgb` fields are presentation-only (UCSC track visualization).
- **Zero-length intervals are common**: many rows have `start == end` (e.g., `chr11 700235 700235`). Hail's `hl.import_bed(skip_invalid_intervals=True)` filters these. In the test fixture, 21 raw data rows produced 17 Hail rows (4 zero-length skipped).

**Gap 1 (real builder limitation, deferred): PPI track identity is dropped.** `hl.import_bed` ignores the `track name=...` metadata lines. After ingestion, each row's `target` field is `.` (the BED name column) — no field carries the PPI identity. Then `_create_interactome_tb` calls `.distinct()` which collapses overlapping intervals from different PPIs into one row. Net effect: the Table answers "does any PPI interface touch this position?" but **not** "which PPI(s)?". Real information loss. Follow-up PR would need to parse track metadata manually (split the BED on `track` lines, assign a `ppi_id: str` column per track, then union, then optionally aggregate `array<str>` of PPI IDs per interval instead of `.distinct()`). Cite: `hvantk/tables/table_builders.py:238-243`.

**Gap 2 (fixed inline in this PR): wrong filename in the catalog entry.** Prior `INSIDER_v1.0` catalog entry listed `insider_interaction_sites.tsv` which does not exist in any INSIDER distribution. The two real products are `Whole_Human_Interactome_Interface_hg38.bed` (this skill) and `H_sapiens_interfacesALL.txt` (separate skill). The filename, format, and size_bytes are corrected in the same PR; the entry is rewritten to describe the actual content. Cite: `hvantk/resources/registry/genomics/datasets.json` (entry `INSIDER_v1.0`).

**Gap 3 (documented, not a bug): >1 GB file → manual acquisition.** Acknowledged by the downloader strategy. No skill-side downloader.

## 5. Output contract

- **Object:** `hl.Table` checkpointed to `output_path` (a `.ht` directory).
- **Key:** `[interval]` (`interval<locus<GRCh38>>`).
- **Globals:** `hvantk_metadata` set by `_create_table_base`.
- **Fields:**
  - `interval: interval<locus<GRCh38>>` — half-open `[start, end)` per BED.
  - `target: str` — BED column 4 (`name`). For this file, always `"."` (the data carries no per-row name; PPI identity is in the track header, which is dropped — see Gap 1).
- **Reference genome:** GRCh38. Required argument to `create_interactome_tb`.
- **Snapshot key form:** intervals serialize as `{"start": "<contig>:<pos>", "end": "<contig>:<pos>"}` (the new branch in `_snapshot_utils`).

After `.distinct()`, intervals are unique-in-table. Test inlines sample keys.

## 6. hvantk integration points

- **Builder:** `create_interactome_tb` in `hvantk/tables/table_builders.py:200`. Uses `_create_table_base()` with `import_func = lambda: hl.import_bed(...)` and `transform_func = lambda ht: ht.repartition(100).distinct()`.
- **Registry:** `TABLE_BUILDERS["interactome"] = create_table_adapter("hvantk.tables.table_builders", "create_interactome_tb")` in `hvantk/tables/registry.py:250`.
- **CLI:** `mktable_interactome` in `hvantk/commands/make_table_cli.py:168` (command name `interactome`). Standard input/output/overwrite/export options.
- **Snapshot util branch:** `hvantk/tests/_snapshot_utils.py` — new `hl.tinterval` handlers added in this PR.
- **Downloader:** out of scope (manual acquisition).

## 7. Workflow steps

1. **Acquire** the BED file from http://interactomeinsider.yulab.org/downloads.html — manual download, no skill-side acquisition (file is >1 GB).
2. **Build the Hail Table:**
   ```bash
   hvantk mktable interactome \
       --raw-input /path/to/Whole_Human_Interactome_Interface_hg38.bed \
       --output-ht /path/to/insider.ht \
       --ref-genome GRCh38
   ```
3. **Internal flow** (already implemented in `create_interactome_tb`):
   - `hl.import_bed(input_path, skip_invalid_intervals=True, reference_genome="GRCh38")` — parses BED, ignoring browser/track metadata. Outputs key `interval`, field `target`.
   - `.repartition(100).distinct()` — collapses duplicate intervals.
   - `_create_table_base` handles checkpoint, globals, optional TSV export.
4. **Use for variant annotation:** semi-join a variant Table on `interval.contains(variant.locus)` or equivalent interval-overlap operator (`hl.is_defined(insider_ht[variant.locus])` after appropriate keying).

## 8. Update playbook

INSIDER updates are irregular. To onboard a new release:

1. **Acquire** the new BED. Update `path` / `size_bytes` / `last_updated` in the catalog entry; bump `accession` if the release version changes (`INSIDER_v1.0` → `INSIDER_v1.x`).
2. **Re-run round-trip (§ 9).** If the BED format is unchanged (still 9-column UCSC-style with `track` headers), no builder change.
3. **Resolve Gap 1 (track identity loss)** — follow-up PR. The skill anchor will need rewriting if track-id preservation changes the output schema (likely adds `ppi_id: array<str>` and re-keys).
4. **Consider onboarding the `.txt` product** as a sibling skill. It carries Source provenance (ECLAIR / PDB / I3D) which the BED does not; downstream filtering on confidence level requires the TXT. Builder would need range-notation parsing for `*_IRES` arrays (e.g., `[1-11,13-14,...]`).

## 9. Validation contract

Per `_conventions` § 9:

- **fixture:** `hvantk/tests/testdata/raw/insider/insider_sample.bed`. 5 PPI tracks (~21 raw data rows; 17 valid after `skip_invalid_intervals`). ~1.9 KB. Sliced via `local/planning/skills-insider-fixture-slicer.py` (gitignored) which preserves track headers (a custom BED-aware slicer rather than `head -N` to maintain structural integrity).
- **schema_snapshot:** `hvantk/tests/snapshots/insider/schema.json`.
- **row_snapshot:** `hvantk/tests/snapshots/insider/sample_rows.json`. Intervals are unique-in-table after `.distinct()`; test inlines 3-4 sample keys (per `_conventions` § 9 post-#101 rule — unique-key skills inline).
- **test_command:** `pytest hvantk/tests/test_insider_builder.py -m hail`.

Round-trip test asserts: builder idempotent with `overwrite=True`; checkpointed schema matches `schema.json`; deterministic sample-row slice matches `sample_rows.json`. The test exercises the new `hl.tinterval` handling in `_snapshot_utils` — if that branch breaks, this test breaks.

Regenerate via `--regenerate-snapshots` when:
- Gap 1 (track identity) is fixed — the output schema gains `ppi_id` and the row snapshot rebases.
- A new INSIDER release changes the BED column shape (currently unchanged for v1.0).
