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

- **Status:** provisional. The builder and catalog entry pre-exist this skill (the catalog filename is corrected in the same PR — see § 4 Gap 2). This skill is the design contract.
- **Anchored variant:** `Whole_Human_Interactome_Interface_hg38.bed` — the genomic projection product. UCSC-style BED with browser/track metadata; 18.6M data rows across **208,448 named PPI tracks**.
- **In scope:** building an `interval`-keyed Hail Table from the BED file for variant-interval intersection (e.g., "does this variant fall in any predicted interface residue?").
- **Out of scope:**
  - The complementary `H_sapiens_interfacesALL.txt` product (protein-pair keyed; encodes per-protein interface residue arrays with Source = ECLAIR / PDB / I3D). This is a separately-onboardable skill — same resource, different shape and different builder. **Not in this PR.**
  - Downloader. The BED is >1 GB; per `_conventions` § 11 and the downloader strategy (CLAUDE.md), acquisition is manual.

This skill is the **first interval-keyed skill** in hvantk. Conventions § 3 declares `interval` keying as valid; this skill anchors it.

## 2. Source identity

- **Provider:** Yu lab (Cornell). Wei et al., *Nat Methods* 2017, PMID 29036289.
- **Distribution:** http://interactomeinsider.yulab.org/downloads.html
- **License:** Academic use (per the existing catalog entry).
- **Catalog entry:** `INSIDER_v1.0` in `hvantk/resources/registry/genomics/datasets.json`. **Filename and metadata corrected in the same PR that adds this skill** — the prior entry listed `insider_interaction_sites.tsv` which is not a real INSIDER distribution product (see § 4 Gap 2).

Stable note (not in catalog): INSIDER releases two complementary products from one source. This skill anchors the genomic BED only. The protein-residue TXT is documented in the catalog as a follow-up; see § 8.

## 3. Backend choice + reasoning

**`backend: hail`, `domain: protein`.** Per `_conventions` § 3, protein-level annotations key on `interval` or `protein_id`. The genomic projection BED is naturally interval-keyed: each row is a genomic span (typically 1-3 bp) corresponding to a codon at a predicted interface residue. The builder produces an `interval<locus<rg>>`-keyed Table with a `ppi_ids: array<str>` field per interval, aggregating PPI identifiers across the `track name=...` directives.

**Ingestion path: custom track-aware parser**, not `hl.import_bed` directly. `hl.import_bed` silently skips `track name=...` directives, dropping PPI identity (see § 4 historical note on the prior implementation). This builder instead pre-processes the BED line-by-line in Python: it tracks the current PPI from each `track name=<P1>_ppi_<P2>` directive, writes a 4-column TSV (`contig\tstart\tend\tppi_id`) to a Hail temp file, then imports via `hl.import_table`, constructs intervals with `hl.locus_interval(contig, start+1, end+1, ...)` to match `hl.import_bed`'s 0-based-BED → 1-based-Hail conversion, and aggregates `group_by(interval).aggregate(ppi_ids=hl.agg.collect_as_set(ppi_id))`. The `ppi_ids` array is then `hl.sorted(hl.array(...))` for deterministic snapshots.

**Snapshot utility upgrade.** The `_snapshot_utils._jsonable_to_hail_python` helper did not previously handle `hl.tinterval`. PR #105 added an interval branch so interval-keyed snapshots round-trip cleanly. Format: `{"start": "<contig>:<pos>", "end": "<contig>:<pos>"}` — symmetric with locus serialization, half-open `[start, end)` per the BED convention. See `hvantk/tests/_snapshot_utils.py` (the interval branch).

## 4. Raw format & gotchas

UCSC-style BED with browser/track metadata. Three structural lines to be aware of:

```
browser hide all
track name=A0A0A0MS80_ppi_P56705 description="Interface Residues for A0A0A0MS80_P56705 (SOURCE: Predicted Interface)" visibility=dense itemRgb="On"
chr11   700235   700235        .   0   +   700235   700235   247,176,91
```

- The `browser` directive (line 1) is ignored.
- `track name=<P1>_ppi_<P2> description="..."` lines name each PPI. **The builder parses these and assigns the parsed name as the row's `ppi_id` for the subsequent data block**, preserving PPI identity (this is the fix to the historical Gap 1; see below).
- BED data rows: `chr, start, end, name, score, strand, thickStart, thickEnd, itemRgb`. Only the first three columns are read; column 4 (`name`) is always `.` in this file (PPI identity is in the track header instead).
- **Zero-length intervals are common**: many rows have `start == end` (e.g., `chr11 700235 700235`). The custom parser skips these (matching `hl.import_bed(skip_invalid_intervals=True)`'s historical behavior). In the test fixture, 21 raw data rows produced 17 valid intervals (4 zero-length skipped); the aggregation then collapses overlapping intervals across PPIs, yielding 17 unique-interval rows here (each with a singleton `ppi_ids` since the 5-track fixture doesn't include cross-PPI overlaps).
- **Multi-PPI intervals.** In the full 208,448-track file, a genomic position can fall on the interface of multiple PPIs (a residue in a hub protein that participates in many complexes). The aggregation `group_by(interval).aggregate(ppi_ids=collect_as_set(ppi_id))` produces a length-N array per such position. The fixture doesn't exercise this case (each interval has a singleton array), so the multi-PPI behavior is by inspection only — not covered by the round-trip test.

**Historical note (Gap 1, fixed in PR #105 via the track-aware parser):** prior implementation used `hl.import_bed(...).distinct()`, which silently dropped `track name=...` headers and then collapsed overlapping intervals from different PPIs. The output table answered "does any PPI interface touch this position?" but **not** "which PPI(s)?". The current builder restores that identity via the custom parser described in § 3 (now `_parse_insider_bed_to_temp_tsv` in `hvantk/skills/insider/builder.py`).

**Gap 2 (fixed in PR #105): wrong filename in the catalog entry.** Prior `INSIDER_v1.0` catalog entry listed `insider_interaction_sites.tsv` which does not exist in any INSIDER distribution. The two real products are `Whole_Human_Interactome_Interface_hg38.bed` (this skill) and `H_sapiens_interfacesALL.txt` (separate skill). The filename, format, and size_bytes were corrected in PR #105.

**Gap 3 (documented, not a bug): >1 GB file → manual acquisition.** Acknowledged by the downloader strategy. No skill-side downloader.

## 5. Output contract

- **Object:** an `AnnotationTable` artifact wrapping a lazy `hl.Table` (materialized to a `.ht` directory on `artifact.save()`).
- **Key:** `[interval]` (`interval<locus<GRCh38>>`).
- **Provenance:** stamped via `ctx.provenance(schema_id="insider-variants-v1")`; persisted as a sidecar `.provenance.json`.
- **Fields:**
  - `interval: interval<locus<GRCh38>>` — half-open `[start, end)`. Constructed via `hl.locus_interval(contig, start+1, end+1, ...)` to match `hl.import_bed`'s 0-based-BED → 1-based-Hail conversion, so semantic compatibility with prior interval-based variant annotation downstream is preserved.
  - `ppi_ids: array<str>` — sorted, deduplicated PPI identifiers (`<P1>_ppi_<P2>` format) from `track name=...` directives whose data rows cover this interval. Length ≥ 1; arrays of length > 1 indicate the position is shared between multiple PPI interfaces.
- **Reference genome:** GRCh38. Passed as the `reference_genome` keyword param to `build_insider_interactome` (default `"GRCh38"`).
- **Snapshot key form:** intervals serialize as `{"start": "<contig>:<pos>", "end": "<contig>:<pos>"}` (the branch in `_snapshot_utils` added in PR #105).

After aggregation, intervals are unique-in-table. Test inlines sample keys.

## 6. hvantk integration points

- **Builder:** `build_insider_interactome(parsed_input, ctx, *, reference_genome="GRCh38") -> AnnotationTable` in `hvantk/skills/insider/builder.py`. Builds the table inline: calls `_parse_insider_bed_to_temp_tsv(str(parsed_input))` (track-aware Python pre-processor), then `hl.import_table` + `hl.locus_interval`, then `group_by(interval).aggregate(ppi_ids=collect_as_set(ppi_id))` + `hl.sorted(hl.array(...))`, and wraps the result via `AnnotationTable.from_hail(grouped, provenance=ctx.provenance(schema_id="insider-variants-v1"))`. The only shared helper it imports is `cleanup_temp_file` from `hvantk/core/utils/hail_helpers.py` (used on the error path). It does **not** use `create_table_base`.
- **Track parser helper:** `_parse_insider_bed_to_temp_tsv` (private) lives in `hvantk/skills/insider/builder.py` (only consumer is this plugin). Reads the BED, tracks `current_ppi_id` from `track name=...` headers, skips zero-length and malformed rows, writes a 4-column TSV to `hl.utils.new_temp_file(extension="tsv")`.
- **Registry / loader:** plugin-driven; the in-tree plugin manifest at `hvantk/skills/insider/plugin.yaml` declares dataset `variants` with `builder.function: build_insider_interactome`. The plugin loader (`hvantk/core/plugin/loader.py`) auto-resolves the dataset under compound key `insider:variants` via `get_registry().get_dataset("insider:variants")`; the top-level build runs through `run_builder_for_spec` (`hvantk/core/plugin/run_builder.py`). There is no `TABLE_BUILDERS` registry or `registry.py`.
- **CLI:** `hvantk reprocess insider:variants --raw-dir <dir> --output <out> [--plugin-arg reference_genome=GRCh38]`. `reprocess` passes the `--raw-dir` directory to the builder as `parsed_input` (insider has no `lifecycle.parse`); `_parse_insider_bed_to_temp_tsv` consumes it as the BED path.
- **Snapshot util branch:** `hvantk/tests/_snapshot_utils.py` — `hl.tinterval` handlers added in PR #105.
- **Downloader:** out of scope (manual acquisition; the BED is >1 GB, so no `lifecycle.download` entry in `plugin.yaml`).

## 7. Workflow steps

1. **Acquire** the BED file from http://interactomeinsider.yulab.org/downloads.html — manual download, no skill-side acquisition (file is >1 GB).
2. **Build the Hail Table** via the reprocess CLI (point `--raw-dir` at the directory containing the BED):
   ```bash
   hvantk reprocess insider:variants \
       --raw-dir /path/to/insider_raw_dir \
       --output /path/to/insider.ht \
       --plugin-arg reference_genome=GRCh38
   ```
3. **Internal flow** (implemented in `build_insider_interactome`):
   - `_parse_insider_bed_to_temp_tsv(input_path)` — Python-side BED iterator: tracks current PPI from `track name=...`, skips browser lines, skips zero-length and malformed rows, writes `contig\tstart\tend\tppi_id` to a Hail temp TSV.
   - `hl.import_table(tmp_tsv, types={"start": tint32, "end": tint32})` → annotate with `hl.locus_interval(contig, start+1, end+1, reference_genome=...)`. The `+1` matches `hl.import_bed`'s 0-based-BED → 1-based-Hail conversion so intervals are consistent with prior `hl.import_bed`-derived outputs.
   - `group_by(interval).aggregate(ppi_ids=hl.agg.collect_as_set(ppi_id))` → `grouped.annotate(ppi_ids=hl.sorted(hl.array(ppi_ids)))` — collapse overlapping intervals while preserving the union of PPI IDs, sorted for deterministic output.
   - `AnnotationTable.from_hail(grouped, provenance=ctx.provenance(schema_id="insider-variants-v1"))` wraps the lazy Table; checkpointing happens on `artifact.save()`.
4. **Use for variant annotation:** semi-join a variant Table on `interval.contains(variant.locus)` or equivalent interval-overlap operator (`hl.is_defined(insider_ht[variant.locus])` after appropriate keying). Downstream consumers can filter / explode on `ppi_ids` to attribute hits to specific PPIs.

## 8. Update playbook

INSIDER updates are irregular. To onboard a new release:

1. **Acquire** the new BED. Update `path` / `size_bytes` / `last_updated` in the catalog entry; bump `accession` if the release version changes (`INSIDER_v1.0` → `INSIDER_v1.x`).
2. **Re-run round-trip (§ 9).** If the BED format is unchanged (still 9-column UCSC-style with `track name=...` directives in the established `<P1>_ppi_<P2>` shape), no builder change. If a release changes the track naming pattern (e.g., adds a third underscore-separated field), the `_TRACK_NAME_RE` regex in `hvantk/skills/insider/builder.py` and the `<P1>_ppi_<P2>` convention in the catalog description need updating.
3. **Consider onboarding the `.txt` product** as a sibling skill. It carries Source provenance (ECLAIR / PDB / I3D) which the BED does not; downstream filtering on confidence level requires the TXT. Builder would need range-notation parsing for `*_IRES` arrays (e.g., `[1-11,13-14,...]`).

## 9. Validation contract

Per `_conventions` § 9:

- **fixture:** `hvantk/skills/insider/tests/testdata/raw/insider/insider_sample.bed`. 5 PPI tracks (~21 raw data rows; 17 valid after zero-length filtering). ~1.9 KB. Sliced from the full BED by a track-aware sub-sampler (keeps the `browser` directive plus the first N `track` blocks, each header paired with its data rows) — a plain `head -N` would split a track block and produce an invalid BED.
- **schema_snapshot:** `hvantk/skills/insider/tests/snapshots/schema.json`. Records the `{interval, ppi_ids: array<str>}` shape.
- **row_snapshot:** `hvantk/skills/insider/tests/snapshots/sample_rows.json`. Intervals are unique-in-table after the aggregation; test inlines 3 sample keys (per `_conventions` § 9 post-#101 rule — unique-key skills inline).
- **test_command:** `pytest hvantk/skills/insider/tests -m hail`.

Round-trip test (`hvantk/skills/insider/tests/test_builder.py`, via `phase_b_snapshot_adapter(build_insider_interactome, "insider:variants")`) asserts: checkpointed schema matches `schema.json`; deterministic sample-row slice matches `sample_rows.json`. The test exercises the `hl.tinterval` handling in `_snapshot_utils` (added in PR #105) — if that branch breaks, this test breaks.

Regenerate via `--regenerate-snapshots` when:
- The builder adds / removes / renames fields (e.g., if a future PR adds a `source: array<str>` derived from `track description="..."` Source values).
- A new INSIDER release changes the BED column shape (currently unchanged for v1.0).
