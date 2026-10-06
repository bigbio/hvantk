---
name: hvantk:resource-alphagenome
description: AlphaGenome variant-effect scores, ingested from the SDK's tidy_scores() parquet (produced by the user through the API or local weights) and summarised into one Hail Table row per variant.
status: provisional
backend: hail
domain: genomics
---

# alphagenome

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

Provisional. Covers the `alphagenome:predictions` builder, which **only ingests** AlphaGenome scores the user has already produced: a directory of parquet files written from the AlphaGenome SDK's `variant_scorers.tidy_scores()` long format. hvantk makes **no AlphaGenome API calls** and needs no credentials; it does not import the AlphaGenome SDK at all.

The builder summarises the long format into one row per variant (that still has rows after the filters), keyed by `(locus, alleles)`, with one struct per known variant scorer (§ 5).

Out of scope:
- Producing the scores. That is the user's step, with the SDK (§ 2), under the AlphaGenome Output Terms. The manifest declares `acquisition: {mode: byo, reason: credentialed}`, so there is no `lifecycle.download`.
- Track-level output. Only the per-variant, per-scorer summary is kept; the per-track rows stay in the user's parquet.
- Interval scorers (`tidy_scores()` output with an `interval_scorer` column and no `variant_id`).

## 2. Source identity

AlphaGenome is a DNA sequence model from Google DeepMind that predicts variant effects on gene expression, splicing, chromatin and contact maps. Its SDK is the `alphagenome` package on PyPI (docs: https://www.alphagenomedocs.com). `source.catalog_ref: alphagenome` is declared in `plugin.yaml`.

**Licence: non-commercial.** AlphaGenome outputs are subject to the AlphaGenome Output Terms of Use (https://deepmind.google.com/science/alphagenome/output-terms): non-commercial use only, never used to train machine-learning models, and redistributed only with a conspicuous notice. The committed fixture carries that notice in `tests/testdata/raw/alphagenome/NOTICE.md`; keep it beside the parquet. Scores a user builds into a Table stay under the same terms.

**What produces the input.** The SDK's `score_variant()` returns one AnnData per scorer; `variant_scorers.tidy_scores()` turns them into the long DataFrame hvantk reads. The model can be the hosted API (`dna_client.create(api_key)`) or local weights (`alphagenome_research.model.dna_model.create_from_kaggle("all_folds")`); both expose the same `score_variant`, so the rest is identical. The calls below are taken from the SDK documentation's batch-scoring notebook and match the run that produced the fixture's source shard:

```python
import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

model = dna_client.create(api_key)
scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())  # the 19 hvantk maps
cap = dna_client.MAX_VARIANT_SCORERS_PER_REQUEST

frames = []
for variant in variants:  # e.g. genome.Variant(chromosome="chr3", position=39408741,
    #                                           reference_bases="T", alternate_bases="C")
    interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    scores = []
    for i in range(0, len(scorers), cap):
        scores.extend(
            model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=scorers[i : i + cap],
                organism=dna_client.Organism.HOMO_SAPIENS,
            )
        )
    frames.append(variant_scorers.tidy_scores(scores))

df = pd.concat(frames, ignore_index=True)
# tidy_scores() stores the genome.Variant / genome.Interval objects themselves;
# parquet needs strings. str(Variant) is "chr3:39408741:T>C".
for column in ("variant_id", "scored_interval"):
    df[column] = df[column].astype(str)
df.to_parquet("scores.parquet", index=False)
```

## 3. Backend choice + reasoning

`hail` (`plugin.yaml`: `backend: hail`). The output is a variant-level annotation keyed by `(locus, alleles)` (`_conventions` § 3), joinable to every other variant table in the toolkit. The input is large: the real shard behind the fixture holds 16,642,316 rows for 411 variants (about 40,000 rows per variant), so the aggregation runs in Hail/Spark rather than pandas. On 8 cores that shard builds in about 20 seconds.

The parquet is read with Spark (`SparkSession.builder.getOrCreate()` → `spark.read.parquet(...)` → `hl.Table.from_spark`), the same bridge `gtex_eqtl` uses. The builder calls `init_hail()` first, so the Spark session is Hail's. Spark also checks and filters the rows in one pass before the conversion, so a bad input fails before Hail's shuffle; Hail only aggregates.

## 4. Raw format & gotchas

- **Long format.** One row per variant × scorer × track, and × gene (gene-level scorers) or × junction (`SpliceJunctionScorer()`). The builder reads `variant_id`, `output_type`, `variant_scorer`, `track_name`, `track_strand`, `gene_id`, `raw_score`, `quantile_score` (all required in every file; the SDK always writes `gene_id` and `track_strand`), `junction_Start`/`junction_End` when present, and `ontology_curie` (only for the `ontology_curies` filter). Every other column is ignored, including `Assay title` (with a space), `scored_interval`, `gtex_tissue`, and an `input_variant_id` column some producers add (identical to `variant_id` in the fixture's source). A required column missing from any input file raises `ValueError` naming the file and the columns.
- **Only `*.parquet` files directly in `--raw-dir` are read**, so `NOTICE.md` or a script beside them is ignored. Every one of them must be a tidy-scores file: a derived parquet left in the same folder (say, a per-variant summary written next to the shards) fails the build with that `ValueError`. Spark gets absolute paths, as it requires.
- **Each file is read with its own schema**, checked, and cast to strings and float64 scores before the files are combined with `unionByName`. Spark reading several files at once applies one file's schema to all of them, which fails on a column typed differently in another file (parquet `null` vs string `gene_id`, float32 vs float64 scores) and reads a column the schema-giving file lacks as null everywhere. The price is one read per file, so many small files are slow: concatenate them into a few large ones first.
- **`variant_id` is `chrom:pos:ref>alt`**, e.g. `chr3:39408741:T>C`: the SDK's `str(genome.Variant)`, **1-based** (checked against hg38), on GRCh38 contigs with the `chr` prefix. The builder splits it into `locus` and `alleles` without renaming contigs. Any other form (e.g. `chr3:39408741:T:C`, bases outside ACGTN, a null) raises `ValueError` naming an example; the check mirrors the SDK's own default format, where an allele may be empty.
- **`output_type` is always read from its column, never parsed from `variant_scorer`.** Three scorer strings name no output type: `PolyadenylationScorer()` (RNA_SEQ), `ContactMapScorer()` (CONTACT_MAPS) and `SpliceJunctionScorer()` (SPLICE_JUNCTIONS). The `output_types` parameter filters on this column and only accepts the 11 values in `OUTPUT_TYPES`.
- **19 scorers, matched by exact string.** `variant_scorer` holds `str(scorer)`, e.g. `CenterMaskScorer(requested_output=ATAC, width=501, aggregation_type=ACTIVE_SUM)`. `SCORER_FIELDS` in `shared/constants.py` maps the 19 strings of the SDK's `RECOMMENDED_VARIANT_SCORERS` to output fields. An unknown string raises `ValueError` naming it instead of being dropped or guessed: a new SDK can add scorers or change their parameters (and so their strings). The check covers the whole input, before any filter, so a filtered build fails on it too.
- **Scores are float32** in the parquet. They are cast to float64, which is exact, so `max_abs_raw` equals a float32 maximum computed elsewhere bit for bit. The real shard has no null or NaN score. The builder treats NaN as missing (NaN is a value, not a null, so it is tested separately) and drops only a row with neither score: a row with just a raw score still counts in `n_rows` and `max_abs_raw` but cannot be the `top_*` pick, and a row with just a quantile score can be the pick, with `top_raw` missing. Every dropped row and every such half-scored row is counted in a warning.
- **`gene_id` is null on the 13 track-level scorers** (the 12 `CenterMaskScorer` variants and `ContactMapScorer()`) and set on the 6 gene-level ones, already stripped of its Ensembl version by `tidy_scores()`. So `top_gene_id` is always missing on track-level fields. A file with only track-level scorers stores `gene_id` with parquet type `null`; the per-file cast turns it into null strings.
- **`junction_Start`/`junction_End`** are set only on SPLICE_JUNCTIONS rows (int64, nullable).
- **One row per variant, scorer, track, strand, gene and junction** (`ROW_KEY` in `builder.py`). The same key twice means the same scores were passed twice, a file given twice or two scoring runs mixed, and would be counted twice, so the build fails with `ValueError` naming one duplicated row. On the real 16.6M-row shard every row has its own key (checked 2026-10-07); without `track_strand` 2.1M rows would collide, without the junction columns 0.37M, without `gene_id` 12.9M. Keep one copy of each scoring run in `--raw-dir`, and keep the junction columns in a parquet that holds SPLICE_JUNCTIONS rows.
- **Rows per variant.** Track-level counts are fixed per output type (both aggregations together): ATAC 334, CAGE 1,092, CHIP_HISTONE 2,232, CHIP_TF 3,234, CONTACT_MAPS 28, DNASE 610, PROCAP 24. Gene-level counts follow gene density: RNA_SEQ 3,564–110,880, SPLICE_JUNCTIONS 367–19,818, SPLICE_SITE_USAGE 367–1,101, SPLICE_SITES 2–6. `PolyadenylationScorer()` scored only 307 of the 411 variants, so `rna_seq_polyadenylation` is often missing.
- **Ties at the top |quantile_score| are common**, because quantiles saturate near ±1 (0.99999994). In the real shard 862 of 7,705 (variant, scorer) groups have two or more rows at the maximum, always with different tracks or values. `(track_name, gene_id)` does not identify a row either: it repeats within the CAGE, PROCAP, RNA_SEQ and SPLICE_JUNCTIONS scorers. The builder therefore orders by |quantile_score| descending, then `track_name`, `gene_id` (missing last), `quantile_score` and `raw_score` ascending. Rows that tie on all five print identically, so the result is deterministic.
- **`gtex_tissue` is the empty string on non-GTEx tracks, not missing data.** 85% of SPLICE_SITE_USAGE rows have it empty; all of them are ENCODE tracks (`data_source == "encode"`), the rest are GTEx.
- **Filtering on the two heart curies is not a GTEx-only "Heart" filter.** UBERON:0006631 also tags an ENCODE SPLICE_SITE_USAGE track, `usage_UBERON:0006631 total RNA-seq` (right atrium auricular region). So `ontology_curies=["UBERON:0006566", "UBERON:0006631"]` selects three tracks: the GTEx `Heart_Left_Ventricle` and `Heart_Atrial_Appendage` tracks and that ENCODE one. A GTEx-only feature, such as `ag_heart_ssu` in an independent per-variant ClinVar benchmark computed by the maintainer (`gtex_tissue` containing "Heart"), needs a filter on `gtex_tissue`, which the builder does not expose. The two differ for 319 of 339 ClinVar variants in the real shard, where the heart build's `splice_site_usage.max_abs_raw` is larger; it is never smaller.
- **An ontology filter also drops every track with no ontology term.** SPLICE_SITES tracks (`donor`, `acceptor`) carry none, so `splice_sites` is always missing from an `ontology_curies` build; on the real shard the heart build keeps 11 of the 19 fields.
- **`--plugin-arg` lists.** `--plugin-arg ontology_curies=UBERON:0006566,UBERON:0006631` arrives as a list; a single value arrives as a plain string, which the builder wraps.
- **Filters that leave nothing fail.** If no row survives (no score at all, or `output_types` / `ontology_curies` match nothing), the builder raises `ValueError` listing the row count after each filter, e.g. `2 rows, 2 with a raw or quantile score, 0 with output_types ['DNASE']`, instead of writing an empty table. A filter that keeps some rows logs the same counts at INFO.

## 5. Output contract

- **Object:** an `AnnotationTable` wrapping a Hail Table, from `AnnotationTable.from_hail(ht, provenance=ctx.provenance(schema_id="alphagenome-v2"))`. `hvantk reprocess` writes it to `--output` with its provenance sidecar.
- **Key:** `locus: locus<GRCh38>`, `alleles: array<str>`. One row per distinct `variant_id` that still has rows after the filters: a variant with no row of the requested `output_types`, or no track with the requested `ontology_curies`, is absent.
- **`variant_id: str`**, as in the input.
- **19 summary fields, one per scorer**, each `struct{top_raw: float64, top_quantile: float64, top_track: str, top_gene_id: str, max_abs_raw: float64, n_rows: int64}`:
  - `top_raw`, `top_quantile`, `top_track`, `top_gene_id`: the row with the largest |quantile_score|, sign kept, ties broken as in § 4, among the rows that have a quantile score (all four missing when none has one);
  - `max_abs_raw`: the largest |raw_score| over the rows that have a raw score;
  - `n_rows`: the rows aggregated, after the filters.
  - The struct is missing when the variant has no row for that scorer, or none survived the filters.
- **Field names** (`SCORER_FIELDS` in `shared/constants.py` holds the exact scorer strings):

  | Output type | Fields |
  | --- | --- |
  | ATAC, CAGE, CHIP_HISTONE, CHIP_TF, DNASE, PROCAP (`CenterMaskScorer`) | `<type>_active_sum`, `<type>_diff_log2_sum`, e.g. `atac_active_sum`, `chip_histone_diff_log2_sum` |
  | CONTACT_MAPS (`ContactMapScorer()`) | `contact_maps` |
  | RNA_SEQ (`GeneMaskActiveScorer`, `GeneMaskLFCScorer`, `PolyadenylationScorer()`) | `rna_seq_active`, `rna_seq_lfc`, `rna_seq_polyadenylation` |
  | SPLICE_JUNCTIONS (`SpliceJunctionScorer()`) | `splice_junctions` |
  | SPLICE_SITE_USAGE, SPLICE_SITES (`GeneMaskSplicingScorer`) | `splice_site_usage`, `splice_sites` |

- **Parameters** (`--plugin-arg`): `reference_genome` (default `GRCh38`), `output_types` (keep only rows of these `output_type` values; default all) and `ontology_curies` (keep only rows of tracks with these `ontology_curie` values, before aggregating; default all). `ontology_curies` matches the curie alone, so it cannot keep GTEx tracks of a tissue while dropping ENCODE ones; the two heart curies select an ENCODE track too (§ 4).
- **Schema snapshot:** `tests/snapshots/schema.json`. `schema_id` `alphagenome-v2` replaced `alphagenome-v1`, whose builder returned only the input variants.

## 6. hvantk integration points

- Plugin manifest: `hvantk/skills/alphagenome/plugin.yaml` (dataset `predictions`, `artifact_type: AnnotationTable`, `schema_id: alphagenome-v2`), resolved by `get_registry().get_dataset("alphagenome:predictions")`.
- Builder: `build_alphagenome_predictions(parsed_input, ctx, *, reference_genome="GRCh38", output_types=None, ontology_curies=None)` in `hvantk/skills/alphagenome/builder.py`; its module docstring states the input contract.
- Constants: `SCORER_FIELDS` and `OUTPUT_TYPES` in `hvantk/skills/alphagenome/shared/constants.py`.
- Drift probe: `fetch_fingerprint` in `hvantk/skills/alphagenome/drift_probe.py`, run by `hvantk drift alphagenome:predictions`. It fingerprints the SDK's PyPI release stream, the event that can change the `tidy_scores()` columns or the scorer strings.
- Tests: `hvantk/skills/alphagenome/tests/` -- `test_builder.py` (round trip on the real-subset fixture, `hail`-marked), `test_builder_guards.py` (every check that fails the build, missing-score handling, files of different float types, string-valued filters; `hail`-marked), `test_alphagenome.py` (registration), `test_drift_probe.py` (offline, `requests_mock`). Run with the `command` in § 9.
- Build with: `python -m hvantk reprocess alphagenome:predictions --raw-dir <dir-with-parquet> --output <out.ht>`. No `lifecycle` or `cli:` block is declared.

## 7. Workflow steps

1. **Produce the scores** with the AlphaGenome SDK (§ 2), through the API or local weights, and write `tidy_scores()` output to parquet with `variant_id` and `scored_interval` as strings. Several files are fine: one per batch or shard.
2. **Put only those parquet files in one directory.** Derived parquet files must live elsewhere (§ 4).
3. **Choose filters, if any.** `--plugin-arg output_types=SPLICE_SITES,SPLICE_SITE_USAGE` keeps only those output types; `--plugin-arg ontology_curies=UBERON:0006566,UBERON:0006631` keeps the tracks tagged with the two heart curies, GTEx and ENCODE alike. That is not a GTEx-only heart feature: the builder has no `gtex_tissue` filter (§ 4).
4. **Build.** `--skip-download` is not needed: the manifest declares `acquisition.mode: byo`.

```bash
python -m hvantk reprocess alphagenome:predictions --raw-dir <dir-with-parquet> --output <out.ht>
```

## 8. Update playbook

Triggered by a new AlphaGenome SDK release, which the drift probe detects (`hvantk drift alphagenome:predictions`, or the fortnightly drift workflow).

1. **Re-check the `tidy_scores()` columns** in the new SDK against `REQUIRED_COLUMNS` in `builder.py`, and against `ontology_curie` for the filter. A renamed column fails the build with `ValueError`; fix the builder, not the data.
2. **Re-check the scorer list**: compare `str(s)` for every `s` in `variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()` against the keys of `SCORER_FIELDS`. A new or reprinted scorer fails the build with `ValueError` naming it. Map it to a field and bump `schema_id` (`alphagenome-v3`) and the manifest `version`.
3. **Regenerate the fixture and snapshots** when the format changes: rerun `tests/testdata/raw/alphagenome/make_fixture.py` on a shard written by the new SDK (it takes the source shard as its argument; its docstring records the selection and the AlphaGenome licence constraints), update `NOTICE.md`, then rerun the round trip with `--regenerate-snapshots` and read the diff.
4. **Regenerate the drift fingerprint** with `hvantk drift --regenerate alphagenome:predictions` once the change is handled, not in the same PR as a behaviour change (`_conventions` § 12).

What the probe cannot see: a server-side model update shipped without an SDK release. No unauthenticated probe can detect it; record the SDK version and model used when producing scores.

## 9. Validation contract

Declared in `plugin.yaml`'s `tests:` block (paths relative to `hvantk/skills/alphagenome/`):

- `fixture`: `tests/testdata/raw/alphagenome` -- `clinvar-subset.parquet`, 150 real rows of the 16.6M-row shard for 3 ClinVar variants (`chr3:39408741:T>C`, `chr6:112216367:C>A`, `chrX:153694448:T>G`), all 19 scorers each. Per (variant, scorer) it keeps the top-|quantile_score| row under the builder's tie-break, the top-|raw_score| row and one random row, so every `top_*` and `max_abs_raw` equals its full-shard value; plus every heart SPLICE_SITE_USAGE row. `make_fixture.py` records the cut. Non-commercial data under the AlphaGenome Output Terms: `NOTICE.md` must stay beside it.
- `schema_snapshot`: `tests/snapshots/schema.json`
- `row_snapshot`: `tests/snapshots/sample_rows.json` -- all three variants (one row each, so the keys are inlined in `test_builder.py`).
- `drift_fingerprint`: `tests/drift_fingerprint.json`
- `command`: `pytest hvantk/skills/alphagenome/tests -m hail` -- selects the `hail`-marked round trip and guard tests; the registration and drift-probe tests run in the default suite.

`test_builder.py::test_alphagenome_round_trip` builds the fixture with default parameters, compares the schema and sample rows to the snapshots, and checks an independent oracle: `splice_sites.max_abs_raw` equals `ag_splice_sites` (max |raw_score| over SPLICE_SITES rows, from an independent per-variant ClinVar benchmark computed by the maintainer) for each fixture variant. On the full shard the builder reproduces it for all 339 variants that have it. A second build with `ontology_curies=["UBERON:0006566", "UBERON:0006631"]` pins `splice_site_usage.max_abs_raw` to its curie-defined values, cross-checked pandas against Hail on the full shard. They differ by design from that benchmark's GTEx-only `ag_heart_ssu` for 2 of the 3 variants, because the ENCODE heart track scores higher there (§ 4).
