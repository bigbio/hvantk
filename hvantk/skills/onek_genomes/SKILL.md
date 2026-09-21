---
name: hvantk:resource-onek-genomes
description: 1000 Genomes high-coverage (NYGC/CCDG) callset — per-chromosome VCFs as a VariantMatrix, plus IGSR sample metadata as an AnnotationTable.
status: provisional
backend: hail
domain: genomics
---

# 1000 Genomes (high-coverage NYGC/CCDG callset)

Read `hvantk/skills/_conventions/SKILL.md` first. This skill assumes its repository map, helpers, keying conventions, builder pattern, and validation contract.

## 1. Status & scope

Provisional. This plugin ships **two datasets**, both `domain: genomics` / `backend: hail` per `plugin.yaml`:

- `onek-genomes:variants` — a `VariantMatrix` built from per-chromosome bgzipped VCFs of the 1000 Genomes Project high-coverage callset (release `1000G_2504_high_coverage`, 2504 samples, GRCh38). **BYO data**: no downloader ships for this dataset (see § 2).
- `onek-genomes:samples` — an `AnnotationTable` of IGSR canonical sample metadata (population, super-population, sex). Auto-downloaded via a lifecycle stage.

Out of scope for this skill:
- Joining `:samples` onto `:variants`. This is a deliberate post-load, user-side step (not baked into the `:variants` build) — see the worked example in § 7.
- Downloading the ~1.5 TB genotype VCF set programmatically — out of the "build a downloader" tier per `CLAUDE.md`'s size threshold; the plugin has no `cli.py` and no `lifecycle.download` for `:variants`.

## 2. Source identity

1000 Genomes Project high-coverage NYGC/CCDG callset, release directory `1000G_2504_high_coverage` (2504 samples, GRCh38).

- Variant VCFs: no downloader is shipped, so no single URL is authoritative for acquisition — the previous revision of this skill pointed users at the NCBI FTP mirror (`https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/`). The plugin's own code (the drift probe) instead anchors against the EBI mirror: `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz` (`hvantk/skills/onek_genomes/drift_probe.py`). Both are mirrors of the same NYGC/CCDG release; either is usable to populate `--raw-dir`. Per-chromosome files follow the pattern `..._chr<N>.filtered.shapeit2-duohmm-phased.vcf.gz` with a matching `.tbi`.
- Sample panel: `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel` (`hvantk/skills/onek_genomes/samples_download.py`, `_IGSR_URL`; matches `drift_probe.py`'s `_SAMPLES_URL`).
- No `catalog/datasets.json` exists for this plugin (verified: no `catalog:` key in `plugin.yaml`, no `catalog/` folder, and no reference to `onek-genomes`/`onek_genomes` anywhere under `hvantk/resources/`) — this is a catalog gap in the same sense as `hgnc`'s (see `hvantk/skills/hgnc/SKILL.md` § 2).

## 3. Backend choice + reasoning

`hail` for both datasets. `:variants` is keyed `(locus, alleles)` per `_conventions` § 3's variant-domain convention. `:samples` is a small lookup table that could in principle be pandas, but is built as a Hail Table keyed by `sample` because it is designed to be joined directly onto `:variants`' Hail `MatrixTable` column key (`s`) — see the worked join in § 7, carried over from the plugin's own documentation.

## 4. Raw format & gotchas

**`:variants`** (`hvantk/skills/onek_genomes/builder.py`):
- Input is a directory of bgzipped VCFs matching `*.vcf.gz` (`_discover_vcf_files`); each must have a sibling `.tbi` index, or the builder raises `FileNotFoundError` listing every file missing one.
- Chromosome token extraction uses `_CHROM_PATTERN = re.compile(r"(?i)chr[_]?(\d+|X|Y|M|MT)\b")` against the basename — matches `chr22`, `chr_1`, `CHR1`, etc., case-insensitively. Files are then re-sorted into biological order (autosomes 1–22, then X, Y, M/MT) via `_chrom_sort_key`, not left in glob/alphabetical order (which would otherwise put `chr10` before `chr2`).
- `chromosomes` plugin-arg accepts either a `list[str]` or a comma-separated string; a string is split on `,` — an explicit workaround for issue #119 in the `--plugin-arg` coercion path (`build_onek_genomes_variants`). Filtering matches the extracted chrom token case-insensitively (uppercased).
- If `chromosomes` is not passed, the builder only *warns* (does not fail) when any of the 24 standard chromosomes (`chr1`–`chr22`, `chrX`, `chrY`) are missing from the discovered set — a partial `--raw-dir` silently builds a partial cohort unless the caller checks the log.
- `auto_convert_bgz` (default `False`) is forwarded per-file, in parallel via a `ThreadPoolExecutor`, to `resolve_compression()` (`hvantk.core.utils.file_utils`), which re-bgzips files that are gzip-compressed but not valid BGZF blocks; the resulting `force_bgz` flag (false if *any* file needed non-bgz handling) is passed through to `hl.import_vcf`.
- `hl.import_vcf(..., array_elements_required=False, reference_genome=params.get("reference_genome", "GRCh38"))` — `array_elements_required=False` is needed because the INFO struct's per-population `AC`/`AF`/`AN`/`AC_Het`/`AC_Hom` arrays can be sparse for subpopulations with zero coverage at a site.
- The committed fixture (`tests/testdata/raw/onek_genomes/chr22.vcf.gz`) confirms `#CHROM` values carry the `chr` prefix (e.g. `chr22`), consistent with the `GRCh38` default. The VCF `ID` column, by contrast, uses **un**-prefixed contig names (e.g. `22:10519265:CA:C`) — the builder does not touch `ID`, but a downstream consumer assuming `ID` mirrors the locus's `chr`-prefixed contig will be wrong.
- No field renaming: `filters` (`set<str>`), the `info` struct, and the `GT: call` entry field are exactly what `hl.import_vcf` produces natively.

**`:samples`** (`hvantk/skills/onek_genomes/samples_builder.py`, `samples_download.py`):
- The builder looks for a file named exactly `igsr_samples.tsv` under `parsed_input`/`raw_dir` — hardcoded, not glob-discovered. Missing it raises `FileNotFoundError` with a hint to re-run without `--skip-download` or place the file manually.
- The upstream IGSR panel ships trailing tabs on its header line; `download_igsr_samples` strips trailing whitespace per line before writing (`cleaned = "\n".join(line.rstrip() for line in resp.text.splitlines()) + "\n"`) "so Hail's import_table sees consistent field counts across header and data" (comment in `samples_download.py`). A user who fetches the raw panel by hand instead of via the lifecycle stage must replicate this strip, or `hl.import_table` may see a header/data field-count mismatch.
- Imported with `hl.import_table(..., impute=True, key="sample")`. Despite `impute=True`, none of the four columns (`sample`, `pop`, `super_pop`, `gender`) infer as numeric — all remain `str` (`tests/snapshots/samples_schema.json`).
- The committed fixture ships the full 2504-sample panel (2505 lines including header), matching the release's stated cohort size.

## 5. Output contract

**`:variants`** — `VariantMatrix` wrapping a Hail `MatrixTable`, `schema_id="onek-genomes-variants-v1"`. Row key `locus<GRCh38>` + `alleles: array<str>`; column key `s: str` (sample ID); entry field `GT: call`. Row also carries `qual: float64`, `rsid: str`, `filters: set<str>`, and an `info` struct of 89 fields (`tests/snapshots/variants_schema.json`): population-stratified `AC`/`AF`/`AN`/`AC_Het`/`AC_Hom` (suffixed `_AFR`/`_AMR`/`_EAS`/`_EUR`/`_SAS`, each with an `_unrel` unrelated-subset variant), and site-quality/VQSR fields (`DP`, `QD`, `FS`, `MQ`, `MQ0`, `MQRankSum`, `ReadPosRankSum`, `BaseQRankSum`, `ClippingRankSum`, `SOR`, `HaplotypeScore`, `InbreedingCoeff`, `ExcHet`(+ per-population), `HWE`(+ per-population), `VQSLOD`, `culprit`, `NEGATIVE_TRAIN_SITE`, `POSITIVE_TRAIN_SITE`, `MLEAC`, `MLEAF`, `RAW_MQ`, `VariantType`, `END`, `DS`).

**`:samples`** — `AnnotationTable` wrapping a Hail Table keyed by `sample: str`, `schema_id="onek-genomes-samples-v1"`. Row fields: `pop: str`, `super_pop: str`, `gender: str` (`tests/snapshots/samples_schema.json`).

## 6. hvantk integration points

- Plugin manifest: `hvantk/skills/onek_genomes/plugin.yaml`, dataset keys `onek-genomes:variants` and `onek-genomes:samples`.
- Builders: `build_onek_genomes_variants` in `hvantk/skills/onek_genomes/builder.py`; `build_onek_genomes_samples` in `hvantk/skills/onek_genomes/samples_builder.py`. Both are `(parsed_input, ctx, **params) -> <Artifact>`.
- Downloader: only `:samples` has one — `download_igsr_samples` in `hvantk/skills/onek_genomes/samples_download.py`, wired as `lifecycle.download` for that dataset in `plugin.yaml`. There is no `cli.py` and no `cli:` block for this plugin, so neither dataset appears under `hvantk download`. Supply `:variants` yourself (`acquisition.mode: byo`).
- Drift probes: `fetch_variants_fingerprint` and `fetch_samples_fingerprint`, both in `hvantk/skills/onek_genomes/drift_probe.py` — two independent probes, matching the two separate `drift_probe:` entries in `plugin.yaml`.
- Tests: `hvantk/skills/onek_genomes/tests/test_builder.py` builds both datasets from the shared fixture directory and asserts against per-dataset snapshots via `phase_b_snapshot_adapter` (`hvantk.tests._snapshot_utils`), which bridges the `(parsed_input, ctx, **params)` builder signature to the legacy `(input_path, output_path, **kw)` snapshot-test calling convention.

## 7. Workflow steps

**`:variants`** (BYO data):
1. Populate a directory with per-chromosome `*.vcf.gz` + `.tbi` files from a 1000 Genomes high-coverage mirror (§ 2).
2. Build: `hvantk reprocess onek-genomes:variants --raw-dir /data/1kg/vcfs/ --output /data/1kg.mt --skip-download --plugin-arg reference_genome=GRCh38 [--plugin-arg chromosomes=chr1,chr2,chrX]`.
3. Watch the build log for the "standard chromosomes not found" warning if a full-cohort build was expected (§ 4) — it will not fail the build.

**`:samples`**:
1. Build (download + build in one step): `hvantk reprocess onek-genomes:samples --raw-dir /tmp/igsr --output /data/1kg_samples.ht`.

**Joining samples onto the variant cohort** (post-load, user-side):
```python
import hail as hl
from hvantk.core.models import VariantMatrix, AnnotationTable

vm = VariantMatrix.load("/data/1kg.mt")
ann = AnnotationTable.load("/data/1kg_samples.ht").to_hail()
mt = vm.to_hail_mt()
mt = mt.annotate_cols(sample_annotations=ann[mt.s])
# mt.sample_annotations.{super_pop,pop,gender}
```
This join is intentionally not baked into the `:variants` build: sample metadata is user-provided in many real workflows, and keeping it separate means the `:variants` build-time fingerprint attests only to genotype-release identity.

## 8. Update playbook

Both datasets use **release-identity** drift probes (`drift_probe.py` module docstring): the upstream sources are immutable per release, so nothing should drift *within* a release; a new release requires a code change, not a fingerprint regeneration alone.

1. `:variants` — when a successor release ships (e.g., a larger-cohort expansion), update `_VARIANTS_RELEASE_DIR` / `_VARIANTS_ANCHOR_FILE` / `_VARIANTS_BASE` in `drift_probe.py` to repoint at the new release's chr22 anchor, then regenerate `tests/drift_fingerprint_variants.json` via `hvantk drift --regenerate onek-genomes:variants`.
2. `:samples` — if IGSR republishes the panel at a new URL/version, update `_SAMPLES_URL` / `_SAMPLES_VERSION` in `drift_probe.py` (and `_IGSR_URL` in `samples_download.py`, which must stay in sync manually — the two are not derived from one shared constant), then regenerate `tests/drift_fingerprint_samples.json`.
3. If VCF INFO/format fields or the samples panel's columns change shape, regenerate both snapshot pairs with `pytest hvantk/skills/onek_genomes/tests -m hail --regenerate-snapshots` and review the diff — note this plugin's snapshot filenames (`variants_schema.json`, `samples_sample_rows.json`, …) and flat head-of-N row shape predate `_snapshot_utils.regenerate_snapshots`'s key-matched `{"key": ..., "row": ...}` convention used elsewhere; the regeneration path in `test_builder.py` deliberately preserves the older per-dataset-filename convention rather than migrating it.

## 9. Validation contract

Declared in `plugin.yaml`'s `tests:` block, one per dataset, both under a **plugin-local** fixture directory shared between the two datasets (not the shared repo-level form — contrast with `gnomad_metrics`, § 9 of that skill):

- `onek-genomes:variants`:
  - `fixture`: `tests/testdata/raw/onek_genomes` (contains `chr22.vcf.gz` + `.tbi`)
  - `schema_snapshot`: `tests/snapshots/variants_schema.json`
  - `row_snapshot`: `tests/snapshots/variants_sample_rows.json`
  - `drift_fingerprint`: `tests/drift_fingerprint_variants.json`
  - `command`: `pytest hvantk/skills/onek_genomes/tests -m hail`
- `onek-genomes:samples`:
  - `fixture`: `tests/testdata/raw/onek_genomes` (same directory; contains `igsr_samples.tsv`)
  - `schema_snapshot`: `tests/snapshots/samples_schema.json`
  - `row_snapshot`: `tests/snapshots/samples_sample_rows.json`
  - `drift_fingerprint`: `tests/drift_fingerprint_samples.json`
  - `command`: `pytest hvantk/skills/onek_genomes/tests -m hail`

Per `_conventions` § 12, these two datasets deliberately do **not** share one drift fingerprint: each has its own probe function (`fetch_variants_fingerprint` vs `fetch_samples_fingerprint`) and its own baseline file, because they are independent upstream signals (a VCF release vs. a samples panel) rather than one provider-wide signal fanned out to variants.
