# 1000 Genomes (high-coverage NYGC/CCDG callset)

This plugin ships two datasets from the IGSR-hosted 1000 Genomes Project high-coverage
callset (release `1000G_2504_high_coverage`, 2504 samples, GRCh38):

- `onek-genomes:variants` — a `VariantMatrix` of per-chromosome bgzipped VCFs
  imported into a single Hail MatrixTable. **BYO data** — see below.
- `onek-genomes:samples` — an `AnnotationTable` of IGSR canonical sample metadata
  (population, super-population, sex, family relationships). Downloaded automatically.

## Downloading

### Variants (`:variants`)

The genotype VCFs total ~1.5 TB across 24 per-chromosome files. The plugin does not
fetch them; users supply a pre-populated directory. Canonical source:

    https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/

Each file follows the naming pattern
`20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr<N>.filtered.shapeit2-duohmm-phased.vcf.gz`
with a matching `.tbi` index.

### Samples (`:samples`)

Auto-downloaded from IGSR's canonical samples endpoint (small TSV, ~250 KB).

## Usage

Build the variant cohort (BYO data):

    hvantk reprocess onek-genomes:variants \
        --raw-dir /data/1kg/vcfs/ \
        --output /data/1kg.mt \
        --skip-download \
        --plugin-arg reference_genome=GRCh38 \
        --plugin-arg chromosomes=chr1,chr2,chrX     # optional subset

Build the sample-metadata table:

    hvantk reprocess onek-genomes:samples \
        --raw-dir /tmp/igsr \
        --output /data/1kg_samples.ht

## Post-load: joining samples onto the variant cohort

```python
import hail as hl
from hvantk.core.models import VariantMatrix, AnnotationTable

vm = VariantMatrix.load("/data/1kg.mt")
ann = AnnotationTable.load("/data/1kg_samples.ht").to_hail()
mt = vm.to_hail_mt()
mt = mt.annotate_cols(sample_annotations=ann[mt.s])
# mt now has mt.sample_annotations.super_population, .population, etc.
```

The join is a deliberate post-load step (not baked into the `:variants` artifact)
because (a) sample metadata is user-provided in many real workflows, and
(b) the build-time fingerprint then cleanly attests to just the genotype release
identity rather than mixing in an annotations file.

## Plugin parameters

`:variants`:

| `--plugin-arg` | Type | Default | Notes |
|---|---|---|---|
| `reference_genome` | str | `"GRCh38"` | Passed to `hl.import_vcf` |
| `chromosomes` | str (comma-sep) | (all) | Workaround for #119 — split inside the builder |
| `auto_convert_bgz` | bool | `false` | Re-bgz any non-BGZF-clean VCFs before import |

`:samples`: no parameters.

## Drift behavior

The drift probe for `:variants` HEADs a canonical anchor file (chr22 VCF) at the
1KG release URL and captures the release directory name (`1000G_2504_high_coverage`)
as `source_version`. This is a **release-identity probe** — 1KG releases are
immutable per release, so within a release nothing changes. When a successor
release lands (e.g., a 3202-sample expansion), bump the constants in
`drift_probe.py` to repoint at the new release.

The probe for `:samples` HEADs the IGSR samples TSV URL and captures
Last-Modified + Content-Length headers.
