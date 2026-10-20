# 1000 Genomes Reference Build

Scripts for building a Hail MatrixTable from 1000 Genomes NYGC high-coverage data, for use as a reference panel (e.g., ancestry inference).

## Contents

| Script | Description |
|--------|-------------|
| `build_1kg_nygc.py` | Stages NYGC-pattern VCFs and shells out to `hvantk reprocess onek-genomes:variants` |
| `build_1kg_nygc_cli.py` | Same as `build_1kg_nygc.py`, kept as a separate file for backward compatibility |

## Prerequisites

Download the NYGC/CCDG 1000 Genomes high-coverage per-chromosome phased VCFs
(`*_chr<N>.filtered.shapeit2-duohmm-phased.vcf.gz` + `.tbi`) into a local directory. See
[1000 Genomes](https://bigbio.github.io/hvantk/guide/data-sources/#1000-genomes-nygcccdg-high-coverage-callset)
for the mirror URL — `onek-genomes:variants` is BYO (no downloader ships for it; files are
multi-GB per chromosome).

## Quick start

Both scripts stage the VCFs into a temp directory, then shell out to
`hvantk reprocess onek-genomes:variants --skip-download` as a subprocess — neither calls
a Python build function directly.

```bash
python examples/1k_genome/build_1kg_nygc.py \
  --vcf-dir /data/1kg/vcfs \
  --output-mt /data/1kg/out.mt

# With chromosome selection
python examples/1k_genome/build_1kg_nygc.py \
  --vcf-dir /data/1kg/vcfs \
  --output-mt /data/1kg/out.mt \
  --chromosomes chr1,chr2,chrX \
  --overwrite

# CLI wrapper (same staging logic, exercises the hvantk CLI entry point with
# friendlier error handling if `hvantk` is missing from PATH or exits non-zero)
python examples/1k_genome/build_1kg_nygc_cli.py \
  --vcf-dir /data/1kg/vcfs \
  --output-mt /data/1kg/out.mt
```

`--sample-annotations` is still accepted for compatibility but is **not** passed to the
build (the plugin build step takes no such parameter) — the script logs a warning and
ignores it. To attach sample metadata (population, super-population, sex), build
`onek-genomes:samples` separately and join it post-load:

```bash
hvantk reprocess onek-genomes:samples --raw-dir /tmp/igsr --output /data/1kg_samples.ht
```

```python
from hvantk.core.models import VariantMatrix, AnnotationTable

vm = VariantMatrix.load("/data/1kg/out.mt")
ann = AnnotationTable.load("/data/1kg_samples.ht").to_hail()
mt = vm.to_hail_mt()
mt = mt.annotate_cols(sample_annotations=ann[mt.s])
# mt.sample_annotations.{super_pop,pop,gender}
```

That `super_pop`/`pop` column is the reference panel's ancestry label — pass its name as
`--ancestry-col` / `ancestry_col` to [ancestry inference](https://bigbio.github.io/hvantk/tools/ancestry/).

## Requirements

- Hail 0.2.x with Spark configured
- Sufficient memory and disk for full 1000 Genomes data
