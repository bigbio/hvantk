# 1000 Genomes Reference Build

Scripts for building a Hail MatrixTable from 1000 Genomes NYGC high-coverage data, for use as a reference panel (e.g., ancestry inference).

## Contents

| Script | Description |
|--------|-------------|
| `build_1kg_nygc.py` | Stages NYGC-pattern VCFs and shells out to `hvantk reprocess onek-genomes:variants` |
| `build_1kg_nygc_cli.py` | Same as `build_1kg_nygc.py`, kept as a separate file for backward compatibility |

## Prerequisites

Download the NYGC/CCDG 1000 Genomes (2020) per-chromosome recalibrated VCFs into a local directory.

## Quick Start

```bash
# Python API
python examples/1k_genome/build_1kg_nygc.py \
  --vcf-dir /data/1kg/vcfs \
  --output-mt /data/1kg/out.mt

# With chromosome selection and sample annotations
python examples/1k_genome/build_1kg_nygc.py \
  --vcf-dir /data/1kg/vcfs \
  --output-mt /data/1kg/out.mt \
  --chromosomes chr1,chr2,chrX \
  --sample-annotations samples.ped \
  --overwrite

# CLI wrapper (tests the hvantk CLI entry point)
python examples/1k_genome/build_1kg_nygc_cli.py \
  --vcf-dir /data/1kg/vcfs \
  --output-mt /data/1kg/out.mt
```

## Requirements

- Hail 0.2.x with Spark configured
- Sufficient memory and disk for full 1000 Genomes data
