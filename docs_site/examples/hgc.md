# HGC Examples

The HGC (Hail Genotype Combiner) pipeline provides joint genotyping for GVCF cohorts with quality control and benchmarking support.

See the [HGC reference](../tools/hgc.md) for full CLI options and architecture details.

## QC Workflow

Compute quality metrics and generate reports for a joint-called cohort:

```bash
# Using the CLI pipeline (recommended)
hvantk hgc pipeline -i /path/to/gvcfs -o /path/to/output

# Or step by step
hvantk hgc gvcf-combine -g /data/gvcfs -o cohort.vds
hvantk hgc vds2mt -i cohort.vds -o cohort.mt --adjust-genotypes
hvantk hgc compute-qc -i cohort.mt -o cohort_qc.mt
hvantk hgc qc-report -i cohort_qc.mt -o report.html
```

### Python API

```python
from hvantk.algorithms.hgc.pipeline import PipelineConfig, PipelineRunner

config = PipelineConfig(
    input_dir="/data/gvcfs",
    output_dir="/data/output",
)
runner = PipelineRunner(config)
runner.run()
```

### Expected Outputs

- `qc_report_*.html` - Static HTML QC report (embedded plots, tables, recommendations)

For a publication figure or any custom plot, export the QC tables
(`QCMetrics.get_sample_metrics_df()` / `get_variant_metrics_df()`) and plot with matplotlib —
see [Plot your own QC from the tables](../tools/hgc.md#plot-your-own-qc-from-the-tables).

## Scalability Benchmarks

Two benchmark suites are included:

**Sample scalability** - measures performance as cohort size grows:

```bash
cd examples/hgc/scalability
bash benchmark.sh --gvcf-dir /path/to/gvcfs --output-dir ./results
```

**CPU scaling** - measures strong scaling (speedup vs. CPU count):

```bash
cd examples/hgc/cpu_scaling
bash benchmark.sh --gvcf-list samples.txt --output-dir ./results
```

Both produce timing metrics (JSON/CSV) and scaling plots (PNG).

## Requirements

- Hail 0.2.x with Spark configured
- Sufficient CPU cores and memory for benchmarks
- MatrixTable with sample and variant data for QC

## Runnable Scripts

See the [`examples/hgc/`](https://github.com/bigbio/hvantk/tree/main/examples/hgc/) directory for all scripts and benchmark utilities.
