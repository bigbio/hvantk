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

- `qc_report_*.html` - Self-contained HTML report with embedded plots

### Custom Plots from Exported QC Tables

There is no standalone plotting/dashboard command. For ad hoc or custom plots,
export the QC metrics as pandas DataFrames and plot with matplotlib directly:

```python
import hail as hl
import matplotlib.pyplot as plt
from hvantk.algorithms.hgc import compute_full_qc

mt = hl.read_matrix_table("cohort_qc.mt")
qc = compute_full_qc(mt)
df = qc.get_sample_metrics_df()

fig, ax = plt.subplots()
ax.hist(df["sample_qc.call_rate"], bins=30)
ax.axvline(0.85, color="red", ls="--", label="min call rate")
ax.set_xlabel("Sample call rate")
ax.legend()
fig.savefig("call_rate.png", dpi=300, bbox_inches="tight")
```

See [Quality Control Functions](../tools/hgc.md#quality-control-functions-post-combination)
in the HGC reference for the full set of exported metrics.

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
