# HGC (Hail Genotype Combiner) Examples

This directory contains examples and benchmarks for the HGC joint genotyping pipeline.

## Contents

### Examples

**`hgc_qc_example.py`** - Quality control workflow for joint-called cohorts

Demonstrates:
- Computing QC metrics for combined cohorts
- Generating QC visualizations
- Creating HTML reports with interactive plots
- Quality-based filtering strategies

```bash
python examples/hgc/hgc_qc_example.py
```

### Benchmarks

**`hgc_cpu_scaling_benchmark.py`** - CPU scaling benchmark

Tests how HGC performance scales with different CPU core counts.

**`hgc_scalability_benchmark.py`** - Sample size scalability benchmark

Tests how HGC performance scales with cohort size (number of samples).

**`test_benchmark_setup.py`** - Benchmark environment validation

Verifies that the environment is properly configured for running benchmarks.

### Supporting Scripts

The `scripts/` directory contains:
- `plot_cpu_scaling_results.py` - Visualization for CPU scaling results
- `plot_scalability_results.py` - Visualization for scalability results
- `extract_timing.py` - Timing data extraction utilities
- Shell scripts for running benchmarks (`run_*.sh`)
- Environment setup scripts (`setup_*.sh`)

## Quick Start

### Run QC Example

```bash
# Activate environment
poetry shell

# Run QC workflow (uses test data)
python examples/hgc/hgc_qc_example.py

# Check outputs
ls examples/hgc/results/
```

### Run Benchmarks

```bash
# Validate environment first
python examples/hgc/test_benchmark_setup.py

# Run CPU scaling benchmark
python examples/hgc/hgc_cpu_scaling_benchmark.py

# Run sample scalability benchmark
python examples/hgc/hgc_scalability_benchmark.py
```

## Expected Outputs

### QC Workflow

- `qc_report_*.html` - Interactive HTML report with:
  - Summary statistics
  - Sample and variant QC metrics
  - Ti/Tv ratio analysis
  - Quality recommendations
- `qc_dashboard_*.png` - Multi-panel QC visualization

### Benchmarks

- Timing metrics (JSON format)
- Performance plots (PNG format)
- Scalability analysis results

## Documentation

For detailed documentation on the HGC pipeline:
- [HGC Documentation](../../docs/tools/hgc.md)
- [Architecture Overview](../../docs/ARCHITECTURE.md)
- [Scalability Guide](./README_scalability.md)

## Customization

### QC Example

Modify `hgc_qc_example.py` to use your own data:

```python
# Change input paths
input_mt = "path/to/your/cohort.mt"
output_dir = "path/to/output/"
```

### Benchmarks

Edit benchmark parameters:

```python
# In hgc_cpu_scaling_benchmark.py
cpu_counts = [2, 4, 8, 16]  # Test different core counts

# In hgc_scalability_benchmark.py
sample_sizes = [100, 500, 1000, 5000]  # Test different cohort sizes
```

## Requirements

- Hail 0.2.x
- Spark configured with appropriate resources
- For benchmarks: sufficient CPU cores and memory
- For QC: MatrixTable with sample and variant data

## Tips

1. **QC workflow**: Start with small test datasets to verify the pipeline
2. **Benchmarks**: Ensure adequate system resources (CPU, memory, disk)
3. **Results**: Check `results/` directory for all outputs
4. **Customization**: Copy scripts and modify for your specific needs
