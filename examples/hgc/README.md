# HGC (Hail Genotype Combiner) Examples

This directory contains examples and benchmarks for the HGC joint genotyping pipeline.

## Directory Structure

```
examples/hgc/
├── README.md                  # This file
├── qc/                        # Quality control examples
│   └── hgc_qc_example.py     # QC workflow for joint-called cohorts
├── scalability/               # Sample size scalability benchmark
│   ├── README.md             # Scalability benchmark documentation
│   ├── benchmark.py          # Python workflow runner
│   ├── benchmark.sh          # Shell orchestration script
│   ├── plot_results.py       # Results visualization
│   └── run_example.sh        # Quick start script
├── cpu_scaling/               # CPU scaling benchmark
│   ├── benchmark.py          # Python workflow runner
│   ├── benchmark.sh          # Shell orchestration script
│   ├── plot_results.py       # Results visualization
│   └── run_example.sh        # Quick start script
├── common/                    # Shared utilities
│   ├── test_benchmark_setup.py  # Environment validation
│   ├── setup_hvantk_env.sh      # Environment setup
│   ├── extract_timing.py        # Timing utilities
│   └── generate_sample_list.sh  # Sample list generation
└── results/                   # Example output files
```

## Quick Start

### 1. Validate Environment

```bash
python examples/hgc/common/test_benchmark_setup.py
```

### 2. Run QC Example

```bash
# Activate environment
poetry shell

# Run QC workflow (uses test data)
python examples/hgc/qc/hgc_qc_example.py

# Check outputs
ls examples/hgc/results/
```

### 3. Run Benchmarks

**Scalability benchmark** (tests performance vs. sample count):

```bash
cd examples/hgc/scalability
bash run_example.sh
# Or run directly:
bash benchmark.sh --gvcf-dir /path/to/gvcfs --output-dir ./results
```

**CPU scaling benchmark** (tests performance vs. CPU count):

```bash
cd examples/hgc/cpu_scaling
bash run_example.sh
# Or run directly:
bash benchmark.sh --gvcf-list samples.txt --output-dir ./results
```

## Components

### QC Workflow (`qc/`)

Demonstrates quality control for joint-called cohorts:
- Computing QC metrics
- Generating visualizations
- Creating HTML reports
- Quality-based filtering strategies

### Scalability Benchmark (`scalability/`)

Tests how HGC performance scales with cohort size:
- Runs HGC workflow with varying sample counts
- Measures timing for each step
- Generates scaling plots

See [scalability/README.md](scalability/README.md) for detailed documentation.

### CPU Scaling Benchmark (`cpu_scaling/`)

Tests strong scaling (speedup vs. CPU count):
- Runs HGC workflow with fixed cohort size
- Varies CPU core count
- Measures speedup and efficiency

### Common Utilities (`common/`)

Shared scripts for all benchmarks:
- **test_benchmark_setup.py** - Validates environment setup
- **setup_hvantk_env.sh** - Creates conda environment
- **extract_timing.py** - Extracts timing data from logs
- **generate_sample_list.sh** - Generates sample lists from GVCF directory

## Expected Outputs

### QC Workflow

- `qc_report_*.html` - Interactive HTML report
- `qc_dashboard_*.png` - Multi-panel QC visualization

### Benchmarks

- Timing metrics (JSON/CSV format)
- Performance plots (PNG format)
- Scalability analysis summaries

## Documentation

- [HGC Documentation](../../docs/tools/hgc.md)
- [Architecture Overview](../../docs/ARCHITECTURE.md)
- [Scalability Guide](scalability/README.md)

## Requirements

- Hail 0.2.x
- Spark configured with appropriate resources
- For benchmarks: sufficient CPU cores and memory
- For QC: MatrixTable with sample and variant data
