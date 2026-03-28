# HGC Examples

Examples and benchmarks for the HGC joint genotyping pipeline.

For full documentation, see the [HGC docs](https://bigbio.github.io/hvantk/tools/hgc/) and [HGC examples guide](https://bigbio.github.io/hvantk/examples/hgc/).

## Directory Structure

```
examples/hgc/
├── qc/                        # Quality control examples
│   └── hgc_qc_example.py     # QC workflow for joint-called cohorts
├── scalability/               # Sample size scalability benchmark
│   ├── README.md
│   ├── benchmark.py
│   ├── benchmark.sh
│   ├── plot_results.py
│   └── run_example.sh
├── cpu_scaling/               # CPU scaling benchmark
│   ├── benchmark.py
│   ├── benchmark.sh
│   ├── plot_results.py
│   └── run_example.sh
├── common/                    # Shared utilities
│   ├── test_benchmark_setup.py
│   ├── setup_hvantk_env.sh
│   ├── extract_timing.py
│   └── generate_sample_list.sh
└── results/                   # Example output files
```

## Quick Start

```bash
# Validate environment
python examples/hgc/common/test_benchmark_setup.py

# Run QC workflow (uses test data)
python examples/hgc/qc/hgc_qc_example.py

# Run scalability benchmark
cd examples/hgc/scalability && bash run_example.sh

# Run CPU scaling benchmark
cd examples/hgc/cpu_scaling && bash run_example.sh
```

## Expected Outputs

- `qc_report_*.html` - Interactive HTML report
- `qc_dashboard_*.png` - Multi-panel QC visualization
- Timing metrics and scaling plots (benchmarks)

## Requirements

- Hail 0.2.x with Spark configured
- For benchmarks: sufficient CPU cores and memory
- For QC: MatrixTable with sample and variant data
