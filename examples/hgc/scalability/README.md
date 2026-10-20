# HGC Scalability & Performance Benchmarks

This directory contains scripts to benchmark the HGC (Hail-based Genotype Combiner) workflow performance.

## Available Benchmarks

### 1. Sample Size Scalability (Weak Scaling)
**Purpose**: Measure how runtime scales with increasing cohort size  
**Varies**: Sample count (20, 50, 100, 250, 500, 750, 1000)  
**Fixed**: CPU count (144 CPUs)  
**Script**: `benchmark.sh`

### 2. CPU Scaling (Strong Scaling)
**Purpose**: Measure speedup and efficiency with varying CPU count  
**Varies**: CPU count (16, 24, 32, 48, 64, 96, 144)  
**Fixed**: Sample count (e.g., 500 samples)  
**Script**: `../cpu_scaling/benchmark.sh`

---

## Overview

The benchmarks measure the complete end-to-end HGC workflow performance:

1. **GVCF → VDS**: Combine individual GVCF files into a Variant DataSet
2. **VDS → MatrixTable**: Convert VDS to analysis-ready MatrixTable
3. **Compute QC**: Calculate sample and variant quality control metrics
4. **Export VCF**: Generate final cohort VCF file

### Metrics Collected

- **Runtime**: Wall-clock time for each step and total workflow
- **Memory**: Peak memory usage per run
- **Scalability**: How performance scales with sample size or CPU count
- **Speedup**: Performance gain from additional CPUs (CPU scaling only)
- **Efficiency**: CPU utilization effectiveness (CPU scaling only)

## Files

### Sample Size Scalability
```text
├── benchmark.sh       # Orchestration script (sample size scaling)
├── benchmark.py       # Python workflow runner
├── plot_results.py        # Visualization script
└── README.md                          # This file
```

### CPU Scaling
```text
├── ../cpu_scaling/benchmark.sh       # Orchestration script (CPU scaling)
├── ../cpu_scaling/benchmark.py       # Python workflow runner with CPU config
├── ../cpu_scaling/plot_results.py        # CPU scaling visualization
├── ../common/generate_sample_list.sh            # Helper to create fixed sample lists
└── ../cpu_scaling/run_example.sh         # Complete example workflow
```

## Requirements

- Python 3.8+
- hvantk library installed
- Hail installed and configured
- Required Python packages: matplotlib, seaborn, pandas, numpy, scipy

## Usage

### Sample Size Scalability Benchmark

#### 1. Run the Benchmark

```bash
cd examples/hgc/scalability

# Run with default settings (the GVCF directory is required)
bash benchmark.sh --gvcf-dir /path/to/gvcfs

# Or customize parameters
bash benchmark.sh \
    --gvcf-dir /path/to/gvcfs \
    --output-dir ./my_results \
    --sample-sizes 20,50,100,250,500,750,1000 \
    --reference GRCh38 \
    --seed 42
```

**Default Configuration:**
- GVCF directory: none; pass `--gvcf-dir` or export `GVCF_DIR`
- Output directory: `./scalability_results`
- Sample sizes: 20, 50, 100, 250, 500, 750, 1000
- Reference genome: GRCh38
- Random seed: 42 (for reproducible sampling)

**Options:**
- `--gvcf-dir DIR`: Directory containing GVCF files
- `--output-dir DIR`: Output directory for results
- `--sample-sizes SIZES`: Comma-separated sample sizes to test
- `--reference REF`: Reference genome (GRCh38 or GRCh37)
- `--seed SEED`: Random seed for GVCF sampling
- `--resume`: Resume from previous run (skip completed runs)
- `--help`: Show help message

#### 2. Generate Plots

After the benchmark completes:

```bash
python plot_results.py --results-dir ./scalability_results
```

This generates:
- Runtime breakdown (stacked bar chart)
- Total runtime vs sample size (with trend line)
- Step-by-step comparison (line plot)
- Scaling efficiency (time per sample)
- Memory usage vs sample size
- Summary report (text file)

---

### CPU Scaling Benchmark

#### 1. Prepare Fixed Sample List

```bash
# Generate a fixed list of 500 GVCFs
bash ../common/generate_sample_list.sh /path/to/gvcfs 500 samples_500.txt
```

#### 2. Run CPU Scaling Benchmark

```bash
# Run with default CPU counts (16,24,32,48,64,96,144)
bash ../cpu_scaling/benchmark.sh \
    --gvcf-list samples_500.txt \
    --output-dir ./cpu_scaling_500 \
    --sample-size 500 \
    --cpu-counts 16,24,32,48,64,96,144

# Or use the example script
bash ../cpu_scaling/run_example.sh
```

**Default Configuration:**
- Sample size: 500 (FIXED)
- CPU counts: 16, 24, 32, 48, 64, 96, 144 (VARIABLE)
- Driver memory: 128GB (FIXED)
- Shuffle partitions: 2048 (FIXED)
- Max partition bytes: 128MB
- Reference genome: GRCh38

**Options:**
- `--gvcf-list FILE`: File with GVCF paths (one per line)
- `--output-dir DIR`: Output directory for results
- `--sample-size N`: Fixed sample size (for naming/logging)
- `--cpu-counts COUNTS`: Comma-separated CPU counts to test
- `--driver-memory GB`: Driver memory in GB, kept fixed across runs (default: 128)
- `--shuffle-partitions N`: Shuffle partition count, kept fixed across runs (default: 2048)
- `--max-partition-bytes MB`: Max partition bytes for input reading (default: 128)
- `--local-dir DIR`: Spark local directory for shuffle spill (recommend fast NVMe)
- `--aqe-enabled`: Enable Adaptive Query Execution (default: off, for pure scaling)
- `--reference REF`: Reference genome
- `--conda-env ENV`: Conda environment name (default: auto-detect, or `hvantk`)
- `--resume`: Skip already completed runs
- `--help`: Show help message

#### 3. View CPU Scaling Results

```bash
# View summary
cat ./cpu_scaling_500/cpu_scaling_summary.txt

# View comprehensive plot
open ./cpu_scaling_500/cpu_scaling_comprehensive.png
```

**Generated Plots:**
- `cpu_scaling_runtime.png` - Runtime vs CPU count
- `cpu_scaling_speedup.png` - Speedup analysis with ideal reference
- `cpu_scaling_efficiency.png` - Parallel efficiency curve
- `cpu_scaling_step_breakdown.png` - Per-step scaling
- `cpu_scaling_comprehensive.png` - 4-panel overview (publication-ready)

**Key Metrics:**
- **Speedup**: How much faster with more CPUs (ideal = linear)
- **Efficiency**: CPU utilization percentage (100% = perfect)
- **Sweet spot**: Optimal CPU count for efficiency/cost tradeoff

---

## Sample Size Benchmark Output Structure

```text
scalability_results/
├── timings.csv                      # Aggregated timing results
├── memory_usage.csv                 # Memory usage data
├── all_gvcfs.txt                    # List of all available GVCFs
├── hgc_scalability_report.txt       # Summary report
├── sample_sets/                     # Sampled GVCF lists per size
│   ├── samples_20.txt
│   ├── samples_50.txt
│   └── ...
├── run_logs/                        # Detailed logs per run
│   ├── run_20.log
│   ├── time_20.txt
│   └── ...
├── run_20/                          # Results for 20 samples
│   ├── combined_20.vds/             # Variant DataSet
│   ├── analysis_20.mt/              # MatrixTable
│   ├── sample_qc_20.ht/             # Sample QC table
│   ├── variant_qc_20.ht/            # Variant QC table
│   ├── cohort_20.vcf.bgz             # Final cohort VCF
│   ├── timing_20.json               # Timing breakdown
│   └── workflow_20.log              # Workflow log
├── run_50/
│   └── ...
└── plots/                           # Generated visualizations
    ├── hgc_scalability_runtime_breakdown.png
    ├── hgc_scalability_total_runtime.png
    ├── hgc_scalability_step_comparison.png
    ├── hgc_scalability_scaling_efficiency.png
    └── hgc_scalability_memory_usage.png
```

## Understanding the Results

### CSV Formats

**timings.csv** (column format; values are illustrative, not a recorded run):
```csv
sample_size,gvcf_combine_sec,vds_to_mt_sec,compute_qc_sec,mt_to_vcf_sec,total_sec
20,<seconds>,<seconds>,<seconds>,<seconds>,<seconds>
50,<seconds>,<seconds>,<seconds>,<seconds>,<seconds>
```

**memory_usage.csv** (column format):
```csv
sample_size,peak_memory_mb,run_time_sec
20,<mb>,<seconds>
50,<mb>,<seconds>
```

### Interpreting Plots

1. **Runtime Breakdown**: Shows which steps dominate processing time
2. **Total Runtime**: Linear trend line indicates O(n) scaling
3. **Step Comparison**: Identifies bottlenecks across sample sizes
4. **Scaling Efficiency**: Flat line = perfect linear scaling, increasing = superlinear
5. **Memory Usage**: Shows memory requirements for capacity planning

### Scaling Metrics

The summary report includes:
- **Linear fit**: `y = ax + b` where `a` is seconds added per sample
- **R-squared**: How well data fits linear model (1.0 = perfect)
- **Doubling analysis**: How runtime changes when doubling sample size
- **Memory per sample**: Average memory overhead per sample

## Workflow Details

### Step 1: GVCF Combination
- Uses Hail's `new_combiner()` to merge individual GVCFs
- Creates a Variant DataSet (VDS) optimized for storage
- Most computationally intensive step

### Step 2: VDS to MatrixTable
- Splits multi-allelic variants at the VDS level (the split already produces biallelic
  `GT`/`AD`; no separate LGT→GT downcoding is needed)
- Densifies the VDS to a dense MatrixTable
- Annotates adjusted genotypes (`adj`)
- Keys and sorts columns by sample ID

### Step 3: Compute QC
- Calculates sample-level QC metrics (call rate, Ti/Tv, etc.)
- Calculates variant-level QC metrics (allele frequency, HWE, etc.)
- Exports QC tables as Hail Tables (`.ht` format)

### Step 4: Export VCF
- Filters to adjusted genotypes
- Adds standard INFO fields (AF, AC, AN, call_rate)
- Exports clean multi-sample VCF (`.vcf.bgz`)

## Notes

### Memory Measurement
- On macOS: Uses `/usr/bin/time -l` to capture peak resident set size
- Reports maximum memory used by the Python process
- May not capture all Spark worker memory

### GVCF Sampling
- `benchmark.sh` samples each size with Python's `random.shuffle()`, seeded by `--seed`
  (default 42), so the same seed and GVCF list reproduce the same sample sets
- `../common/generate_sample_list.sh` (used to build the CPU-scaling benchmark's fixed
  list) instead uses `shuf` (Linux), `gshuf` (macOS, via `brew install coreutils`), or
  `sort -R` as a fallback

### No Cleanup
- All intermediate files are kept by default
- Allows for debugging and reproducibility
- Can be cleaned up manually after analysis

### Resume Capability
- Use `--resume` flag to skip already-completed runs
- Checks for existence of `timing_*.json` files
- Useful if benchmark is interrupted

## Example Workflow

```bash
# 1. Run benchmark (takes several hours for 7 sample sizes)
bash benchmark.sh --gvcf-dir /path/to/gvcfs

# 2. Generate plots
python plot_results.py --results-dir ./scalability_results

# 3. View results
cat scalability_results/hgc_scalability_report.txt
open scalability_results/plots/hgc_scalability_total_runtime.png

# 4. Clean up intermediate files (optional)
rm -rf scalability_results/run_*/combined_*.vds
rm -rf scalability_results/run_*/analysis_*.mt
```

## Troubleshooting

### Out of Memory
- Reduce sample sizes or test fewer sizes
- Increase Hail memory allocation (set `PYSPARK_SUBMIT_ARGS`)
- Use a machine with more RAM

### GVCF Files Not Found
- Check the `--gvcf-dir` path is correct
- Verify GVCF files have expected extensions (`.g.vcf.gz`, `.gvcf.gz`, etc.)
- Check file permissions

### Hail Initialization Errors
- Ensure Hail is properly installed: `pip install hail`
- Check Java is installed: `java -version`
- Review Hail logs in `run_*/hail_*.log`

### Resume Not Working
- Ensure `timing_*.json` files exist in run directories
- Use `--resume` flag explicitly
- Delete incomplete run directories to force re-run

## Advanced Usage

### Custom Sample Sizes

```bash
bash benchmark.sh \
    --gvcf-dir /path/to/gvcfs \
    --sample-sizes 10,25,50,100,200,400,800 \
    --output-dir ./custom_results
```

### Test on Subset First

```bash
# Quick test with small sizes
bash benchmark.sh \
    --gvcf-dir /path/to/gvcfs \
    --sample-sizes 20,50 \
    --output-dir ./test_run
```

### Different Reference Genome

```bash
bash benchmark.sh \
    --gvcf-dir /path/to/gvcfs \
    --reference GRCh37
```

## Performance Expectations

Runtime depends on hardware, storage, and GVCF size, so no fixed numbers are given
here. Run the smallest sample size first (20 samples) to get a real baseline for your
own setup, then use `hgc_scalability_report.txt`'s linear fit (`y = ax + b`) to
extrapolate to the larger sizes before committing to the full sweep — all 7 default
sizes together can take a long time.

## Citation

If you use these benchmarking scripts in your research, please cite:
- hvantk: <https://github.com/bigbio/hvantk>
- Hail: <https://hail.is/>

## Support

For issues or questions, open an issue on GitHub: <https://github.com/bigbio/hvantk/issues>

