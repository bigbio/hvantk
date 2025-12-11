# HGC Scalability Benchmark

This directory contains scripts to benchmark the scalability of the HGC (Hail-based Genotype Combiner) workflow at different sample sizes.

## Overview

The benchmark measures the complete end-to-end HGC workflow performance:

1. **GVCF → VDS**: Combine individual GVCF files into a Variant DataSet
2. **VDS → MatrixTable**: Convert VDS to analysis-ready MatrixTable
3. **Compute QC**: Calculate sample and variant quality control metrics
4. **Export VCF**: Generate final cohort VCF file

### Metrics Collected

- **Runtime**: Wall-clock time for each step and total workflow
- **Memory**: Peak memory usage per run
- **Scalability**: How performance scales with sample size

## Files

```
analysis/
├── hgc_scalability_benchmark.sh       # Main orchestration script (bash)
├── hgc_scalability_benchmark.py       # Python workflow runner
├── plot_scalability_results.py        # Results visualization script
└── README_scalability.md              # This file
```

## Requirements

- Python 3.8+
- hvantk library installed
- Hail installed and configured
- Required Python packages: matplotlib, seaborn, pandas, numpy, scipy

## Usage

### 1. Run the Benchmark

```bash
cd analysis

# Run with default settings
bash hgc_scalability_benchmark.sh

# Or customize parameters
bash hgc_scalability_benchmark.sh \
    --gvcf-dir /path/to/gvcfs \
    --output-dir ./my_results \
    --sample-sizes 20,50,100,250,500,750,1000 \
    --reference GRCh38 \
    --seed 42
```

**Default Configuration:**
- GVCF directory: `/mnt/nfs/KOL_UOL/projects/CHD_1000WGS/variant_calling/split_vcfs/chr20`
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

### 2. Generate Plots

After the benchmark completes:

```bash
python plot_scalability_results.py --results-dir ./scalability_results
```

This generates:
- Runtime breakdown (stacked bar chart)
- Total runtime vs sample size (with trend line)
- Step-by-step comparison (line plot)
- Scaling efficiency (time per sample)
- Memory usage vs sample size
- Summary report (text file)

## Output Structure

```
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
│   ├── cohort_20.vcf.gz             # Final cohort VCF
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

**timings.csv:**
```csv
sample_size,gvcf_combine_sec,vds_to_mt_sec,compute_qc_sec,mt_to_vcf_sec,total_sec
20,45.2,12.3,8.5,5.1,71.1
50,89.5,25.7,18.2,9.3,142.7
...
```

**memory_usage.csv:**
```csv
sample_size,peak_memory_mb,run_time_sec
20,2048.5,71
50,3072.8,143
...
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
- Converts VDS to dense MatrixTable format
- Splits multi-allelic variants
- Annotates adjusted genotypes
- Converts LGT to GT

### Step 3: Compute QC
- Calculates sample-level QC metrics (call rate, Ti/Tv, etc.)
- Calculates variant-level QC metrics (allele frequency, HWE, etc.)
- Exports QC tables as Hail Tables (`.ht` format)

### Step 4: Export VCF
- Filters to adjusted genotypes
- Adds standard INFO fields (AF, AC, AN)
- Exports clean multi-sample VCF

## Notes

### Memory Measurement
- On macOS: Uses `/usr/bin/time -l` to capture peak resident set size
- Reports maximum memory used by the Python process
- May not capture all Spark worker memory

### GVCF Sampling
- Reproducible with fixed random seed
- Uses `shuf` (Linux) or `gshuf` (macOS) for random sampling
- Falls back to `sort -R` if neither available

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
bash hgc_scalability_benchmark.sh

# 2. Generate plots
python plot_scalability_results.py --results-dir ./scalability_results

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
bash hgc_scalability_benchmark.sh \
    --sample-sizes 10,25,50,100,200,400,800 \
    --output-dir ./custom_results
```

### Test on Subset First

```bash
# Quick test with small sizes
bash hgc_scalability_benchmark.sh \
    --sample-sizes 20,50 \
    --output-dir ./test_run
```

### Different Reference Genome

```bash
bash hgc_scalability_benchmark.sh \
    --reference GRCh37
```

## Performance Expectations

Typical runtime (approximate, varies by hardware):
- 20 samples: ~1 hour
- 50 samples: ~2 hours
- 100 samples: ~4 hours
- 250 samples: ~10 hours
- 500 samples: ~20 hours
- 750 samples: ~30 hours
- 1000 samples: ~40 hours

**Full benchmark**: ~100+ hours total

## Citation

If you use these benchmarking scripts in your research, please cite:
- hvantk: https://github.com/bigbio/hvantk
- Hail: https://hail.is/

## Support

For issues or questions:
- Open an issue on GitHub: https://github.com/bigbio/hvantk/issues
- Contact the hvantk development team

---

**Last Updated**: December 11, 2025

