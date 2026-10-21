# What a run writes

The chain the template runs (`gvcf-combine`, `vds2mt`, `compute-qc`, `mt2vcf`) writes,
under the output directory with the prefix you chose:

| path | what it is |
| --- | --- |
| `<prefix>.vds/` | the lossless VDS (`reference_data/`, `variant_data/`, each with `_SUCCESS`) |
| `<prefix>.mt/` | the dense MatrixTable, columns keyed and sorted by sample |
| `qc/<prefix>_sample_qc.csv` | per-sample QC, full, comma-separated |
| `qc/<prefix>_variant_qc.tsv` | per-variant QC, full, tab-separated; star alleles dropped |
| `<prefix>.vcf.bgz` | biallelic cohort VCF, `adj` filtered, `min_ac` as given; no index |
| `recipe.txt` | the four commands, for resuming or for a recipient |

The orchestrated `hvantk hgc pipeline` writes the layout below when it runs; with
`<prefix>` from `--output-prefix` (default `cohort`):

| path | what it is | notes |
| --- | --- | --- |
| `<prefix>.vds/` | the lossless VDS: `reference_data/` and `variant_data/`, each with `_SUCCESS` | the artifact to archive; a subset of it equals a fresh combine |
| `<prefix>.mt/` | the dense MatrixTable, columns keyed and sorted by sample | rows include gVCF placeholder alleles until QC removes them |
| `<prefix>_filtered.mt/` | only with `--apply-qc-filters`: the MatrixTable the VCF is exported from | here even when `--mt-path` points elsewhere |
| `qc/<prefix>_sample_qc.ht`, `.csv` | per-sample QC; the CSV is the full table | one row per sample |
| `qc/<prefix>_variant_qc.ht`, `.tsv` | per-variant QC; the TSV, tab-separated, holds the first 10,000 variants | the `.ht` is complete |
| `<prefix>.vcf.bgz` | biallelic cohort VCF, `adj` filtered, `min_ac` 1, multiallelics split | no index; Hail writes none |
| `qc/<prefix>_qc_report.html` | only with `--generate-qc-report` | |
| `.pipeline_state.json` | stage start, end and output paths | read it before resuming |
| `logs/pipeline_<timestamp>.log` | the run's log from every HGC module, at the console's level: warnings by default, progress with `hvantk -v` | the Slurm log is elsewhere |

Resuming: every stage has a `--skip-*` flag and the two artifacts have `--vds-path` and
`--mt-path` overrides; the state file says which stages completed.

Standalone commands write differently: `compute-qc -o <dir>` writes
`<prefix>_sample_qc.csv`, `<prefix>_variant_qc.tsv` (full, tab-separated) and
`<prefix>_with_qc.mt`; `qc-summary` reads a QC directory and writes markdown, CSV or JSON.

Counting: the MatrixTable row count is not the variant count while placeholder alleles
are present; the QC tables and the VCF drop them. Record the MatrixTable row count as
produced by the pipeline, the VCF record count, and the sample set, so a run is
compared with a run, not with a notion of "variants".
