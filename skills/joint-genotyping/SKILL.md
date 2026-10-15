---
name: joint-genotyping
description: 'Use when a set of single-sample gVCFs must become a joint callset with hvantk HGC pipeline on a Slurm cluster: combining gVCFs into a VDS, densifying to a MatrixTable, sample and variant QC, and exporting a cohort VCF, for a handful of samples or for thousands per contig. Also when a combine job runs out of memory or time, when counts look inflated by placeholder alleles, when the exported VCF has no index, or when contig names do not match the build.'
user-invocable: true
argument-hint: "<gVCF directory> <output directory> [--prefix NAME] [--build GRCh38|GRCh37] [--per-contig]"
metadata:
  stage: cohort
  status: provisional
  scale: node
---

# Joint genotyping with the HGC pipeline

Read `skills/_conventions/SKILL.md` first; this skill assumes its site-profile rule,
its section contract and its measurement contract. Jobs run through `hpc-hail-job`.

## 1. Scope and outcome

From one gVCF per sample, bgzipped and indexed, to a cohort callset: the lossless
VDS (`<prefix>.vds`), the dense MatrixTable (`<prefix>.mt`), per-sample and per-variant
QC tables (`qc/`), and a biallelic cohort VCF (`<prefix>.vcf.bgz`), all under one output
directory with a state file and logs. The outcome also includes the counts a recipient
needs: samples, variant rows, VCF records, and where the QC summary is.

Use it for a cohort of any size. A few samples run as one job on one node; hundreds or
thousands run one job per contig (section 4), then the per-contig outputs are kept as
the callset. Not this skill: subsetting an existing callset (`callset-subset-export`),
variant calling itself (gVCFs are the input), and multi-node Spark.

## 2. Preflight

1. Site profile: `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`;
   stop and ask for any missing value. Every later command that reads a profile value
   starts with the same `source "$HVANTK_SITE_PROFILE" &&`.
2. Checkout and image as in `hpc-hail-job` (commit, digest, run notes).
3. Inputs: one gVCF per sample, `*.g.vcf.gz` with a `.tbi` next to each, readable from
   compute nodes. Check the contig names in a header (`zcat <file> | grep -m1 '^##contig'`)
   against the build you will pass: `chr1` for GRCh38, `1` for GRCh37; HGC does no
   recoding. Count the files: that is the sample count every size decision uses.
4. Output: a directory under `$HVANTK_WORK`; the pipeline refuses to overwrite unless
   told. Never under the home directory.
5. Size the job from the sample count with `references/sizing.md`; for more than a few
   hundred samples plan per contig from the start.

## 3. Procedure

1. Dry-run the plan: the pipeline prints its stages and paths with `--dry-run` and
   exits. Run it through the job template so the paths resolve on a compute node:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/joint-genotyping/templates/pipeline-job.sh" \
     --gvcf-dir "$HVANTK_WORK/gvcfs" --out "$HVANTK_WORK/callsets/<name>" --prefix <name> \
     --build GRCh38 --cpus 8 --mem-gb 32 --time 0-00:30:00 --name <run>-plan --dry-run-pipeline
   ```
2. Submit the real run with the same template without `--dry-run-pipeline`; add
   `--wait` to block, or run it in the background and poll the job id (`hpc-hail-job`
   section 3). The template passes `--reference-genome`, `--output-prefix` and, when
   you give it, `--n-partitions` for the dense MatrixTable.
3. Read `<out>/.pipeline_state.json` and `<out>/logs/`: each stage records its start,
   end and output path. A stage that failed leaves the earlier outputs in place; resume
   with `--skip-combine-gvcfs --vds-path <out>/<prefix>.vds` (and the other `--skip-*`
   flags the state file justifies) rather than starting over.
4. QC: `hvantk hgc qc-summary -d <out>/qc -o <out>/qc/summary.md` as a small job, and
   `hvantk hgc qc-report -i <out>/<prefix>.mt -o <out>/qc/report.html` when a report is
   wanted. `references/outputs.md` says what each QC file holds and what it leaves out.
5. Index the VCF if the site has `bcftools` or `tabix` (Hail writes none): a cluster
   module or a user install; record which. Without one, say so in the run notes.
6. Verify (section 6) and write the run notes with the counts.

## 4. Scale-up

- Per contig. Split the gVCFs by contig (one directory per contig, each file indexed),
  then `templates/per-contig-array.sbatch` runs `hvantk hgc gvcf-combine` as one array
  task per contig with the same container environment as a single job, followed by
  `vds2mt` and `mt2vcf` per contig. Keep the per-contig outputs; a genome-wide VCF, if
  needed, is a `bcftools concat` of the per-contig files in contig order.
- Memory grows with the sample count and the interval size; time grows with both. Start
  from the table in `references/sizing.md` and double memory before adding cores.
- `--n-partitions` on the dense MatrixTable reduces part files on the parallel
  filesystem: aim for a few hundred MB per part.
- Interval size (`--import-interval-size`) trades combiner memory against partitions;
  the default genome intervals are right unless the combine runs out of memory.
- Intermediates stay (`--keep-intermediates` is always on); the VDS is the lossless
  record and the thing to archive, not the VCF.

## 5. Traps

- The combine fails on a contig it cannot find. Cause: gVCF contigs named for the other
  build. Fix: pass the build the files use; never recode gVCFs by hand.
- Row counts look far too high. Cause: gVCF placeholder alleles (`<*>`, `<NON_REF>`,
  `*`) are rows in the dense MatrixTable. Fix: count variants from the QC tables or the
  VCF, which drop them; `compute-qc` drops star alleles by default.
- The variant QC CSV has 10,000 rows for a cohort with millions. Cause: the pipeline
  caps that export. Fix: the full table is `qc/<prefix>_variant_qc.ht`; the standalone
  `compute-qc` writes a full `.tsv`.
- `mt2vcf` refuses with an `adj` error. Cause: the `adj` entry field is missing or
  mostly missing, so filtering on it would delete genotypes. Fix: run `vds2mt` with
  genotype adjustment on (the default), or export with `--no-filter-adj` deliberately.
- The exported file ends in `.vcf.gz` and tools complain. Cause: Hail writes block-gzip;
  the toolkit's extension is `.vcf.bgz`. Fix: name it `.vcf.bgz`, index it yourself.
- No `.tbi` next to the VCF. Cause: Hail writes none. Fix: section 3 step 5.
- `seff` reports more memory than requested on a job that completed. Cause: accounting
  sums per-process RSS. Fix: trust the request and the cgroup (`hpc-hail-job`).
- A combine of thousands of samples times out at the default interval size. Cause:
  too few, too large partitions. Fix: per contig, smaller intervals, more time.
- The per-contig VCFs disagree on sample order. Cause: columns keyed but not sorted in
  older versions. Fix: current `vds2mt` sorts columns; verify with the header before
  concatenating.
- A job writes under the home directory and dies on quota. Cause: a relative output
  path; the job's working directory is the log directory. Fix: absolute paths under
  the work root.

## 6. Verification

1. `<out>/<prefix>.vds/reference_data/_SUCCESS`, `<out>/<prefix>.vds/variant_data/_SUCCESS`
   and `<out>/<prefix>.mt/_SUCCESS` exist; `<out>/qc/<prefix>_sample_qc.csv` has one row
   per sample; `<out>/<prefix>.vcf.bgz` exists.
2. Counts, through the image as a small job (`templates/pipeline-job.sh --describe <out> --prefix <name>`):
   the sample set of the MatrixTable, its row count, and the VCF record count. Put all
   three in the run notes; a recipient compares against them.
3. The VCF header's sample columns equal the gVCF sample names (suffixes included).
4. The run notes carry commit, image digest, every job id, the build, and the counts.
