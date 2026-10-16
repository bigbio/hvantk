---
name: callset-subset-export
description: 'Use when a sample subset of a joint callset must be delivered as a verified VCF: subsetting the lossless VDS by sample list, densifying and exporting with hvantk, verifying header samples and record counts, writing checksums and a manifest, per contig for large callsets. Also when sample names in a list do not match the callset, when per-contig VCFs refuse to concatenate, or when a delivery must go to a read-only share.'
user-invocable: true
argument-hint: "<callset VDS> <sample list> <output directory> [--prefix NAME] [--per-contig]"
metadata:
  stage: cohort
  status: provisional
  scale: node
---

# Subset a joint callset and deliver a verified VCF

Read the `hvantk:workflow-conventions` skill (`skills/workflow-conventions/SKILL.md`)
first; this skill assumes its site-profile rule, its section contract and its
measurement contract. Jobs run through `hpc-hail-job`.

## 1. Scope and outcome

A recipient asks for some samples of a joint callset. The outcome is a directory with a
biallelic VCF of exactly those samples (`<prefix>.vcf.bgz`), the subset's own VDS and
MatrixTable, a `SHA256SUMS` file that validates, and a `README.md` naming the source
callset, the commit and image that produced the files, the sample list, the commands,
and the counts. The subset is taken from the lossless VDS, so it equals a fresh joint
call of those samples: allele counts and sites are recomputed for the subset, not
inherited from the full cohort.

Use it for any subset size. Not this skill: making the callset (`joint-genotyping`),
subsetting a cohort VCF with bcftools (that keeps the cohort's sites and stale counts),
and sharing controlled-access data, which is a human decision this skill never makes.

## 2. Preflight

1. Site profile: `source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/_site/check-site-profile.sh"`;
   stop and ask for any missing value. Every later command that reads a profile value
   starts with the same `source "$HVANTK_SITE_PROFILE" &&`.
2. Checkout and image as in `hpc-hail-job` (commit, digest, run notes).
3. The callset's own sample spelling: list it as a small job,
   `templates/subset-job.sh --list-samples <callset.vds>`; callsets often carry
   suffixes or doubled identifiers that a human list lacks. Map the requested list onto
   the callset's spelling and keep both lists in the delivery; never guess a mapping.
4. Output under `$HVANTK_WORK`; the delivery target, if it is a read-only mount on
   compute nodes (`$HVANTK_READONLY_MOUNTS`), receives a copy from the login node at the
   end, never a job's direct write.
5. Size from the callset: a subset of a few samples on one contig slice runs in minutes
   on 8 cores; a genome-wide callset of hundreds of samples goes per contig (section 4).

## 3. Procedure

1. Submit the chain as four jobs from the login node with `templates/subset-job.sh`:
   ```bash
   source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/callset-subset-export/templates/subset-job.sh" \
     --vds "$HVANTK_WORK/callsets/<name>/<name>.vds" --samples "$HVANTK_WORK/deliveries/<id>/samples.txt" \
     --out "$HVANTK_WORK/deliveries/<id>" --prefix <id> --cpus 8 --mem-gb 32 --time 0-00:30:00 --name <run>-subset
   ```
   The chain is: `subset-vds.py` (keeps the listed samples, drops alleles no listed
   sample carries, refuses to overwrite, verifies the written sample set), `hvantk hgc
   vds2mt`, `hvantk hgc mt2vcf` (biallelic, `adj`-filtered, `min_ac` 1), then
   `verify-export.py`, which checks the header samples against the list, counts
   records, writes `SHA256SUMS` and the `README.md`. Each step waits for the previous
   one and stops the chain on failure.
2. Index the VCF when the site has `bcftools` or `tabix` (a module or a user install);
   Hail writes no index. Record which tool, or that none was available; a `.tbi` added
   afterwards is not in `SHA256SUMS` unless you append it.
3. Read `<out>/verify.json` (`references/verification.md` lists the checks): the sample set, the record count and the checks, which the
   README repeats for the recipient.
4. Deliver: copy the directory to the target from the login node, throttled when the
   target is a shared mount, then validate the checksums at the target
   (`sha256sum -c SHA256SUMS`) before telling anyone it is there.
5. Write the run notes with the commit, the image digest, every job id, the sample
   list's path, and the counts.

## 4. Scale-up

- Per contig for large callsets: `templates/per-contig-subset-array.sbatch` runs the
  subset, densify and export for one contig per array task, with the same container
  environment as a single job; a genome-wide VCF is a `bcftools concat` in contig
  order, re-indexed and re-checksummed.
- Dense MatrixTable partitions: derive from the VDS partition count divided by ten with
  a floor of fifty, or pass `--n-partitions` explicitly.
- Memory is per node and grows with the subset's size times the sites; the subset step
  is the heavy one. Thirty-two cores and 160 GB handle a chromosome of a thousand-sample
  callset in minutes for a small subset.
- Copying a VDS or a MatrixTable tree moves millions of part files; copy with a
  parallel rsync per tree, or deliver only the VCF, which is one file.

## 5. Traps

- The subset refuses the list: samples "not in the callset". Cause: spelling; the
  callset's names carry suffixes or doubled identifiers. Fix: preflight step 3.
- A `zcat ... | head` pipeline exits 141 under `set -o pipefail`. Cause: SIGPIPE. Fix:
  `{ zcat f || true; } | head`.
- `bcftools concat --naive` refuses per-contig files. Cause: headers differ, often
  only by a `bcftools_viewCommand` line. Fix: full `concat` (recompresses), or make the
  headers identical before the naive path.
- `bcftools index -n` cannot find the index of a `.vcf.bgz`. Cause: the tool expects
  `.vcf.gz` naming. Fix: pass the index explicitly (`file##idx##file.tbi`).
- The exported VCF has no `.tbi`. Cause: Hail writes none. Fix: section 3 step 2.
- `seff` says the subset used far more memory than requested. Cause: accounting sums
  shared pages. Fix: trust the request and the cgroup.
- Per-contig VCFs disagree on sample order. Cause: columns not sorted in older
  versions. Fix: current `vds2mt` sorts columns; check headers before concatenating.
- The delivery copy fails with "read-only file system" from a job. Cause: the target is
  read-only on compute nodes. Fix: copy from the login node (section 3 step 4).
- The VCF is a subset of the cohort VCF rather than of the VDS. Cause: the shortcut.
  Fix: subset the VDS; a cohort-VCF subset keeps sites and counts of the full cohort.

## 6. Verification

1. `<out>/verify.json`: `samples_match` true, `n_samples` equal to the list, `records`
   the count, `biallelic` true.
2. `sha256sum -c <out>/SHA256SUMS` passes where the files are, and again at the
   delivery target.
3. The README names the source callset, the commit and image digest, the sample list
   (both spellings), the commands, and the counts; a recipient can repeat the subset.
4. The run notes carry the commit, image digest, every job id and the counts.
