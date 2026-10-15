# Verifying a delivered VCF

What `templates/verify-export.py` checks and writes, and how to repeat it by hand.

| check | how | passes when |
| --- | --- | --- |
| header samples | the `#CHROM` line's columns from the tenth on | equal to the sample list, as a set and in sorted order |
| record count | lines not starting with `#` | the number the README states; equal across copies |
| biallelic | `ALT` has no comma | every record |
| checksums | `SHA256SUMS` over every visible file in the delivery except itself, README included; a `.tbi` added later must be appended | `sha256sum -c SHA256SUMS` passes |

By hand, with bcftools where the site has it:

```bash
bcftools query -l <prefix>.vcf.bgz | sort | diff - <(sort samples_callset.txt)   # no output
bcftools index -t <prefix>.vcf.bgz && bcftools index -n '<prefix>.vcf.bgz##idx##<prefix>.vcf.bgz.tbi'
sha256sum -c SHA256SUMS
```

Without bcftools, the same through the image:

```bash
source "$HVANTK_SITE_PROFILE" && bash "$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh" --name <run>-verify --cpus 2 --mem-gb 8 --time <enough to read the whole VCF once> --wait -- \
  python "$HVANTK_REPO/skills/callset-subset-export/templates/verify-export.py" <out> --prefix <prefix> --samples <list>
```

The README that `verify-export.py` writes carries: the source VDS path, the commit and
image digest given on the command line, the sample list in both spellings, the exact
commands of the chain, the counts, and the checksum file's name. A recipient who has the
README and the files can repeat the subset and compare.
