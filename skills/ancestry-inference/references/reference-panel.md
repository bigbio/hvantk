# Building a 1000 Genomes reference

## The public files

- Genotypes: the 1000 Genomes high-coverage phased release, one bgzipped VCF with index
  per chromosome, GRCh38, 3202 samples. The directory and file pattern are recorded in
  `docs_site/guide/data-sources.md` (the `onek-genomes:variants` dataset), and the
  plugin's drift probe names the chromosome 22 file.
- Labels: the 2504-sample panel (`sample`, `pop`, `super_pop`, `gender`), downloaded and
  built by `hvantk reprocess onek-genomes:samples --raw-dir <dir> --output <samples.ht>`
  as a job (the build needs Hail; the download needs egress from wherever it runs).

Citation and reuse terms are in `THIRD_PARTY_DATA.md`.

## From a chromosome VCF to a reference MatrixTable

`templates/prepare-1kg-reference.py` does, as one job:

1. `hl.import_vcf` of the chromosome file (block-gzipped, so `force_bgz`), GRCh38.
2. Biallelic SNPs only; minor allele frequency above the threshold (default 0.05),
   computed from the genotypes, so related and unlabelled samples do not bias it after
   the join.
3. Thinning to a target count (default 20,000) by taking every k-th remaining variant in
   genomic order, which is deterministic and spreads the variants along the chromosome.
4. Join the sample table on `s`; samples without a label are dropped (the 698 related
   samples of the 3202 are not in the 2504 panel); `super_pop` becomes a flat column.
5. Write the MatrixTable.

Several chromosomes: run it per chromosome with the same thresholds and union the
results (`hl.MatrixTable.union_rows`) in a small job, or prepare from a few large
chromosomes; continental ancestry does not need the genome.

## Sizes

Chromosome 22 of the phased release is a file of about a gigabyte with about a million
variants; importing and filtering it takes a 16-core, 120 GB job under an hour. The
reference written from it is small (thousands of samples by 20,000 variants).
