# Sizing a joint-genotyping run

Memory is set by the combine stage (samples times interval size); time by the number of
partitions and the node's cores. Starting points for one node of a CPU partition:

| samples | layout | cores | memory | time | notes |
| --- | --- | --- | --- | --- | --- |
| up to 10 | one job | 8 | 32 GB | 30 min | the eval runs three samples on one contig slice in a few minutes |
| up to 100 | one job | 16 | 120 GB | 4 h | genome-wide; `--n-partitions` 200 |
| up to 500 | per contig | 16 | 100 GB | 2 h per contig | an array over contigs, 6 to 8 concurrent |
| 1000 and more | per contig | 32 | 160 GB | 2 h per contig | chromosome 1 and 2 need the time; small contigs finish in minutes |

Rules:

- Double the memory before adding cores; the combiner's memory is per node, not per
  core, and a cgroup kill wastes the whole stage.
- Time is elapsed plus half; a timeout at the combine stage loses it entirely.
- Interval size: the default genome intervals are right for whole genomes; exome data
  uses `--use-exome-default-intervals`; a combine that runs out of memory can halve the
  interval size at the cost of more partitions.
- Dense MatrixTable partitions: aim for a few hundred MB per part file; set
  `--n-partitions` to about two to four times the cores for a single job, more per
  contig for a large cohort.
- Per contig, the dense-MT partition count can be derived from the VDS partition count
  divided by ten, with a floor of fifty, which is what large runs used.

What the eval's three-sample slice needs: the first budget row. Its combine finishes in
under two minutes; the rest is JVM start and the VCF export.
