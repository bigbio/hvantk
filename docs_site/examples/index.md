# Examples

Tutorials and walkthroughs for hvantk tools and workflows. Each page covers concepts, CLI and Python API usage, and expected outputs.

Runnable scripts are available in the [`examples/`](https://github.com/bigbio/hvantk/tree/main/examples) directory on GitHub.

## Tools

| Example | Description |
|---------|-------------|
| [HGC](hgc.md) | Joint genotyping QC and scalability benchmarks |
| [EnrichEx](enrichex.md) | Gene set overlap enrichment and burden testing |
| [PS-ROC](psroc.md) | Pathogenicity score ROC evaluation |
| [Ancestry](ancestry.md) | Ancestry inference with PCA + Random Forest |
| [QTL Cascade](qtlcascade.md) | Molecular QTL integration (eQTL + pQTL) |
| [PTM](ptm.md) | Post-translational modification site analysis |
| [Re-rank](rerank.md) | Re-rank a gene burden by multi-omic credibility, on synthetic fixtures |

## Data Sources

| Example | Description |
|---------|-------------|
| [ClinVar Streaming](clinvar.md) | Stream and filter ClinVar variant annotations |
| [ClinGen](clingen.md) | Gene-disease validity queries and ontology categorization |


## Runnable examples in the repository

Every example above also ships as a runnable directory under
[`examples/`](https://github.com/bigbio/hvantk/tree/main/examples), with fixtures and its
own `README.md`. One has no page on this site:

| Example | What it runs | Needs |
|---------|--------------|-------|
| [1k_genome](https://github.com/bigbio/hvantk/tree/main/examples/1k_genome) | Building a Hail MatrixTable from 1000 Genomes NYGC high-coverage data, for use as an ancestry reference panel | Hail/Spark; large download |

Note that the runnable examples need more than the base install:

- `rerank` wants `--extras ml` (scikit-learn), and starts no Hail/Spark runtime.
- `psroc` and `ancestry` want `--extras psroc` / `--extras ancestry` (scikit-learn, plus
  plotting libraries).
- `enrichex` and the `hgc` benchmarks both plot, so they want `--extras enrichex` /
  `--extras hgc`.
- The `enrichex` burden scripts, `hgc`, `psroc`, `ancestry`, `clinvar`, `clingen` and
  `1k_genome` start a Hail/Spark runtime (the `enrichex` overlap scripts do not).
  `clinvar`, `clingen` and `1k_genome` also need network access, or a local copy of the
  raw data with the plugin's `--skip-download`.
