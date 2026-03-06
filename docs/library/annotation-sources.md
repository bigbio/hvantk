# hvantk

Hail-based multiomics variant annotation toolkit.

List of raw sources and URL to download the data.

## File format note

Hail supports standard gzip (`.gz`) and uncompressed files but processes them single-threaded. **Block gzip (BGZF)** compressed files (`.bgz`) enable parallel import and are strongly recommended. If your downloaded files are `.gz`, either:

1. Use `--auto-convert-bgz` in supported commands to convert on-the-fly
2. Run `hvantk convert-bgz input.gz` to convert before import

See the [Usage Guide](usage.md#file-format-conversion) for details.

For download instructions (built-in downloaders and manual download steps), see the [Data Acquisition Guide](data-acquisition.md).

## Annotation sources

- Variants and genomic regions

  - **Missense variants prediction scores (from dbNSFP v4.9a)**

    Description: A database of functional predictions scores for human missense variants.\
    URL: https://sites.google.com/site/jpopgen/dbNSFP

  - **ClinVar annotations VCF (GRCh38)**

    Description: A database of clinically relevant variants and their annotations (e.g. Pathogenic, Benign, VUS).\
    URL: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/

  - **gnomAD annotations (gene level)**

    Description: Gene level constraint metric annotations (e.g. lof and missense) from the Genome Aggregation Database (gnomAD).\
    URL: https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv

  - **Protein-protein interaction site (INSIDER)**

    Description: Protein-protein interaction sites from INSIDER database.\
    URL: http://interactomeinsider.yulab.org/downloads.html

  - **Ensemble gene annotations**

    Description: Ensemble gene annotations (e.g. gene name, gene id, biotype, transcript ID).\
    URL: https://www.ensembl.org/info/data/ftp/index.html

  - **GeVIR score**

    Description: Gene variation intolerance rank.\
    URL: https://www.nature.com/articles/s41588-019-0560-2

  - Coding-constrained region (CCR) score

    Description:\*\* highly constrained coding regions (CCRs) in the human genome.\
    URL: https://www.nature.com/articles/s41588-018-0294-6

- Bulk RNA-seq data

  - **Human tissue expression E-MTAB-6814**

    Description: Human tissue gene expression (brain, heart, liver, kidney), multiple developmental time points.\
    URL: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6814

- Single-cell RNA-seq data

  - **Human heart single-cell RNA-seq data (Asp 2019)**

    Description: Embryonic human heart single-cell RNA-seq data 6.5 wpc (PMID:31835037).\
    URL: https://data.mendeley.com/datasets/mbvhhf8m62/2

  - **Human heart single-cell RNA-seq data (Farah 2024)**

    Description: Single-cell RNA-seq data of the developing human heart, 9-15 wpc (PMID:31835037).\
    URL: https://cells.ucsc.edu/?bp=heart&ds=hoc

  - **Human heart cell atlas (HCA)**
    Description: Adult human heart cell atlas (https://doi.org/10.1038/s41586-020-2797-4).\
    URL: https://cells.ucsc.edu/?bp=heart&ds=heart-cell-atlas

