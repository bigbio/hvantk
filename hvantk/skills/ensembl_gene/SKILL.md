# ensembl-gene

The `ensembl-gene` provider ships a single dataset, `ensembl-gene:structure` —
the canonical per-gene Ensembl table (gene ID, gene name, biotype,
coordinates, CDS length, coding-exon count, transcript count, and MANE
Select) parsed from the pinned-release GTF. It is the gene spine the
annotation pipeline builds on.

See `structure/SKILL.md` for the dataset's build, schema, and test details.

Upstream: https://www.ensembl.org/info/data/ftp/index.html
