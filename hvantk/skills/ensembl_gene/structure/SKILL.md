# ensembl-gene:structure

The canonical per-gene Ensembl table (location, name, biotype, and structure) parsed from the pinned GTF.

**Key:** `gene_id` (Ensembl, version suffix stripped)

| Column | Type | Meaning |
|---|---|---|
| `gene_id` | str | Ensembl gene ID, unversioned |
| `gene_name` | str | HGNC/Ensembl gene symbol |
| `chromosome` | str | Sequence name (e.g. `1`, `X`) |
| `gene_start` | int | Gene start (1-based) |
| `gene_end` | int | Gene end (1-based, inclusive) |
| `gene_biotype` | str | e.g. `protein_coding`, `lncRNA` |
| `mane_select` | str | MANE Select transcript ID, `""` if the gene has none |
| `cds_transcript` | str | Representative coding transcript: MANE if present, else longest CDS |
| `cds_length` | int | Summed CDS bp of `cds_transcript`; `0` for non-coding genes |
| `n_coding_exons` | int | CDS feature count of `cds_transcript` |
| `n_transcripts` | int | Distinct transcripts annotated for the gene |

**Release pin:** `hvantk/resources/ensembl_release.py`. The same pin governs the
PTM coordinate mapper, so both read one gene model.

**Params:** `protein_coding_only` (bool, default `False`).

**Why it exists:** `ensembl-gene:genes` carries coordinates only. CDS length is required
to normalise PTM site density, and gene length mechanically confounds rare-variant burden
counts, so it must be available as a nuisance covariate.
