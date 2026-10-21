# Third-Party Data Attribution

This document lists third-party datasets committed to this repository as test fixtures. These are real excerpts from public data sources, committed to support reproducible testing. Synthetic fixtures (fabricated values for testing purposes) are not listed here. Test data is excluded from built packages via the `exclude` list in `pyproject.toml`, so this document concerns the git repository only. A pull request that commits a real excerpt as a fixture adds its entry here.

## Expression Atlas

### Fixture: `hvantk/skills/expression_atlas/tests/testdata/raw/expression-atlas/`

**Source**: Expression Atlas experiment E-MTAB-6798, "Mouse RNA-seq time-series of the development of seven major organs"

**Version**: no upstream release version; derived from upstream files retrieved before 2025-05-12 (a full copy was committed in bf2b310c and later removed as unused).

**Files**:
- `E-MTAB-6798-transcripts-tpms.tsv`: Transcript-level TPM expression matrix
- `E-MTAB-6798.condensed-sdrf.tsv`: Sample and data relationship format metadata

**Licence**: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

**Attribution**: Expression Atlas, EMBL-EBI — https://www.ebi.ac.uk/gxa/experiments/E-MTAB-6798. Citation: "Expression Atlas in 2026: enabling FAIR and open expression data through community collaboration and integration" (*Nucleic Acids Research*, 10 December 2025).

**Modifications**:
- Expression matrix: first 21 transcript rows (20 distinct genes; one gene's second transcript retained to exercise the fixture's multi-transcript-per-gene shape)
- Samples: first 4 of 317 samples (ERR2588382, ERR2588384, ERR2588383, ERR2588399 in header order)
- SDRF: lines for the 4 retained samples only, preserving the original ragged field layout

## AlphaGenome

### Fixture: `hvantk/skills/alphagenome/tests/testdata/raw/alphagenome/`

**Source**: AlphaGenome (Google DeepMind) variant-effect predictions for 411 ClinVar canonical-splice variants, from an AlphaGenome SDK `score_variant` → `variant_scorers.tidy_scores()` run on 2026-06-24 (RECOMMENDED scorers, 1 Mb interval, organism human, API backend).

**Version**: the AlphaGenome SDK and hosted model as of that run; the exact SDK version was not recorded (the run's manifest records only the backend, scorer set, interval and organism). Source shard `scores.shard-00000-of-00001.parquet`, 16,642,316 rows, sha256 `261cf87273f1b36d20db640305f414b2e9bdc55824cddff9e4cc4bcc4c6745a1`.

**Files**:
- `clinvar-subset.parquet`: AlphaGenome scores in the `tidy_scores()` long format
- `NOTICE.md`: the AlphaGenome Output Terms notice that must stay with the data
- `make_fixture.py`: the script that selected the rows (code, not data)

**Other copies of these outputs** (same source and terms):
- `hvantk/skills/alphagenome/tests/snapshots/sample_rows.json`: the builder's summary of the fixture, with scores and track names verbatim; its notice is `tests/snapshots/NOTICE.md` beside it
- `hvantk/skills/alphagenome/tests/test_builder.py`: the oracle constants `AG_SPLICE_SITES` and `HEART_SPLICE_SITE_USAGE`, marked by a comment

**Licence**: AlphaGenome Output Terms of Use (https://deepmind.google.com/science/alphagenome/output-terms): non-commercial use only, and the outputs must not be used to train machine-learning models. Not covered by this repository's MIT licence; the conspicuous notice the terms require is `NOTICE.md` in the fixture directory, with a second one beside the snapshot.

**Attribution**: Google DeepMind, AlphaGenome. ClinVar variant identities are public domain.

**Modifications**:
- Rows only: 3 of the 411 variants (`chr3:39408741:T>C`, `chr6:112216367:C>A`, `chrX:153694448:T>G`) and 150 of the 16,642,316 rows, selected by `make_fixture.py` with seed 20261006 and written in reverse source order. For each variant and scorer it keeps the row with the largest |quantile_score|, the row with the largest |raw_score| and one random row, plus every heart SPLICE_SITE_USAGE row
- Columns: unchanged

## dbNSFP

### Fixture: `hvantk/tests/testdata/raw/dbnsfp/`

**Source**: dbNSFP, academic branch, release v4.9a: functional predictions and annotations for all potential human non-synonymous and splice-site SNVs, distributed through https://www.dbnsfp.org after academic registration.

**Version**: v4.9a (academic). Committed to this repository in efd6c3bd (2025-08-21).

**Files**:
- `dbNSFP4_v49a_example_variants.bgz`: rows of the dbNSFP variant table in its distributed tab-separated layout, BGZF-compressed
- `NOTICE.md`: the licence notice that stays with the data

**Other copies of these data** (same source and terms):
- `hvantk/skills/dbnsfp/tests/snapshots/sample_rows.json`: six fixture rows as the dbNSFP builder outputs them, with the values parsed into typed fields; its notice is `tests/snapshots/NOTICE.md` beside it. `schema.json` there holds only field names and types.

**Licence**: CC BY-NC-ND 4.0 (https://creativecommons.org/licenses/by-nc-nd/4.0/), the licence of the dbNSFP academic branch (https://www.dbnsfp.org/license/): non-commercial use only, with attribution, and no sharing of modified versions. Commercial use requires dbNSFP's paid commercial licence, and the CADD, VEST, M-CAP, MutScore, PolyPhen-2, PrimateAI and RGC Million Exome scores in the academic branch also need commercial licences from their authors. Not covered by this repository's MIT licence; the notice is `NOTICE.md` in the fixture directory, with a second one beside the snapshot.

**Attribution**: dbNSFP, © 2024–2026 Genos Bioinformatics LLC — https://www.dbnsfp.org. Citation: Liu X, Li C, Mou C, Dong Y, Tu Y. "dbNSFP v4: a comprehensive database of transcript-specific functional predictions and annotations for human nonsynonymous and splice-site SNVs." *Genome Medicine* 12, 103 (2020), https://doi.org/10.1186/s13073-020-00803-9.

**Modifications**:
- None to the content: the header (all 458 columns) and the first 4,999 variants on chromosome 10 of the v4.9a variant table, byte-identical to the release once decompressed, stored BGZF-compressed
- To regenerate it, extract whole lines unedited: trimming columns or editing values would make a modified version, which the licence does not allow sharing

## ClinVar

### Fixture: `hvantk/skills/clinvar/tests/testdata/raw/clinvar/clinvar_20220403_chr20.vcf.bgz`

**Source**: ClinVar, NCBI/NLM — public archive of relationships among human variants and phenotypes. Chromosome 20 slice of the 2022-04-03 VCF release (`##fileDate=2022-04-03`, `##source=ClinVar`, `##reference=GRCh38`).

**Version**: 2022-04-03 monthly release, GRCh38. Committed alongside the gnomAD fixture below in the same early commit (`57c79da2`/`4cd76369`, "add test data"/"added test data").

**Files**:
- `clinvar_20220403_chr20.vcf.bgz`: 24,554 chromosome-20 variant records in ClinVar's distributed VCF layout (`#CHROM POS ID REF ALT QUAL FILTER INFO`), BGZF-compressed

**Licence**: Public domain. NCBI/NLM: "Information that is created by or for the US government on this site is within the public domain. Public domain information on the National Library of Medicine (NLM) Web pages may be freely distributed and copied." (https://www.ncbi.nlm.nih.gov/home/about/policies/). ClinVar's own page asks for attribution when data is redistributed: "If you distribute or copy data from ClinVar, we ask that you provide attribution to ClinVar as a data source in publications and websites." (https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/).

**Attribution**: ClinVar, NCBI/NLM — https://www.ncbi.nlm.nih.gov/clinvar/. Citation: Landrum MJ et al., PMID 29165669.

**Modifications**:
- None to the content: whole chromosome-20 records from the 2022-04-03 release, unedited

## gnomAD

### Fixture: `hvantk/tests/testdata/raw/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.chr20.tsv.bgz`

**Source**: gnomAD (Genome Aggregation Database), Broad Institute — per-gene loss-of-function constraint metrics (pLI, LOEUF/`oe_lof`, missense/synonymous observed-vs-expected), v2.1.1, chromosome 20 subset. Shared fixture: referenced in place by `gnomad_metrics` (`hvantk/skills/gnomad_metrics/tests/`) and `hvantk/tests/test_plugin_conformance.py`.

**Version**: v2.1.1 `lof_metrics.by_gene` table. Committed in the same early commit as the ClinVar fixture above (`57c79da2`/`4cd76369`).

**Files**:
- `gnomad.v2.1.1.lof_metrics.by_gene.chr20.tsv.bgz`: 542 gene rows from the v2.1.1 by-gene constraint table (e.g. `TOP1`, `CHD6`, `TAF4`), BGZF-compressed

**Licence**: CC0. gnomAD's own policies/terms page: "The primary data from the gnomAD exomes and genomes are available free of restrictions under the Creative Commons Zero Public Domain Dedication." (https://gnomad.broadinstitute.org/terms). Corroborated by the AWS Open Data Registry's mirror of the same statement: "The summary data provided here are released for the benefit of the wider scientific community without restriction on use." (https://registry.opendata.aws/broad-gnomad/).

**Attribution**: gnomAD, Broad Institute — https://gnomad.broadinstitute.org. Citation: Karczewski KJ et al., "The mutational constraint spectrum quantified from variation in 141,456 humans", *Nature* 581, 434–443 (2020).

**Modifications**:
- None to the content: whole gene rows for chromosome 20 from the v2.1.1 release, unedited

## GTEx

### Fixture: `hvantk/skills/gtex_eqtl/tests/testdata/raw/gtex-eqtl/Liver.v11.eQTLs.signif_pairs.parquet`

**Source**: GTEx (Genotype-Tissue Expression) Consortium release v11, cis-eQTL significant variant-gene pairs, Liver tissue (the smallest v11 tissue file, ~28 MB full-size).

**Version**: v11. 800 rows sliced from the full Liver `signif_pairs` parquet via `pyarrow.parquet.read_table(...).slice(0, 800)` (deterministic — row-group order is preserved), keeping the parquet binary format.

**Files**:
- `Liver.v11.eQTLs.signif_pairs.parquet`: 800 significant cis-eQTL variant-gene pairs (`phenotype_id`, `variant_id`, `pval_nominal`, `slope`, …), ~36 KB

**Licence**: open-access, but not unconditional — two distinct GTEx rules both apply, and the dbGaP line below is an embargo rule, not the whole licence. (1) **Embargo**: NIH dbGaP's GTEx study page (phs000424) states: "All GTEx releases since Version 5 follow the NIH Genomic Data Sharing (GDS) Policy whereby there are no restrictions on use or publication after release." (https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424) — this says a post-release embargo does not apply, not that there is no licence at all. (2) **GTEx's own data licence** (https://gtexportal.org/home/license, checked 2026-10-13) requires redistributors to "make known in a clear and conspicuous manner that [they] do not reflect the most current/accurate data" and to acknowledge the source, e.g. "the GTEx Portal on MM/DD/YY". See `README.md` beside this fixture for the notice text to carry forward with any redistribution.

**Attribution**: GTEx Consortium — https://gtexportal.org. Citation: GTEx Consortium, "The GTEx Consortium atlas of genetic regulatory effects across human tissues", *Science* 369, 1318–1330 (2020). Acknowledge per GTEx's own license: "the GTEx Portal on [retrieval date]".

**Modifications**:
- None to the content: the first 800 rows of the full Liver v11 `signif_pairs` file, unedited

## GWAS Catalog

### Fixture: `hvantk/skills/gwas_catalog/tests/testdata/raw/gwas-catalog/gwas-catalog-sample.tsv`

**Source**: NHGRI-EBI GWAS Catalog, full-associations download (v1.0 schema, 34 columns), release `2026-04-27` (`e115`).

**Version**: `e115_r2026-04-27`. 59 rows sliced from the 556 MB full-associations TSV, selected to exercise multi-trait-per-variant rows, scientific-notation p-values, and both populated and empty `OR or BETA` values.

**Files**:
- `gwas-catalog-sample.tsv`: 59 SNP-trait association rows in the live v1.0 column layout

**Licence**: EBI Services Terms of Use, not CC0 — CC0 at this source covers only the separate "GWAS Catalog summary statistics" product, not the curated associations file this fixture is sliced from. EBI's FAQ: "Summary statistics: GWAS Catalog summary statistics are released under the Creative Commons CC0 1.0 Public Domain Dedication, unless otherwise specified… Other GWAS Catalog data: Other data in the GWAS Catalog can be used under the EBI Services Terms of Use." (https://www.ebi.ac.uk/gwas/docs/faq). The EBI Services Terms of Use themselves state: "EMBL-EBI itself places no additional restrictions on the use or redistribution of the data available via its Data Resources and Tools other than those provided by the original data owners, unless otherwise specified in these Terms of Use." (https://www.ebi.ac.uk/about/terms-of-use) — i.e. EBI imposes no restriction of its own on this compiled associations data; responsibility for any restriction rests with the original contributing studies.

**Attribution**: NHGRI-EBI GWAS Catalog — https://www.ebi.ac.uk/gwas/. Citation: Sollis E et al., "The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource", *Nucleic Acids Research* 51(D1), D977–D985 (2023).

**Modifications**:
- None to the content: whole rows sliced from the full release, unedited

## HGNC

### Fixture: `hvantk/skills/hgnc/tests/testdata/raw/hgnc/hgnc_test_sample.tsv`

**Source**: HGNC (HUGO Gene Nomenclature Committee) complete-set TSV — approved human gene symbols, IDs, and curated cross-references.

**Version**: no dated upstream release (HGNC publishes a single rolling "complete set"); retrieved before the fixture's introducing commit (`4478cbbd`). The header was updated in `2491ff1b` (#381) to drop `location_sortable`, which upstream itself removed (#355).

**Files**:
- `hgnc_test_sample.tsv`: 53 columns, matching the current live complete-set header, plus 6 gene rows — 5 real (`BRCA1`/`HGNC:1100`, `BRCA2`/`HGNC:1101`, `GSTT1`/`HGNC:4641`, `GSTT1-AS1`/`HGNC:51839`, `FAM8A5P`/`HGNC:16376`, with their real Ensembl, UniProt, and MANE Select cross-references) and 1 fabricated (`HGNC:99999`/`FAKE1`, not a real HGNC id or symbol)

**Licence**: CC0. HGNC's own licence page: "all data is released under the Creative Commons Public Domain (CC0) License. This means that any form of reuse of the content is permitted." (https://www.genenames.org/about/license/).

**Attribution**: HGNC — https://www.genenames.org (not mandatory per the CC0 licence, but requested). Citation: Seal RL et al., "Genenames.org: the HGNC resources in 2023", *Nucleic Acids Research* (Database issue).

**Modifications**:
- `location_sortable` column dropped (#381), tracking an upstream schema change (#355)
- 1 of the 6 rows (`HGNC:99999`/`FAKE1`) is fabricated, not a real HGNC record

## 1000 Genomes / IGSR (onek_genomes)

### Fixture: `hvantk/skills/onek_genomes/tests/testdata/raw/onek_genomes/`

**Source**: 1000 Genomes Project high-coverage NYGC/CCDG callset (`1000G_2504_high_coverage`, GRCh38), distributed via the IGSR (International Genome Sample Resource) EBI mirror.

**Version**: NYGC/CCDG high-coverage release, chromosome 22 phased VCF, plus the `20130502` IGSR sample panel. Committed in `2560684e` ("add Phase B snapshot fixtures"), described there as "generated from a real fixture build".

**Files**:
- `chr22.vcf.gz` (+ `.tbi`): 50 variants sliced from the live chr22 high-coverage VCF
- `igsr_samples.tsv`: the full 2,504-row IGSR sample panel (sample ID, population, super-population, sex), unedited except for trimming trailing whitespace per line (the live panel's header declares 6 fields where the data has 4)

**Licence**: no embargo / no restriction, but with required acknowledgements. IGSR's disclaimer page: "Data from the 1000 Genomes Project is now available without embargo, following the final publication from the project. Use of the data should be cited in the usual way…" (https://www.internationalgenome.org/IGSR_disclaimer). IGSR is part of EMBL-EBI and is additionally covered by the general "Terms of Use for EMBL-EBI Services". The high-coverage callset's own reuse README (`20200526_1000G_2504plus698_high_cov_data_reuse_README.txt`, in the same FTP directory as the VCFs) requires, for the 2504-genome callset this fixture is sliced from: "The following cell lines/DNA samples were obtained from the NIGMS Human Genetic Cell Repository at the Coriell Institute for Medical Research… These data were generated at the New York Genome Center with funds provided by NHGRI Grant 3UM1HG008901-03S1."

**Attribution**: 1000 Genomes Project / IGSR, EMBL-EBI — https://www.internationalgenome.org. Citation: The 1000 Genomes Project Consortium, "A global reference for human genetic variation", *Nature* 526, 68–74 (2015); Byrska-Bishop M et al., "High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort", *Cell* 185, 3426–3440 (2022). Acknowledge NIGMS (Coriell Institute cell lines/DNA samples) and NHGRI Grant 3UM1HG008901-03S1 (New York Genome Center) per the callset's reuse README above.

**Modifications**:
- `chr22.vcf.gz`: whole records, unedited
- `igsr_samples.tsv`: whitespace-trimmed per line only (no field values changed)

## UCSC Cell Browser

### Fixture: `hvantk/skills/ucsc_cellbrowser/tests/testdata/raw/ucsc-cellbrowser/` (`default` dataset)

**Source**: **Corrected (review round)** — this fixture was previously attributed to a UCSC Cell Browser accession, `heart-dev-asp2019`, that does not exist (`https://cells.ucsc.edu/heart-dev-asp2019/dataset.json` 404s, and no dataset of that name appears in the live Cell Browser catalog). The actual source is Mendeley Data: "Developmental heart - filtered and unfiltered count matrices and meta tables" (Michaela Asp, 2021), the data deposit for Asp M et al., "A Spatiotemporal Organ-Wide Gene Expression and Cell Atlas of the Developing Human Heart", *Cell* 179(7), 1647–1660 (2019), PMID 31835037 — the paper citation itself was always correct, only the hosting accession was wrong.

**Version**: v2 (DOI `10.17632/mbvhhf8m62.2`, published 2021-02-27), as already cited in `docs_site/guide/data-sources.md`. Cell-type-level pseudo-bulk summary, computed via `hvantk expression summarize` and committed in `c48a8b3c` ("add sliced fixture from asp_2019 cell-type summary").

**Files**:
- `expression_matrix.tsv`: 250 genes x 12 cell-type pseudo-cells (e.g. Atrial cardiomyocytes, Capillary endothelium), ~22 KB
- `metadata.tsv`: per-cell-type `n_cells` counts

**Licence**: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/). Verified directly via the Mendeley Data API (`https://data.mendeley.com/public-api/datasets/mbvhhf8m62`): `"data_licence"` is `{"short_name": "CC BY 4.0", "full_name": "Creative Commons Attribution 4.0 International", "url": "http://creativecommons.org/licenses/by/4.0", "description": "You can share, copy and modify this dataset so long as you give appropriate credit, provide a link to the CC BY license, and indicate if changes were made, but you may not do so in a way that suggests the rights holder has endorsed you or your use of the dataset."}`.

**Attribution**: Michaela Asp, Mendeley Data — https://data.mendeley.com/datasets/mbvhhf8m62/2. Citation: Asp M et al., "A Spatiotemporal Organ-Wide Gene Expression and Cell Atlas of the Developing Human Heart", *Cell* 179(7), 1647–1660 (2019).

**Modifications** (CC BY 4.0 requires indicating changes):
- Values are per-cell-type pseudo-bulk **means** computed by `hvantk expression summarize` from the underlying single-cell data, not the original per-cell values
- Reduced from the full gene/cell-type set to 250 genes x 12 cell-type columns
- Three Fibroblast-like cell subtypes in the source were merged into the single `Fibroblast-like` column here
- Cell-type column labels may be abbreviated/truncated forms of the source's original cluster names
- Not independently re-verified byte-for-byte against the raw Mendeley count matrix in this pass (that would require downloading and reprocessing the full single-cell dataset); the source identity, DOI, version and licence above are independently confirmed

`adult-ctx` and `dev-ctx` fixtures are synthetic (see their SKILL.md section); no entry applies.

## MSigDB

### Fixture: `hvantk/skills/msigdb/tests/testdata/raw/msigdb/c2.cp-sample.gmt`

**Source**: Broad Institute MSigDB, collection C2 (curated gene sets), sub-collection CP (Canonical Pathways), human gene-symbol variant, release v2026.1.

**Version**: v2026.1. 11 gene sets sampled by line index from the live 4,115-row C2:CP source (#426 removed 9 sets originally in this fixture — 4 BioCarta and 5 legacy-KEGG sets — carved out of MSigDB's general licence grant; see Licence below).

**Files**:
- `c2.cp-sample.gmt`: 11 gene sets in GMT format (`SA_B_CELL_RECEPTOR_COMPLEXES`, `SA_CASPASE_CASCADE`, `PID_TCR_PATHWAY`, `PID_P53_DOWNSTREAM_PATHWAY`, `PID_HES_HEY_PATHWAY`, `REACTOME_MITOTIC_G1_PHASE_AND_G1_S_TRANSITION`, `REACTOME_DNA_REPAIR`, `REACTOME_CELL_CYCLE`, `REACTOME_POST_TRANSLATIONAL_PROTEIN_MODIFICATION`, `WP_APOPTOSIS`, `WP_TP53_NETWORK`), ranging 19-1,497 genes per set

**Licence**: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/) for the sets kept here. MSigDB's licence terms page: "The contents of MSigDB are protected by copyright… subject to the terms and conditions of the Creative Commons Attribution 4.0 International License. Some MSigDB gene sets are distributed with further restrictions." (https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp). That page names exactly four carve-out categories — KEGG MEDICUS (`KEGG_MEDICUS_` prefix, CC BY-SA 4.0), other/legacy KEGG (`KEGG_` prefix without `MEDICUS`, "provided for use in the MSigDB collection under license to the Broad Institute… qualified permission to include in this release"), BioCarta (`BIOCARTA_` prefix, "protected by copyright… provided under license to the Broad Institute… with qualified permission to include in this release"), and AAAS/STKE (`ST_` prefix, removed from MSigDB entirely since v7.0) — and none of the four sub-collections kept in this fixture (`SA`, `PID`, `REACTOME`, `WP`) is named among them, so they fall under the general CC BY 4.0 grant. The 9 removed sets (4 `BIOCARTA_*`, 5 legacy `KEGG_*`) fell under the legacy-KEGG/BioCarta carve-outs and are no longer part of this fixture.

**Attribution**: Broad Institute MSigDB — https://www.gsea-msigdb.org/gsea/msigdb/. Citation: Subramanian A et al., "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles", *PNAS* 102(43), 15545–15550 (2005); Liberzon A et al., "The Molecular Signatures Database (MSigDB) hallmark gene set collection", *Cell Systems* 1(6), 417–425 (2015).

**Modifications**:
- None to the content of the 11 kept rows: whole real gene-set rows, unedited
- 9 rows (4 BioCarta, 5 legacy KEGG) removed for licence reasons (#426); no replacement rows added
