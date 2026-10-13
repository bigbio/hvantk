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

**Version**: the AlphaGenome SDK and hosted model as of that run; the exact SDK version was not recorded (the run's manifest records the backend, scorer set, interval and organism, but no SDK version). Source shard `scores.shard-00000-of-00001.parquet`, 16,642,316 rows, sha256 `261cf87273f1b36d20db640305f414b2e9bdc55824cddff9e4cc4bcc4c6745a1`.

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

**Licence**: no restriction. The AWS Open Data Registry entry for gnomAD, maintained alongside the Broad Institute's own distribution, states: "The summary data provided here are released for the benefit of the wider scientific community without restriction on use." (https://registry.opendata.aws/broad-gnomad/, which also links the Broad's own `gnomad.broadinstitute.org/terms` for the full terms).

**Attribution**: gnomAD, Broad Institute — https://gnomad.broadinstitute.org. Citation: Karczewski KJ et al., "The mutational constraint spectrum quantified from variation in 141,456 humans", *Nature* 581, 434–443 (2020).

**Modifications**:
- None to the content: whole gene rows for chromosome 20 from the v2.1.1 release, unedited

## GTEx

### Fixture: `hvantk/skills/gtex_eqtl/tests/testdata/raw/gtex-eqtl/Liver.v11.eQTLs.signif_pairs.parquet`

**Source**: GTEx (Genotype-Tissue Expression) Consortium release v11, cis-eQTL significant variant-gene pairs, Liver tissue (the smallest v11 tissue file, ~28 MB full-size).

**Version**: v11. 800 rows sliced from the full Liver `signif_pairs` parquet via `pyarrow.parquet.read_table(...).slice(0, 800)` (deterministic — row-group order is preserved), keeping the parquet binary format.

**Files**:
- `Liver.v11.eQTLs.signif_pairs.parquet`: 800 significant cis-eQTL variant-gene pairs (`phenotype_id`, `variant_id`, `pval_nominal`, `slope`, …), ~36 KB

**Licence**: open, no restriction on use or publication, for this open-access summary-statistics tier (as distinct from dbGaP-controlled individual-level genotypes, which this fixture is not). NIH dbGaP's own GTEx study page (phs000424) states: "All GTEx releases since Version 5 follow the NIH Genomic Data Sharing (GDS) Policy whereby there are no restrictions on use or publication after release." (https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424).

**Attribution**: GTEx Consortium — https://gtexportal.org. Citation: GTEx Consortium, "The GTEx Consortium atlas of genetic regulatory effects across human tissues", *Science* 369, 1318–1330 (2020).

**Modifications**:
- None to the content: the first 800 rows of the full Liver v11 `signif_pairs` file, unedited

## GWAS Catalog

### Fixture: `hvantk/skills/gwas_catalog/tests/testdata/raw/gwas-catalog/gwas-catalog-sample.tsv`

**Source**: NHGRI-EBI GWAS Catalog, full-associations download (v1.0 schema, 34 columns), release `2026-04-27` (`e115`).

**Version**: `e115_r2026-04-27`. 59 rows sliced from the 556 MB full-associations TSV, selected to exercise multi-trait-per-variant rows, scientific-notation p-values, and both populated and empty `OR or BETA` values.

**Files**:
- `gwas-catalog-sample.tsv`: 59 SNP-trait association rows in the live v1.0 column layout

**Licence**: CC0 1.0. EBI's own page: "The NHGRI-EBI GWAS Catalog and all its contents are available under the general Terms of Use for EMBL-EBI Services. Summary statistics are made available under CC0 unless otherwise stated." (https://www.ebi.ac.uk/gwas/docs/about).

**Attribution**: NHGRI-EBI GWAS Catalog — https://www.ebi.ac.uk/gwas/. Citation: Sollis E et al., "The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource", *Nucleic Acids Research* 51(D1), D977–D985 (2023).

**Modifications**:
- None to the content: whole rows sliced from the full release, unedited

## HGNC

### Fixture: `hvantk/skills/hgnc/tests/testdata/raw/hgnc/hgnc_test_sample.tsv`

**Source**: HGNC (HUGO Gene Nomenclature Committee) complete-set TSV — approved human gene symbols, IDs, and curated cross-references.

**Version**: no dated upstream release (HGNC publishes a single rolling "complete set"); retrieved before the fixture's introducing commit (`4478cbbd`).

**Files**:
- `hgnc_test_sample.tsv`: all 52 columns of the live complete-set header, plus 6 real gene rows (including `BRCA1`/`HGNC:1100` and `BRCA2`/`HGNC:1101`, with their real Ensembl, UniProt, and MANE Select cross-references)

**Licence**: CC0. HGNC's own licence page: "all data is released under the Creative Commons Public Domain (CC0) License. This means that any form of reuse of the content is permitted." (https://www.genenames.org/about/license/).

**Attribution**: HGNC — https://www.genenames.org (not mandatory per the CC0 licence, but requested). Citation: Seal RL et al., "Genenames.org: the HGNC resources in 2023", *Nucleic Acids Research* (Database issue).

**Modifications**:
- None to the content: whole real rows, unedited

## UniProt (uniprot_ptm)

### Fixture: `hvantk/skills/uniprot_ptm/tests/testdata/raw/uniprot-ptm/ptm_sites.tsv`

**Source**: UniProtKB MOD_RES (post-translational modification) features, mapped to genomic coordinates by `hvantk.algorithms.ptm.pipeline.map_ptm_sites`.

**Version**: no archival version at the query level (UniProt's live search endpoint); the version date is the date of download. Introduced in the plugin's original commit.

**Files**:
- `ptm_sites.tsv`: 6 rows — 5 real phosphosite records on real UniProt accessions (`P04637`/TP53 x2, `P51587`/BRCA2, `Q9Y6K9`/NEMO, `P03886`/MT-ND1), plus 1 fabricated edge-case row (`X9XXX9`/`UNPLACED` on contig `GL000009.2`) exercising unplaced-contig handling

**Licence**: CC BY 4.0. UniProt's licence page states the Creative Commons Attribution 4.0 International License applies to all copyrightable parts of the UniProt databases (https://www.uniprot.org/help/license).

**Attribution**: UniProt Consortium — https://www.uniprot.org. Citation: UniProt Consortium, "UniProt: the Universal Protein Knowledgebase in 2025", *Nucleic Acids Research* 53(D1), D609–D617 (2025).

**Modifications**:
- The 5 real rows are unedited UniProt feature values; the 6th (`X9XXX9`/`UNPLACED`) row is fabricated to exercise an edge case and contains no real UniProt accession

## 1000 Genomes / IGSR (onek_genomes)

### Fixture: `hvantk/skills/onek_genomes/tests/testdata/raw/onek_genomes/`

**Source**: 1000 Genomes Project high-coverage NYGC/CCDG callset (`1000G_2504_high_coverage`, GRCh38), distributed via the IGSR (International Genome Sample Resource) EBI mirror.

**Version**: NYGC/CCDG high-coverage release, chromosome 22 phased VCF, plus the `20130502` IGSR sample panel. Committed in `2560684e` ("add Phase B snapshot fixtures"), described there as "generated from a real fixture build".

**Files**:
- `chr22.vcf.gz` (+ `.tbi`): 50 variants sliced from the live chr22 high-coverage VCF
- `igsr_samples.tsv`: the full 2,504-row IGSR sample panel (sample ID, population, super-population, sex), unedited except for trimming trailing whitespace per line (the live panel's header declares 6 fields where the data has 4)

**Licence**: no embargo / no restriction. IGSR's disclaimer page: "Data from the 1000 Genomes Project is now available without embargo, following the final publication from the project. Use of the data should be cited in the usual way…" (https://www.internationalgenome.org/IGSR_disclaimer). IGSR is part of EMBL-EBI and is additionally covered by the general "Terms of Use for EMBL-EBI Services".

**Attribution**: 1000 Genomes Project / IGSR, EMBL-EBI — https://www.internationalgenome.org. Citation: The 1000 Genomes Project Consortium, "A global reference for human genetic variation", *Nature* 526, 68–74 (2015); Byrska-Bishop M et al., "High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort", *Cell* 185, 3426–3440 (2022).

**Modifications**:
- `chr22.vcf.gz`: whole records, unedited
- `igsr_samples.tsv`: whitespace-trimmed per line only (no field values changed)

## UCSC Cell Browser

### Fixture: `hvantk/skills/ucsc_cellbrowser/tests/testdata/raw/ucsc-cellbrowser/` (`default` dataset)

**Source**: UCSC Cell Browser, "Asp 2019 Developing Heart" collection (`heart-dev-asp2019`), derived from Asp et al. 2019 developing human heart scRNA-seq (PMID 31835037).

**Version**: cell-type-level pseudo-bulk summary, committed in `c48a8b3c` ("add sliced fixture from asp_2019 cell-type summary").

**Files**:
- `expression_matrix.tsv`: 250 genes x 12 cell-type pseudo-cells (e.g. Atrial cardiomyocytes, Capillary endothelium), ~22 KB
- `metadata.tsv`: per-cell-type `n_cells` counts

**Licence**: no dataset-specific restriction found on the UCSC Cell Browser dataset page; treated under the general UCSC data policy, which "Projects" (including UCSC Cell Browser) inherit: "No license is needed for the data files and database tables used by the Genome Browser… the raw table data and binary files used to create the graphics by the browser is freely available for both public and commercial use." (https://genome.ucsc.edu/license/). That page names specific restricted sources (HGMD, LOVD, OMIM, Decipher, Genomenon, Genehancer, COSMIC, GISAID-derived data); none apply here.

**Attribution**: UCSC Cell Browser — https://cells.ucsc.edu/?ds=heart-dev-asp2019. Citation: Asp M et al., "A Spatiotemporal Organ-Wide Gene Expression and Cell Atlas of the Developing Human Heart", *Cell* 179(7), 1647–1660 (2019).

**Modifications**:
- Reduced to 250 genes x 12 cell types from the full collection; values are the collection's own pseudo-bulk summary statistics, unedited

### Fixture: `hvantk/skills/ucsc_cellbrowser/tests/testdata/raw/ucsc-cellbrowser-adult-ctx/` (`adult-ctx` dataset)

**Source**: UCSC Cell Browser, "Adult Cortex Meta-Atlas" collection (`adult-ctx-meta-atlas`).

**Version**: class-level pseudo-bulk summary, committed in `02583310` ("tier2 D1: onboard adult_ctx_meta_atlas_class_summary").

**Files**:
- `expression_matrix.tsv`: 250 genes x 19 cell-class pseudo-cells (e.g. Astrocyte, IT, PVALB), ~22 KB
- `metadata.tsv`: per-class `n_cells` counts

**Licence**: as above (UCSC Cell Browser general data policy, https://genome.ucsc.edu/license/).

**Attribution**: UCSC Cell Browser — https://cells.ucsc.edu/?ds=adult-ctx-meta-atlas.

**Modifications**:
- Reduced to 250 genes x 19 classes from the full collection; values unedited

### Fixture: `hvantk/skills/ucsc_cellbrowser/tests/testdata/raw/ucsc-cellbrowser-dev-ctx/` (`dev-ctx` dataset)

**Source**: UCSC Cell Browser, "Developing Cortex Meta-Atlas" collection (`dev-ctx-meta-atlas`).

**Version**: Type_v2-level pseudo-bulk summary, committed in `ff0805a2` ("tier2 D2: onboard dev_ctx_meta_atlas_type_v2_summary").

**Files**:
- `expression_matrix.tsv`: 250 genes x 15 Type_v2 pseudo-cells (e.g. Astrocyte, Excitatory Neuron, Microglia, RG), ~22 KB
- `metadata.tsv`: per-type `n_cells` counts

**Licence**: as above (UCSC Cell Browser general data policy, https://genome.ucsc.edu/license/).

**Attribution**: UCSC Cell Browser — https://cells.ucsc.edu/?ds=dev-ctx-meta-atlas.

**Modifications**:
- Reduced to 250 genes x 15 types from the full collection; values unedited

## MSigDB

### Fixture: `hvantk/skills/msigdb/tests/testdata/raw/msigdb/c2.cp-sample.gmt`

**Source**: Broad Institute MSigDB, collection C2 (curated gene sets), sub-collection CP (Canonical Pathways), human gene-symbol variant, release v2026.1.

**Version**: v2026.1. 11 gene sets sampled by line index from the live 4,115-row C2:CP source (#426 removed 9 sets originally in this fixture — 4 BioCarta and 5 legacy-KEGG sets — carved out of MSigDB's general licence grant; see Licence below).

**Files**:
- `c2.cp-sample.gmt`: 11 gene sets in GMT format (`SA_B_CELL_RECEPTOR_COMPLEXES`, `SA_CASPASE_CASCADE`, `PID_TCR_PATHWAY`, `PID_P53_DOWNSTREAM_PATHWAY`, `PID_HES_HEY_PATHWAY`, `REACTOME_MITOTIC_G1_PHASE_AND_G1_S_TRANSITION`, `REACTOME_DNA_REPAIR`, `REACTOME_CELL_CYCLE`, `REACTOME_POST_TRANSLATIONAL_PROTEIN_MODIFICATION`, `WP_APOPTOSIS`, `WP_TP53_NETWORK`), ranging 19-1,497 genes per set

**Licence**: CC BY 4.0 for the sets kept here. MSigDB's licence terms page: "The contents of MSigDB are protected by copyright… subject to the terms and conditions of the Creative Commons Attribution 4.0 International License. Some MSigDB gene sets are distributed with further restrictions." (https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp). That page names exactly four carve-out categories — KEGG MEDICUS (`KEGG_MEDICUS_` prefix, CC BY-SA 4.0), other/legacy KEGG (`KEGG_` prefix without `MEDICUS`, "provided for use in the MSigDB collection under license to the Broad Institute… qualified permission to include in this release"), BioCarta (`BIOCARTA_` prefix, "protected by copyright… provided under license to the Broad Institute… with qualified permission to include in this release"), and AAAS/STKE (`ST_` prefix, removed from MSigDB entirely since v7.0) — and none of the four sub-collections kept in this fixture (`SA`, `PID`, `REACTOME`, `WP`) is named among them, so they fall under the general CC BY 4.0 grant. The 9 removed sets (4 `BIOCARTA_*`, 5 legacy `KEGG_*`) fell under the legacy-KEGG/BioCarta carve-outs and are no longer part of this fixture.

**Attribution**: Broad Institute MSigDB — https://www.gsea-msigdb.org/gsea/msigdb/. Citation: Subramanian A et al., "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles", *PNAS* 102(43), 15545–15550 (2005); Liberzon A et al., "The Molecular Signatures Database (MSigDB) hallmark gene set collection", *Cell Systems* 1(6), 417–425 (2015).

**Modifications**:
- None to the content of the 11 kept rows: whole real gene-set rows, unedited
- 9 rows (4 BioCarta, 5 legacy KEGG) removed for licence reasons (#426); no replacement rows added
