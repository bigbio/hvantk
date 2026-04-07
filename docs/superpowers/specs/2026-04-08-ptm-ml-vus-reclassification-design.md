# PTM-ML: Phosphorylation-Aware VUS Reclassification

**Date**: 2026-04-08
**Status**: Draft
**Scope**: Phosphorylation only (extensible to other PTM types)

## Problem Statement

Millions of Variants of Uncertain Significance (VUS) exist in ClinVar. For the subset of VUS at or near experimentally validated phosphorylation sites, can multi-source phosphoproteomics evidence (UniProt + CPTAC + PeptideAtlas) provide discriminative signal to reclassify them toward pathogenic or benign?

This is not a general-purpose pathogenicity predictor. It is a specialized evidence layer that evaluates whether phosphorylation-site context adds signal on top of existing predictors (REVEL, CADD) for variants within ±7 amino acid residues of phosphorylation sites.

## Framing

- **Task**: Binary classification — (Pathogenic + Likely Pathogenic) vs (Benign + Likely Benign)
- **Instances**: Missense variants within ±7 amino acid residues of experimentally validated phosphorylation sites
- **Target population**: ClinVar VUS near phospho sites (never used in training)
- **Output**: Continuous pathogenicity score (0–1) mapped to ACMG-aligned evidence tiers:
  - Strong evidence toward pathogenic
  - Moderate evidence toward pathogenic
  - Insufficient evidence
  - Moderate evidence toward benign
  - Strong evidence toward benign
- **Tier thresholds**: Derived from model calibration on the P/B training set

## Training Set Construction

### Variant Selection Pipeline

1. Build phospho-site Hail Table from UniProt + CPTAC + PeptideAtlas (existing PTM pipeline, `flanking_codons=7`)
2. Load ClinVar variant table (existing builder: `create_clinvar_tb`)
3. Annotate ClinVar variants against phospho sites (existing: `annotate_variants_with_ptm`)
4. Filter to: `ptm_distance <= 7` AND missense AND clinical significance in {P, LP, B, LB}

### Labels

- Positive (1): Pathogenic + Likely Pathogenic
- Negative (0): Benign + Likely Benign
- VUS: Prediction target only, never in training

### Feasibility Check (Implementation Phase)

- Count P/LP and B/LB variants within ±7 residues of phospho sites. If < 200 per class, approach needs rethinking.
- Assess class imbalance ratio.

### Items Flagged for Revision

- **ClinVar missense filtering**: Need to verify how to identify missense variants from ClinVar VCF without VEP. May need ClinVar's own molecular consequence field or VEP preprocessing.
- **Amino acid change features**: Deferred. Ref/Alt AA and physicochemical change properties depend on ClinVar annotation completeness verification. If available, these are strong candidate features (e.g., S->A eliminates phosphorylation, S->D is a phosphomimetic).

## Feature Set

### Core Features

| Feature | Type | Source | Description |
|---|---|---|---|
| `ptm_distance` | Continuous (0–7) | PTM annotation | Amino acid distance to nearest phospho site |
| `is_ptm_site` | Binary | PTM annotation | Variant directly on the phosphorylated residue |
| `amino_acid` | Categorical (S/T/Y) | PTM annotation | The phosphorylated residue type |
| `source_db` | Categorical / multi-hot | PTM annotation | UniProt / CPTAC / PeptideAtlas. Multi-source presence is itself a signal. |
| `evidence_type` | Categorical | PTM annotation | Curated vs MS-derived |
| `n_observations` | Continuous | PTM annotation | Number of MS observations (evidence strength). Zero for UniProt curated sites — may need binary `has_ms_evidence` flag. |
| `gnomad_af` | Continuous | Variant builder | Population allele frequency |
| `ccr_percentile` | Continuous | Variant builder | Constrained coding region score |
| `revel_score` | Continuous (0–1) | dbNSFP | REVEL pathogenicity score |
| `cadd_phred` | Continuous | dbNSFP | CADD PHRED-scaled score |

### Design Notes

- `source_db` encoded as multi-hot: a site observed in UniProt AND CPTAC AND PeptideAtlas gets [1,1,1].
- `n_observations` is zero for UniProt curated sites (no MS count). Consider separate handling or a binary `has_ms_evidence` flag alongside.
- All features are variant-level. No gene-level features (pLI, LOEUF, etc.) — these provide zero within-gene discrimination and risk gene identity leakage.
- REVEL and CADD may have missing values for some variants. Strategy: impute with median, or use binary `has_revel`/`has_cadd` flag + imputed value. Investigate during implementation.

### Deferred Features

- Ref/Alt amino acid and physicochemical change properties — pending ClinVar annotation verification.

## ML Pipeline

### Algorithm

Random Forest (sklearn `RandomForestClassifier`). Consistent with existing hvantk infrastructure (ancestry classifier pattern).

### Pipeline Steps

1. **Data assembly** — Hail-side joins: ClinVar variants + PTM annotation + gnomAD AF + CCR + dbNSFP (REVEL, CADD). Export to pandas DataFrame.
2. **Preprocessing** — One-hot encode categoricals (`amino_acid`, `evidence_type`). Multi-hot encode `source_db`. Handle missing values for dbNSFP scores.
3. **Training** — Stratified 5-fold cross-validation. Primary metrics: AUC and AUPRC.
4. **Feature importance** — Gini importance + permutation importance per fold.
5. **Ablation** — Repeat CV without PTM features (only gnomAD AF + CCR + REVEL + CADD). Compare AUC. This is the core scientific experiment: does PTM context add signal on top of state-of-the-art scores?
6. **Calibration** — Platt scaling or isotonic regression to produce calibrated probabilities (0–1).
7. **Tier mapping** — Map calibrated scores to ACMG evidence tiers using threshold optimization (Youden's J on training folds, or fixed precision/recall targets).
8. **VUS prediction** — Apply calibrated model to VUS set. Report tier assignments and cross-reference with REVEL/AlphaMissense.

### Output Artifacts

- Trained model (joblib serialization)
- CV performance metrics (AUC, AUPRC, per-fold)
- Feature importance rankings (plot + table)
- Ablation comparison (with vs without PTM features, plot)
- VUS reclassification table (variant, score, tier, REVEL/AlphaMissense comparison)

## Evaluation Strategy

### Layer 1: Internal Validation
Stratified 5-fold cross-validation on P/LP vs B/LB training set. Metrics: AUC, AUPRC, sensitivity, specificity, calibration plots.

### Layer 2: Feature Ablation (Core Scientific Claim)
- **Model A** (baseline): REVEL + CADD + gnomAD AF + CCR — no PTM features
- **Model B** (PTM-aware): Model A + all PTM features
- Compare AUC/AUPRC. Statistical significance via paired test across CV folds.
- If Model B > Model A: PTM context adds independent signal. This is the publishable finding.

### Layer 3: Feature Importance
- Gini importance + permutation importance from Random Forest
- Ranking of PTM features relative to REVEL/CADD
- Per-fold stability of importance rankings

### Layer 4: VUS Application
- Apply Model B to ClinVar VUS within ±7 AA of phospho sites
- Report: number of VUS reclassified per evidence tier at different confidence thresholds
- Cross-reference: for reclassified VUS, what do REVEL and AlphaMissense predict? Agreements and disagreements.
- Disagreements (PTM-ML says pathogenic, REVEL says benign, or vice versa) are the most scientifically interesting cases.

## What hvantk Supports vs What Needs Building

### Already Exists

| Component | Location |
|---|---|
| UniProt + CPTAC + PeptideAtlas phospho download/mapping | `hvantk/ptm/pipeline.py`, `hvantk/datasets/` |
| PTM Hail Table builder with flanking intervals | `hvantk/tables/table_builders.py` |
| Variant-PTM annotation (distance, type, proximity) | `hvantk/ptm/annotate.py` |
| ClinVar table builder | `hvantk/tables/table_builders.py` |
| gnomAD, CCR table builders | `hvantk/tables/table_builders.py` |
| sklearn Random Forest pattern | `hvantk/ancestry/classify.py` |
| ROC evaluation pattern | `hvantk/psroc/roc.py` |

### Needs Building

| Component | Description |
|---|---|
| Training set assembler | Join ClinVar + PTM annotation + gnomAD + CCR + dbNSFP into a single feature matrix. Hail-side joins, export to pandas DataFrame. |
| dbNSFP score extraction | One-time extraction of REVEL/CADD for training variants. Targeted lookup, not a full builder. |
| RF classifier module | Train, evaluate, ablation, calibration. Reuse patterns from ancestry but PTM-ML specific. |
| Tier mapper | Calibrated score to ACMG evidence tier assignment with configurable thresholds. |
| VUS predictor | Apply model to VUS, generate reclassification table with cross-references. |
| Feature importance reporting | Gini + permutation importance plots and tables. |
| CLI commands | `hvantk ptm-ml train`, `hvantk ptm-ml predict`, `hvantk ptm-ml evaluate` (or similar). |

### Flagged for Revision

| Item | Concern |
|---|---|
| ClinVar missense filtering | Verify how to identify missense variants from ClinVar data without VEP dependency |
| Amino acid change features | Deferred — depends on ClinVar annotation completeness |

## Scope Boundaries

### In Scope
- Phosphorylation sites only (S/T/Y)
- ClinVar variants as training/prediction source
- Multi-source phospho evidence (UniProt, CPTAC, PeptideAtlas)
- Random Forest classifier with sklearn
- Feature ablation study (with/without PTM features)
- VUS reclassification with ACMG-aligned tiers

### Out of Scope (Future Work)
- Other PTM types (acetylation, methylation, ubiquitination, glycosylation, sumoylation) — architecture supports extension
- Deep learning approaches (CNN, transformers, protein language models)
- 3D structural features (AlphaFold distance-based proximity)
- General-purpose pathogenicity prediction (competing with REVEL/AlphaMissense)
- VEP integration
- Full dbNSFP builder in hvantk
