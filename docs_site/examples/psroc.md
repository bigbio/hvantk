# PSROC Examples

PSROC evaluates variant pathogenicity prediction scores (CADD, REVEL, MetaLR, etc.) against ClinVar labels, generating ROC curves and performance metrics.

See the [PSROC reference](../tools/psroc.md) for full CLI options and output format details.

## CLI Usage

```bash
hvantk psroc \
  --genes-file test_genes.txt \
  --clinvar-ht clinvar.ht \
  --dbnsfp-ht dbnsfp.ht \
  --scores "CADD_phred,REVEL_score,MetaLR_score" \
  --output-dir results/ \
  --max-missingness 0.3 \
  --threshold-method youden
```

## Python API

```python
from hvantk.psroc.pipeline import PSROCConfig, run_psroc_pipeline

config = PSROCConfig(
    genes=["BRCA1", "BRCA2", "TP53"],
    scores=["CADD_phred", "REVEL_score", "MetaLR_score"],
    max_missingness=0.3,
    threshold_method="youden",
)

run_psroc_pipeline(
    clinvar_ht_path="clinvar.ht",
    dbnsfp_ht_path="dbnsfp.ht",
    output_dir="results/",
    config=config,
)
```

## Pipeline Stages

1. **Load Tables** - Load ClinVar and dbNSFP Hail Tables
2. **Filter ClinVar** - Filter by genes/variants and review status
3. **Assign Labels** - Convert CLNSIG to binary labels (1=pathogenic, 0=benign)
4. **Annotate Scores** - Join prediction scores from dbNSFP
5. **Compute Missingness** - Calculate missing data rates per score
6. **Compute ROC** - Generate ROC curves and AUC metrics
7. **Generate Outputs** - Create plots, metrics JSON, and TSV exports

## Understanding Results

**ROC curves**: Higher AUC = better discriminator (1.0 = perfect, 0.5 = random). Optimal thresholds are marked using Youden's index by default.

**Missingness filtering**: Scores with >30% missing values are automatically excluded:

```
VEST4_score: 42.2% missing -> EXCLUDED
REVEL_score: 2.0% missing  -> INCLUDED
```

## Expected Outputs

- `psroc_*_roc_curves.png` - ROC curves for all included scores
- `psroc_*_auc_comparison.png` - AUC comparison bar chart
- `psroc_*_dashboard.png` - Combined dashboard view
- `psroc_*_metrics.json` - AUC values and optimal thresholds
- `psroc_*_annotated.tsv` - Annotated variants with scores and labels

## Troubleshooting

**No variants after filtering**: Check gene names match ClinVar (HGNC symbols), verify review status threshold, confirm input tables overlap.

**Scores excluded**: Adjust `--max-missingness` threshold. Some scores have high missingness by design.

**Plots not generating**: Ensure matplotlib backend is configured and output directory is writable. Plots are on by default; use `--no-plots` to skip.

## Runnable Scripts

See the [`examples/psroc/`](https://github.com/bigbio/hvantk/tree/main/examples/psroc/) directory for a complete end-to-end example with synthetic test data.
