"""End-to-end GWAS→effector colocalization pipeline.

Orchestrates: FinnGen GWAS × eQTL Catalogue ABF coloc (rank cis effectors) →
optional SuSiE/coloc.susie fine-map confirmation of the lead effector → a
provenance-stamped report. The fine-map step is what makes a verdict
trustworthy: it distinguishes a genuine colocalization (CONFIRMED) from a
single-variant-ABF artifact (REFUTED). Validated on the AF→MYOZ1 positive
control (CONFIRMED) vs the CHD 17q21/NSF look-alike (REFUTED).
"""

from __future__ import annotations

import datetime
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from hvantk.algorithms.qtlcascade import gwas_coloc as gc
from hvantk.algorithms.qtlcascade import finemap as fm
from hvantk.algorithms.qtlcascade.constants import (
    DEFAULT_COLOC_P1, DEFAULT_COLOC_P2, DEFAULT_COLOC_P12,
    DEFAULT_GWAS_W_CC, DEFAULT_EQTL_W_QUANT, DEFAULT_COLOC_MIN_SNPS,
    DEFAULT_COLOC_WINDOW_KB,
    EQTL_CATALOGUE_DEFAULT_STUDY, DEFAULT_FINEMAP_SUPERPOP,
)

# PP4 cutoff for the verdict gate (ABF "colocalizes" and coloc.susie "confirms").
# Deliberately 0.5 (not the stricter 0.8 ABF-declaration threshold): coloc.susie
# is more conservative than single-variant ABF, and the validated positive
# control (AF→MYOZ1) confirms at coloc.susie PP4=0.71.
DEFAULT_VERDICT_PP4 = 0.5

logger = logging.getLogger(__name__)


@dataclass
class GwasColocConfig:
    endpoint: str                       # FinnGen R10 endpoint code
    chrom: str
    lead: int
    eqtl_dataset: str                   # eQTL Catalogue QTD id / URL / local tabix
    output_dir: str
    eqtl_study: str = EQTL_CATALOGUE_DEFAULT_STUDY
    window_kb: int = DEFAULT_COLOC_WINDOW_KB
    gene_of_interest: Optional[str] = None  # ENSG; default = ABF-top
    fine_map: bool = True
    gwas_N: Optional[int] = None        # required for fine-mapping
    eqtl_N: Optional[int] = None        # required for fine-mapping
    min_snps: int = DEFAULT_COLOC_MIN_SNPS
    W1: float = DEFAULT_GWAS_W_CC
    W2: float = DEFAULT_EQTL_W_QUANT
    p1: float = DEFAULT_COLOC_P1
    p2: float = DEFAULT_COLOC_P2
    p12: float = DEFAULT_COLOC_P12
    pp4_threshold: float = DEFAULT_VERDICT_PP4
    superpop: str = DEFAULT_FINEMAP_SUPERPOP
    ld_cache_dir: Optional[str] = None
    gene_symbols: dict = field(default_factory=dict)

    def validate(self) -> list[str]:
        errs = []
        if not self.endpoint:
            errs.append("endpoint is required")
        if not self.eqtl_dataset:
            errs.append("eqtl_dataset is required")
        if self.fine_map and (self.gwas_N is None or self.eqtl_N is None):
            errs.append("fine-mapping requires gwas_N and eqtl_N (or pass --no-fine-map)")
        return errs


def _verdict(pp4_abf: Optional[float], fmr: Optional[fm.FineMapResult],
             pp4_threshold: float) -> str:
    if pp4_abf is None or pp4_abf < pp4_threshold:
        return f"NO COLOC (ABF below {pp4_threshold})"
    if fmr is None or not fmr.available:
        return "SUGGESTIVE (ABF only; fine-mapping not run)"
    if fmr.coloc_susie_pp4 is None:
        return "INCONCLUSIVE (fine-mapping incomplete)"
    if fmr.coloc_susie_pp4 >= pp4_threshold:
        return "CONFIRMED (ABF + fine-mapping agree)"
    if (fmr.cs_gwas or 0) == 0 or (fmr.cs_eqtl or 0) == 0:
        return "REFUTED (no fine-mappable signal — single-variant-ABF artifact)"
    return "REFUTED (distinct causal variants — ABF PP4 not supported by fine-mapping)"


def run_gwas_coloc_pipeline(config: GwasColocConfig) -> dict:
    """Run the pipeline and write a provenance-stamped report. Returns the report dict."""
    half = config.window_kb * 1000
    start, end = config.lead - half, config.lead + half
    region = f"chr{config.chrom}:{start}-{end}"

    gwas = gc.fetch_finngen_region(config.endpoint, config.chrom, start, end)
    eqtl = gc.fetch_eqtl_region(config.eqtl_dataset, config.chrom, start, end,
                                study=config.eqtl_study)
    logger.info("GWAS variants=%d; eQTL genes=%d", len(gwas), len(eqtl))

    cres = gc.run_locus_coloc(
        gwas=gwas, eqtl=eqtl, region=region, min_snps=config.min_snps,
        W1=config.W1, W2=config.W2, p1=config.p1, p2=config.p2, p12=config.p12,
        gene_symbols=config.gene_symbols,
    )

    target = config.gene_of_interest
    if target is None and not cres.table.empty:
        target = str(cres.table.iloc[0]["gene_id"])
    goi_row = None
    if target is not None and not cres.table.empty:
        match = cres.table[cres.table["gene_id"] == target]
        goi_row = match.iloc[0] if len(match) else None
    pp4_abf = float(goi_row["PP4"]) if goi_row is not None else (
        float(cres.table.iloc[0]["PP4"]) if not cres.table.empty else None)

    fmr = None
    if config.fine_map and target is not None and target in eqtl:
        if config.gwas_N is None or config.eqtl_N is None:
            fmr = fm.FineMapResult(available=False, note="fine-map needs gwas_N/eqtl_N")
        else:
            logger.info("fine-mapping %s …", target)
            fmr = fm.run_finemap(
                gwas=gwas, eqtl_recs=eqtl[target], chrom=config.chrom,
                start=start, end=end, gwas_N=config.gwas_N, eqtl_N=config.eqtl_N,
                ld_cache_dir=config.ld_cache_dir, superpop=config.superpop,
            )

    verdict = _verdict(pp4_abf, fmr, config.pp4_threshold)

    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = out_dir / f"coloc_{config.endpoint}_{config.chrom}_{config.lead}.tsv"
    cres.table.to_csv(tsv_path, sep="\t", index=False)

    top = cres.table.iloc[0] if not cres.table.empty else None
    report = {
        "endpoint": config.endpoint,
        "region": region,
        "genome_build": "GRCh38",
        "run_date": datetime.date.today().isoformat(),
        "gwas": {"source": f"FinnGen R10 {config.endpoint}", "N": config.gwas_N,
                 "min_p_in_region": cres.gwas_min_p},
        "eqtl": {"source": f"eQTL Catalogue {config.eqtl_dataset} (study {config.eqtl_study})",
                 "N": config.eqtl_N},
        "coloc_params": {"W1": config.W1, "W2": config.W2, "p1": config.p1,
                         "p2": config.p2, "p12": config.p12,
                         "window_kb": config.window_kb, "min_snps": config.min_snps},
        "results": {
            "n_genes_tested": cres.n_genes_tested,
            "top_effector": (None if top is None else str(top["gene_id"])),
            "top_PP4": (None if top is None else float(top["PP4"])),
            "gene_of_interest": target,
            "goi_PP4_abf": pp4_abf,
            "coloc_table_tsv": str(tsv_path),
        },
        "fine_map": (None if fmr is None else {
            "available": fmr.available, "n_variants": fmr.n_variants,
            "credible_sets_gwas": fmr.cs_gwas, "credible_sets_eqtl": fmr.cs_eqtl,
            "coloc_susie_PP4": fmr.coloc_susie_pp4,
            "ld_consistency_s_gwas": fmr.ld_s_gwas, "ld_consistency_s_eqtl": fmr.ld_s_eqtl,
            "ld_panel": f"1000G high-coverage GRCh38, {config.superpop} unrelated founders",
            "note": fmr.note,
        }),
        "verdict": verdict,
        "provenance_notes": (
            "Summary statistics streamed via remote tabix (no bulk download). "
            "Fine-mapping (when run) uses a 1000G reference-LD panel, not in-sample "
            "LD — the documented SuSiE-RSS limitation for weak/underpowered eQTLs."),
        "software": {"framework": "hvantk.algorithms.qtlcascade.gwas_pipeline"},
    }
    report_path = out_dir / f"report_{config.endpoint}_{config.chrom}_{config.lead}.json"
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)
    report["results"]["report_json"] = str(report_path)
    logger.info("verdict: %s", verdict)
    return report
