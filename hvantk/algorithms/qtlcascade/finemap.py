"""Optional SuSiE-RSS + coloc.susie fine-mapping confirmation for a GWAS×eQTL locus.

The single-variant ABF in :mod:`gwas_coloc` is fast but can over-call when a
strong GWAS meets a weak eQTL. This module fine-maps both traits with SuSiE-RSS
(using a 1000G EUR reference-LD matrix) and runs ``coloc.susie`` to test for a
**shared credible set** — separating genuine colocalization from single-variant
artifacts.

This layer is OPTIONAL: it needs external tools (``R`` with ``susieR``+``coloc``,
``bcftools``) and network access to the 1000G reference. :func:`finemap_available`
reports what's missing so callers can skip it gracefully.

Caveat: reference (not in-sample) LD is the documented SuSiE-RSS limitation —
for weak/underpowered eQTLs the credible sets may be unreliable; the LD-mismatch
diagnostics (``estimate_s_rss``) are reported alongside the result.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from hvantk.algorithms.qtlcascade.constants import (
    KG_PHASED_VCF_URL,
    KG_PANEL_URL,
    DEFAULT_FINEMAP_SUPERPOP,
)

logger = logging.getLogger(__name__)

_R_SCRIPT = Path(__file__).resolve().parent / "resources" / "susie_coloc.R"
_R_LIB_PREFIX = (
    'ul<-Sys.getenv("R_LIBS_USER"); if(nzchar(ul)) .libPaths(c(ul,.libPaths())); '
)


# ---------------------------------------------------------------------------
# Availability
# ---------------------------------------------------------------------------


def finemap_available() -> tuple[bool, list[str]]:
    """Return ``(ok, missing)`` — whether R+packages and bcftools are present."""
    missing: list[str] = []
    if shutil.which("bcftools") is None:
        missing.append("bcftools")
    if shutil.which("curl") is None:
        missing.append("curl")
    if shutil.which("Rscript") is None:
        missing.append("Rscript (R)")
        return False, missing
    try:
        r = subprocess.run(
            ["Rscript", "-e", _R_LIB_PREFIX
             + 'cat(all(sapply(c("susieR","coloc","jsonlite"),requireNamespace,quietly=TRUE)))'],
            capture_output=True, text=True, timeout=120,
        )
        if "TRUE" not in r.stdout:
            missing.append("R packages: susieR, coloc, jsonlite")
    except Exception as exc:  # pragma: no cover
        missing.append(f"R check failed: {exc}")
    return (len(missing) == 0), missing


# ---------------------------------------------------------------------------
# 1000G EUR reference LD
# ---------------------------------------------------------------------------


def _cache_dir(ld_cache_dir: Optional[str]) -> Path:
    d = Path(ld_cache_dir) if ld_cache_dir else Path(
        os.environ.get("HVANTK_LD_CACHE", Path.home() / ".cache" / "hvantk" / "1kg"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def _eur_unrelated_samples(cache: Path, superpop: str) -> Path:
    """Download the 1000G 3202 panel once; derive unrelated founders of ``superpop``."""
    out = cache / f"{superpop.lower()}_unrelated.txt"
    if out.exists():
        return out
    panel = cache / "g1k_3202.panel"
    if not panel.exists():
        subprocess.run(["curl", "-sSL", "--retry", "3", "--max-time", "120",
                        "-o", str(panel), KG_PANEL_URL], check=True)
    samples = []
    with open(panel) as fh:
        header = fh.readline().split()
        idx = {name: i for i, name in enumerate(header)}
        sid, fat, mot, sup = (idx["SampleID"], idx["FatherID"],
                              idx["MotherID"], idx["Superpopulation"])
        for line in fh:
            f = line.split()
            if len(f) > sup and f[sup] == superpop and f[fat] == "0" and f[mot] == "0":
                samples.append(f[sid])
    out.write_text("\n".join(samples) + "\n")
    logger.info("derived %d unrelated %s founders", len(samples), superpop)
    return out


def _kg_index(chrom: str, cache: Path) -> str:
    """Download the per-chrom 1000G tabix index locally (remote index is flaky)."""
    idx = cache / f"kg_chr{chrom}.tbi"
    if not idx.exists():
        url = KG_PHASED_VCF_URL.format(chrom=chrom) + ".tbi"
        subprocess.run(["curl", "-sSL", "--retry", "3", "--max-time", "300",
                        "-o", str(idx), url], check=True)
    return str(idx)


def _eur_dosages(chrom: str, start: int, end: int, cache: Path, superpop: str
                 ) -> dict[tuple[int, str, str], np.ndarray]:
    """Region ALT-dosage matrix over unrelated `superpop` samples (biallelic SNPs).

    Downloads the region subset to a local VCF first (EBI streams flaky BGZF
    mid-pipe), retrying until bcftools exits cleanly, then queries locally.
    """
    samples = _eur_unrelated_samples(cache, superpop)
    idx = _kg_index(chrom, cache)
    url = KG_PHASED_VCF_URL.format(chrom=chrom) + f"##idx##{idx}"
    reg = cache / f"region_chr{chrom}_{start}_{end}_{superpop.lower()}.vcf.gz"
    if not reg.exists():
        ok = False
        for attempt in range(8):
            p = subprocess.run(
                ["bcftools", "view", "-r", f"chr{chrom}:{start}-{end}",
                 "-S", str(samples), "--force-samples", "-m2", "-M2", "-v", "snps",
                 url, "-Oz", "-o", str(reg)],
                capture_output=True, text=True,
            )
            if p.returncode == 0 and reg.exists() and reg.stat().st_size > 1000:
                ok = True
                break
            logger.warning("1000G region fetch attempt %d failed (rc=%d)",
                           attempt + 1, p.returncode)
            reg.unlink(missing_ok=True)
        if not ok:
            raise RuntimeError("1000G region fetch failed after retries")
    out: dict[tuple[int, str, str], np.ndarray] = {}
    q = subprocess.run(["bcftools", "query", "-f", "%POS\t%REF\t%ALT[\t%GT]\n", str(reg)],
                       capture_output=True, text=True)
    if q.returncode != 0:
        raise RuntimeError(f"bcftools query failed: {q.stderr.strip()[-200:]}")
    for line in q.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 4:
            continue
        dos = []
        for gt in f[3:]:
            gt = gt.replace("|", "/")
            dos.append(np.nan if "." in gt else sum(1 for a in gt.split("/") if a == "1"))
        out[(int(f[0]), f[1], f[2])] = np.array(dos, dtype=float)
    return out


# ---------------------------------------------------------------------------
# Fine-map driver
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FineMapResult:
    available: bool
    n_variants: int = 0
    cs_gwas: Optional[int] = None
    cs_eqtl: Optional[int] = None
    coloc_susie_pp4: Optional[float] = None
    ld_s_gwas: Optional[float] = None
    ld_s_eqtl: Optional[float] = None
    note: str = ""


def _build_inputs(gwas, eqtl_recs, chrom, start, end, gwas_N, eqtl_N,
                  cache: Path, work: Path, superpop: str) -> int:
    """Harmonize GWAS ∩ eQTL ∩ 1000G-EUR → write merged.tsv, ld.tsv, meta.json.

    Effects oriented to the GWAS ALT allele; LD computed as ALT-dosage
    correlation so its sign matches the betas. Returns the variant count.
    """
    kg = _eur_dosages(chrom, start, end, cache, superpop)
    eqtl_by_key = {key: (b, s) for key, b, s, _p in eqtl_recs}
    rows, dosages = [], []
    for (pos, ref, alt), (bg, sg, _pg) in gwas.items():
        be = se = None
        if (pos, ref, alt) in eqtl_by_key:
            be, se = eqtl_by_key[(pos, ref, alt)]
        elif (pos, alt, ref) in eqtl_by_key:
            be, se = eqtl_by_key[(pos, alt, ref)]
            be = -be
        if be is None:
            continue
        dv = None
        if (pos, ref, alt) in kg:
            dv = kg[(pos, ref, alt)]
        elif (pos, alt, ref) in kg:
            dv = 2.0 - kg[(pos, alt, ref)]
        if dv is None:
            continue
        if np.isnan(dv).any():
            dv = np.where(np.isnan(dv), np.nanmean(dv), dv)
        if np.nanstd(dv) == 0:
            continue
        rows.append((f"{chrom}:{pos}:{ref}:{alt}", pos, bg, sg, be, se))
        dosages.append(dv)
    if len(rows) < 2:
        return len(rows)
    order = np.argsort([r[1] for r in rows])
    rows = [rows[i] for i in order]
    R = np.corrcoef(np.array([dosages[i] for i in order]))
    with open(work / "merged.tsv", "w") as fh:
        fh.write("snp\tpos\tbeta_gwas\tse_gwas\tbeta_eqtl\tse_eqtl\n")
        for snp, pos, bg, sg, be, se in rows:
            fh.write(f"{snp}\t{pos}\t{bg}\t{sg}\t{be}\t{se}\n")
    np.savetxt(work / "ld.tsv", R, fmt="%.6f", delimiter="\t")
    with open(work / "meta.json", "w") as fh:
        json.dump({"gwas_N": gwas_N, "eqtl_N": eqtl_N}, fh)
    return len(rows)


def run_finemap(
    *,
    gwas: dict[tuple[int, str, str], tuple[float, float, float]],
    eqtl_recs: list[tuple[tuple[int, str, str], float, float, float]],
    chrom: str,
    start: int,
    end: int,
    gwas_N: int,
    eqtl_N: int,
    ld_cache_dir: Optional[str] = None,
    superpop: str = DEFAULT_FINEMAP_SUPERPOP,
    work_dir: Optional[str] = None,
) -> FineMapResult:
    """Fine-map a single GWAS×eQTL locus and run coloc.susie.

    ``eqtl_recs`` is the per-gene record list from
    :func:`gwas_coloc.fetch_eqtl_region` for the gene of interest.
    """
    ok, missing = finemap_available()
    if not ok:
        return FineMapResult(available=False,
                             note="fine-map skipped; missing: " + ", ".join(missing))
    cache = _cache_dir(ld_cache_dir)
    import tempfile
    work = Path(work_dir) if work_dir else Path(tempfile.mkdtemp(prefix="hvantk_finemap_"))
    work.mkdir(parents=True, exist_ok=True)
    try:
        n = _build_inputs(gwas, eqtl_recs, chrom, start, end, gwas_N, eqtl_N,
                          cache, work, superpop)
        if n < 2:
            return FineMapResult(available=True, n_variants=n,
                                 note="too few overlapping variants for fine-mapping")
        p = subprocess.run(["Rscript", str(_R_SCRIPT), str(work)],
                           capture_output=True, text=True, timeout=1800)
        vals: dict[str, str] = {}
        for line in p.stdout.splitlines():
            parts = line.split()
            if len(parts) == 2:
                vals[parts[0]] = parts[1]
        if "DONE" not in p.stdout:
            return FineMapResult(available=True, n_variants=n,
                                 note=f"R fine-map did not complete: {p.stderr.strip()[-200:]}")

        def _f(k):
            v = vals.get(k)
            try:
                return float(v)
            except (TypeError, ValueError):
                return None

        def _i(k):
            v = vals.get(k)
            try:
                return int(v)
            except (TypeError, ValueError):
                return None

        return FineMapResult(
            available=True, n_variants=n,
            cs_gwas=_i("CS_GWAS"), cs_eqtl=_i("CS_EQTL"),
            coloc_susie_pp4=_f("COLOC_SUSIE_PP4"),
            ld_s_gwas=_f("LD_S_GWAS"), ld_s_eqtl=_f("LD_S_EQTL"),
            note="ok",
        )
    finally:
        if work_dir is None:
            shutil.rmtree(work, ignore_errors=True)
