"""Optional SuSiE-RSS + coloc.susie fine-mapping confirmation for a GWAS×eQTL locus.

The single-variant ABF in :mod:`gwas_coloc` is fast but can over-call when a
strong GWAS meets a weak eQTL. This module fine-maps both traits with SuSiE-RSS
(using a 1000G reference-LD matrix for a configurable super-population, default
EUR via ``--superpop``) and runs ``coloc.susie`` to test for a **shared credible
set** — separating genuine colocalization from single-variant artifacts.

Pure-Python and hvantk-native (issue #193): the fine-mapping math lives in
:mod:`hvantk.algorithms.qtlcascade.susie` (NumPy ports of ``susieR::susie_rss``
and ``coloc::coloc.susie``), and the 1000G reference LD is built from genotypes
streamed via **pysam** remote-tabix (a core dependency) — no external ``R``,
``bcftools`` or ``curl``. The LD reference can also be a local VCF (``ld_vcf``)
for offline / reproducible runs.

This layer is OPTIONAL and gracefully degrading: it needs the 1000G reference
(network, or a local ``ld_vcf``); :func:`run_finemap` reports ``available=False``
when the reference cannot be reached so callers fall back to the ABF-only verdict.

Caveat: reference (not in-sample) LD is the documented SuSiE-RSS limitation —
for weak/underpowered eQTLs the credible sets may be unreliable. The Python
SuSiE is not bit-identical to ``susieR`` but reproduces the R credible-set counts
and ``coloc.susie`` PP4 on the project's positive/negative controls (AF→MYOZ1
CONFIRMED; CHD 17q21/NSF REFUTED). The ``estimate_s_rss`` LD-mismatch diagnostic
(R-only, never used in the verdict) is dropped in the pure-Python port.
"""

from __future__ import annotations

import logging
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from hvantk.algorithms.qtlcascade.constants import (
    KG_PHASED_VCF_URL,
    KG_PANEL_URL,
    DEFAULT_FINEMAP_SUPERPOP,
    DEFAULT_SUSIE_L,
)
from hvantk.algorithms.qtlcascade.susie import susie_rss, coloc_susie

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Availability
# ---------------------------------------------------------------------------


def finemap_available() -> tuple[bool, list[str]]:
    """Return ``(ok, missing)``.

    The fine-mapping math is pure-Python (NumPy); the only hard requirement is
    ``pysam`` (a core hvantk dependency) to read the 1000G reference VCF. There
    is no longer any external ``R`` / ``bcftools`` / ``curl`` requirement — the
    actual LD reference (network or a local VCF) is checked at run time, with
    graceful degradation in :func:`run_finemap`.
    """
    missing: list[str] = []
    try:
        import pysam  # noqa: F401
    except Exception as exc:  # pragma: no cover - pysam is a core dependency
        missing.append(f"pysam ({exc})")
    return (len(missing) == 0), missing


# ---------------------------------------------------------------------------
# 1000G reference LD (configurable super-population, default EUR)
# ---------------------------------------------------------------------------


def _cache_dir(ld_cache_dir: Optional[str]) -> Path:
    d = Path(ld_cache_dir) if ld_cache_dir else Path(
        os.environ.get("HVANTK_LD_CACHE", Path.home() / ".cache" / "hvantk" / "1kg"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def _download_file(url: str, dest: Path, timeout: int = 300, retries: int = 3) -> None:
    """Download ``url`` to ``dest`` with retries (replaces the old ``curl`` calls).

    Writes to a ``*.part`` temp and atomically renames on success, validating the
    byte count against ``Content-Length`` when present — so a truncated stream
    never leaves a corrupt cached file that later runs would reuse forever.
    """
    import requests

    part = dest.with_name(dest.name + ".part")
    last: Optional[Exception] = None
    for attempt in range(retries):
        try:
            with requests.get(url, stream=True, timeout=timeout) as r:
                r.raise_for_status()
                expected = r.headers.get("Content-Length")
                with open(part, "wb") as fh:
                    for chunk in r.iter_content(chunk_size=1 << 16):
                        if chunk:
                            fh.write(chunk)
            size = part.stat().st_size
            if size == 0:
                raise IOError("empty download")
            if expected is not None and size != int(expected):
                raise IOError(f"truncated download: {size} != {expected} bytes")
            os.replace(part, dest)            # atomic; never leaves a partial dest
            return
        except Exception as exc:  # network flakiness; retry with backoff
            last = exc
            logger.warning("download attempt %d/%d failed for %s: %s",
                           attempt + 1, retries, url, exc)
            part.unlink(missing_ok=True)
            if attempt + 1 < retries:
                time.sleep(min(2 ** attempt, 8))
    raise RuntimeError(f"failed to download {url}: {last}")


def _unrelated_samples(cache: Path, superpop: str) -> list[str]:
    """Download the 1000G 3202 panel once; derive unrelated founders of ``superpop``."""
    out = cache / f"{superpop.lower()}_unrelated.txt"
    if out.exists():
        return [s.strip() for s in out.read_text().splitlines() if s.strip()]
    panel = cache / "g1k_3202.panel"
    if not panel.exists():
        _download_file(KG_PANEL_URL, panel, timeout=120)
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
    if not samples:
        raise ValueError(
            f"No unrelated founders found for super-population '{superpop}' in the "
            "1000G panel. Check --superpop (one of EUR, AFR, EAS, SAS, AMR)."
        )
    out.write_text("\n".join(samples) + "\n")
    logger.info("derived %d unrelated %s founders", len(samples), superpop)
    return samples


def _kg_index(chrom: str, cache: Path) -> str:
    """Download the per-chrom 1000G tabix index locally (remote index is flaky)."""
    idx = cache / f"kg_chr{chrom}.tbi"
    if not idx.exists():
        _download_file(KG_PHASED_VCF_URL.format(chrom=chrom) + ".tbi", idx, timeout=300)
    return str(idx)


def _resolve_contig(vf, chrom: str) -> str:
    """Resolve a contig name against a VCF's naming (with/without ``chr``)."""
    names = set(vf.header.contigs)
    for cand in (chrom, f"chr{chrom}", chrom[3:] if chrom.startswith("chr") else None):
        if cand and cand in names:
            return cand
    return f"chr{chrom}"


def _panel_dosages(
    chrom: str, start: int, end: int, cache: Path, superpop: str,
    ld_vcf: Optional[str] = None,
) -> dict[tuple[int, str, str], np.ndarray]:
    """Region ALT-dosage matrix over the reference samples (biallelic SNPs).

    Streams genotypes with **pysam** — either from a local ``ld_vcf`` (all its
    samples) or from the remote 1000G phased VCF subset to the unrelated
    ``superpop`` founders. ALT dosage per sample = number of ALT alleles
    (``NaN`` when missing). Retries the remote fetch (EBI streams can drop).
    """
    import pysam

    if ld_vcf:
        sources = [(ld_vcf, None, None)]            # local VCF: use all samples
    else:
        samples = _unrelated_samples(cache, superpop)
        idx = _kg_index(chrom, cache)
        url = KG_PHASED_VCF_URL.format(chrom=chrom)
        sources = [(url, idx, samples)]

    src, index_filename, keep = sources[0]
    last: Optional[Exception] = None
    for attempt in range(8):
        out: dict[tuple[int, str, str], np.ndarray] = {}
        try:
            kwargs = {"index_filename": index_filename} if index_filename else {}
            with pysam.VariantFile(src, **kwargs) as vf:
                if keep is not None:
                    header_samples = set(vf.header.samples)  # build once, not per ID
                    present = [s for s in keep if s in header_samples]
                    if not present:
                        raise RuntimeError(
                            f"none of the {len(keep)} {superpop} samples are in the "
                            "reference VCF header")
                    vf.subset_samples(present)
                contig = _resolve_contig(vf, chrom)
                for rec in vf.fetch(contig, max(0, start - 1), end):
                    if len(rec.ref) != 1 or not rec.alts or len(rec.alts) != 1 \
                            or len(rec.alts[0]) != 1:
                        continue  # biallelic SNPs only (mirrors bcftools -m2 -M2 -v snps)
                    dos = []
                    for s in rec.samples.values():
                        gt = s.get("GT")
                        if gt is None or any(a is None for a in gt):
                            dos.append(np.nan)
                        else:
                            dos.append(float(sum(1 for a in gt if a == 1)))
                    out[(rec.pos, rec.ref, rec.alts[0])] = np.array(dos, dtype=float)
            return out
        except Exception as exc:
            last = exc
            logger.warning("1000G region fetch attempt %d/8 failed: %s", attempt + 1, exc)
            time.sleep(min(2 ** attempt, 8))  # backoff: don't burn all retries on a blip
    raise RuntimeError(f"1000G region fetch failed after retries: {last}")


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


@dataclass(frozen=True)
class _Harmonized:
    snps: list[str]
    z_gwas: np.ndarray
    z_eqtl: np.ndarray
    R: np.ndarray


def _build_inputs(gwas, eqtl_recs, chrom, start, end, cache: Path, superpop: str,
                  ld_vcf: Optional[str] = None) -> _Harmonized:
    """Harmonize GWAS ∩ eQTL ∩ 1000G panel → z-scores + signed LD matrix.

    Effects oriented to the GWAS ALT allele; LD computed as ALT-dosage
    correlation (over the chosen super-population) so its sign matches the betas.
    This ALT-sign handling (flip the eQTL beta and invert the dosage ``2 - d`` on
    a ref/alt swap) is preserved exactly from the R-bridge implementation.
    """
    kg = _panel_dosages(chrom, start, end, cache, superpop, ld_vcf=ld_vcf)
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
        if np.isnan(dv).all():
            continue  # all genotypes missing — would poison the LD correlation
        if np.isnan(dv).any():
            dv = np.where(np.isnan(dv), np.nanmean(dv), dv)
        if not np.isfinite(dv).all() or np.std(dv) == 0:
            continue
        rows.append((f"{chrom}:{pos}:{ref}:{alt}", pos, bg, sg, be, se))
        dosages.append(dv)
    if len(rows) < 2:
        return _Harmonized(snps=[r[0] for r in rows], z_gwas=np.zeros(0),
                           z_eqtl=np.zeros(0), R=np.zeros((0, 0)))
    if len(rows) > 5000:
        logger.warning(
            "%d variants in the LD region; the %dx%d correlation matrix is large "
            "and SuSiE may be slow/memory-heavy. Consider a smaller --window-kb.",
            len(rows), len(rows), len(rows),
        )
    order = np.argsort([r[1] for r in rows])
    rows = [rows[i] for i in order]
    R = np.corrcoef(np.array([dosages[i] for i in order]))
    z_gwas = np.array([r[2] / r[3] for r in rows])
    z_eqtl = np.array([r[4] / r[5] for r in rows])
    return _Harmonized(snps=[r[0] for r in rows], z_gwas=z_gwas, z_eqtl=z_eqtl, R=R)


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
    ld_vcf: Optional[str] = None,
    work_dir: Optional[str] = None,
) -> FineMapResult:
    """Fine-map a single GWAS×eQTL locus and run coloc.susie (pure-Python).

    ``eqtl_recs`` is the per-gene record list from
    :func:`gwas_coloc.fetch_eqtl_region` for the gene of interest. ``ld_vcf``,
    when given, is a local VCF used as the LD reference instead of the remote
    1000G panel (offline / reproducible runs).
    """
    ok, missing = finemap_available()
    if not ok:
        return FineMapResult(available=False,
                             note="fine-map skipped; missing: " + ", ".join(missing))
    cache = _cache_dir(ld_cache_dir)
    try:
        h = _build_inputs(gwas, eqtl_recs, chrom, start, end, cache, superpop,
                          ld_vcf=ld_vcf)
    except Exception as exc:
        # LD reference unobtainable is a *known* graceful-degrade path -> available=False
        # so the verdict falls back to SUGGESTIVE (ABF only). An unexpected SuSiE-kernel
        # crash below is treated differently (available=True, pp4=None -> INCONCLUSIVE).
        logger.warning("fine-map LD reference unavailable, continuing without it: %s", exc)
        return FineMapResult(available=False, note=f"LD reference unavailable: {exc}")

    n = len(h.snps)
    if n < 2:
        return FineMapResult(available=True, n_variants=n,
                             note="too few overlapping variants for fine-mapping")

    if work_dir:  # optional debug dump of the exact SuSiE inputs
        try:
            wd = Path(work_dir)
            wd.mkdir(parents=True, exist_ok=True)
            np.savetxt(wd / "ld.tsv", h.R, fmt="%.6f", delimiter="\t")
            with open(wd / "merged.tsv", "w") as fh:
                fh.write("snp\tz_gwas\tz_eqtl\n")
                for snp, zg, ze in zip(h.snps, h.z_gwas, h.z_eqtl):
                    fh.write(f"{snp}\t{zg}\t{ze}\n")
        except Exception as exc:  # pragma: no cover - debug aid only
            logger.debug("could not write fine-map debug files: %s", exc)

    try:
        sg = susie_rss(h.z_gwas, h.R, gwas_N, L=DEFAULT_SUSIE_L)
        se = susie_rss(h.z_eqtl, h.R, eqtl_N, L=DEFAULT_SUSIE_L)
        cs_g, cs_e = len(sg.cs), len(se.cs)
        pp4 = coloc_susie(sg, se)  # 0.0 when either trait has no credible set
        return FineMapResult(
            available=True, n_variants=n,
            cs_gwas=cs_g, cs_eqtl=cs_e, coloc_susie_pp4=pp4,
            ld_s_gwas=None, ld_s_eqtl=None,  # estimate_s_rss diagnostic dropped (R-only)
            note="ok",
        )
    except Exception as exc:  # fine-mapping is optional: degrade, don't abort the run
        logger.warning("fine-mapping failed, continuing without it: %s", exc)
        return FineMapResult(available=True, n_variants=n,
                             note=f"fine-map error: {exc}")
