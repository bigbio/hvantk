"""Download stage for `onek-genomes:samples` — fetches the 1KG samples panel."""
from __future__ import annotations

import logging
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

# Canonical 1000 Genomes samples panel for the 2504-sample cohort
# (sample / pop / super_pop / gender, tab-separated). The high-coverage
# NYGC callset (1000G_2504_high_coverage) shares these samples, so this
# panel is the correct metadata source.
_IGSR_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)
_TIMEOUT_S = 60


def download_igsr_samples(*, raw_dir, **params) -> None:
    """Fetch the samples panel to ``raw_dir/igsr_samples.tsv``.

    The samples panel is small (~55 KB) and updates infrequently. A new
    download overwrites any prior copy in ``raw_dir``.
    """
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    target = raw_dir / "igsr_samples.tsv"

    logger.info("Downloading 1KG samples panel from %s", _IGSR_URL)
    resp = requests.get(_IGSR_URL, timeout=_TIMEOUT_S)
    resp.raise_for_status()
    # The upstream panel has trailing tabs on the header line — strip them so
    # Hail's import_table sees consistent field counts across header and data.
    cleaned = "\n".join(line.rstrip() for line in resp.text.splitlines()) + "\n"
    target.write_text(cleaned, encoding="utf-8")
    logger.info("Wrote %d bytes to %s", len(cleaned), target)
