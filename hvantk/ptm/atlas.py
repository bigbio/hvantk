"""Phase-2 PTM atlas assembly API.

Thin facade over :func:`hvantk.ptm.pipeline.ptm_build_pipeline` that reproduces
notebook A's ``ptm_sites_combined.tsv.bgz`` output. Defaults match notebook A:
UniProt + PeptideAtlas sources, CPTAC disabled, flanking window of 7 codons.

This module does NOT reimplement the download, coordinate mapping, or
concat logic - it delegates entirely to the shipped build pipeline so the
resulting TSV is byte-identical to the existing workflow.

Example
-------
    >>> from hvantk.ptm.atlas import PTMAtlasConfig, build_atlas
    >>> config = PTMAtlasConfig(
    ...     output_dir="data/ptm/",
    ...     output_ht="data/ptm/ptm_sites.ht",
    ... )
    >>> result = build_atlas(config)
    >>> print(result.combined_tsv, result.n_sites)
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from typing import List, Optional, Sequence

from hvantk.ptm.pipeline import PTMBuildConfig, PTMBuildResult, ptm_build_pipeline

logger = logging.getLogger(__name__)

# Default atlas sources: CPTAC is OFF by default because its download path
# requires the optional ``cptac`` Python package and per-cancer-type fetches.
# Notebook A uses all three, but the minimal reproducible atlas is UniProt +
# PeptideAtlas.
DEFAULT_ATLAS_SOURCES: Sequence[str] = ("uniprot", "peptideatlas")

# Known source identifiers (used for validation).
_KNOWN_SOURCES = frozenset({"uniprot", "peptideatlas", "cptac"})


@dataclass
class PTMAtlasConfig:
    """Configuration for a Phase-2 PTM atlas build.

    Attributes
    ----------
    output_dir : str
        Directory for intermediate + combined TSV output (required).
    output_ht : str
        Path to the final PTM sites Hail Table (required).
    sources : Sequence[str]
        Sources to include; any subset of {"uniprot", "peptideatlas", "cptac"}.
        Defaults to ("uniprot", "peptideatlas") - matches the minimal notebook-A
        atlas that produces ``ptm_sites_combined.tsv.bgz``.
    uniprot_tsv : Optional[str]
        Pre-downloaded UniProt PTM TSV. If None and "uniprot" is selected,
        the pipeline will download it via the REST API.
    peptideatlas_tsv : Optional[str]
        Pre-downloaded PeptideAtlas phospho intermediate TSV.
    cptac_tsv : Optional[str]
        Pre-downloaded CPTAC phospho combined TSV.
    gtf_path : Optional[str]
        Pre-downloaded Ensembl GTF (auto-downloaded if None).
    flanking_codons : int
        Flanking-codon window for proximity intervals. Phase-2 default is 7
        (matches notebook A); the shipped pipeline default is 5.
    overwrite : bool
        Force re-run of downstream steps even if outputs exist.
    """

    output_dir: str = ""
    output_ht: str = ""
    sources: Sequence[str] = field(default_factory=lambda: list(DEFAULT_ATLAS_SOURCES))
    uniprot_tsv: Optional[str] = None
    peptideatlas_tsv: Optional[str] = None
    cptac_tsv: Optional[str] = None
    gtf_path: Optional[str] = None
    flanking_codons: int = 7
    overwrite: bool = False

    def validate(self) -> List[str]:
        """Return a list of configuration errors (empty if valid)."""
        errors: List[str] = []
        if not self.output_dir:
            errors.append("output_dir is required")
        if not self.output_ht:
            errors.append("output_ht is required")
        if not self.sources:
            errors.append("at least one source is required")
        unknown = [s for s in self.sources if s not in _KNOWN_SOURCES]
        if unknown:
            errors.append(
                f"unknown source(s): {unknown}; must be subset of {sorted(_KNOWN_SOURCES)}"
            )
        # Per-source file existence checks mirror PTMBuildConfig.validate().
        if self.uniprot_tsv and not os.path.exists(self.uniprot_tsv):
            errors.append(f"uniprot_tsv not found: {self.uniprot_tsv}")
        if self.peptideatlas_tsv and not os.path.exists(self.peptideatlas_tsv):
            errors.append(f"peptideatlas_tsv not found: {self.peptideatlas_tsv}")
        if self.cptac_tsv and not os.path.exists(self.cptac_tsv):
            errors.append(f"cptac_tsv not found: {self.cptac_tsv}")
        if self.gtf_path and not os.path.exists(self.gtf_path):
            errors.append(f"gtf_path not found: {self.gtf_path}")
        return errors


@dataclass
class PTMAtlasResult:
    """Result of a Phase-2 PTM atlas build.

    Attributes
    ----------
    output_ht : str
        Path to the final PTM sites Hail Table.
    combined_tsv : str
        Path to the combined BGZ TSV produced by the pipeline:
        ``ptm_sites_combined.tsv.bgz`` when CPTAC is disabled, or
        ``ptm_sites_all_combined.tsv.bgz`` when CPTAC is included.
    n_sites : int
        Number of PTM sites successfully mapped (summed across sources).
    sources_used : list[str]
        Sources included in this build (subset of DEFAULT_ATLAS_SOURCES + cptac).
    """

    output_ht: str = ""
    combined_tsv: str = ""
    n_sites: int = 0
    sources_used: List[str] = field(default_factory=list)


def build_atlas(config: PTMAtlasConfig) -> PTMAtlasResult:
    """Build a Phase-2 PTM atlas by delegating to ``ptm_build_pipeline``.

    Reproduces notebook A's ``ptm_sites_combined.tsv.bgz`` exactly - no
    cross-source deduplication is performed (matches the existing pipeline).

    Parameters
    ----------
    config : PTMAtlasConfig
        Build configuration.

    Returns
    -------
    PTMAtlasResult
        Output paths, mapped-site count, and list of sources used.

    Raises
    ------
    ValueError
        If config.validate() returns any errors.
    """
    errors = config.validate()
    if errors:
        raise ValueError(f"Invalid PTMAtlasConfig: {'; '.join(errors)}")

    sources = [s.lower() for s in config.sources]
    include_uniprot = "uniprot" in sources
    include_peptideatlas = "peptideatlas" in sources
    include_cptac = "cptac" in sources

    # Translate the Phase-2 config into the legacy PTMBuildConfig. Unselected
    # sources are passed as None, which disables the corresponding pipeline
    # step. UniProt is always the primary source in the legacy pipeline; if
    # the caller opted it out, we fall back to the REST download anyway
    # (pipeline.ptm_build_pipeline requires ptm_tsv to run).
    if not include_uniprot:
        # The shipped pipeline treats UniProt as the mandatory primary source
        # (step 3 always runs). Warn rather than silently including it.
        logger.warning(
            "uniprot not listed in sources; pipeline still runs UniProt "
            "as the mandatory primary source. Add 'uniprot' to sources "
            "to suppress this warning."
        )

    build_cfg = PTMBuildConfig(
        output_dir=config.output_dir,
        output_ht=config.output_ht,
        gtf_path=config.gtf_path,
        ptm_tsv=config.uniprot_tsv,
        peptideatlas_tsv=config.peptideatlas_tsv if include_peptideatlas else None,
        cptac_tsv=config.cptac_tsv if include_cptac else None,
        flanking_codons=config.flanking_codons,
        overwrite=config.overwrite,
    )

    logger.info(
        "Building Phase-2 PTM atlas (sources=%s, flanking_codons=%d)",
        sources,
        config.flanking_codons,
    )
    build_result: PTMBuildResult = ptm_build_pipeline(build_cfg)

    # The pipeline emits:
    #   - ptm_sites_mapped.tsv.bgz when only UniProt is used
    #   - ptm_sites_combined.tsv.bgz when UniProt + PeptideAtlas
    #   - ptm_sites_all_combined.tsv.bgz when CPTAC is also included
    # Report the final combined TSV path for reproducibility with notebook A.
    combined_tsv = build_result.mapped_tsv_path
    if include_cptac:
        # ptm_build_pipeline points mapped_tsv_path at the combined (no-CPTAC)
        # TSV; the CPTAC-augmented file has a distinct name.
        expected = os.path.join(config.output_dir, "ptm_sites_all_combined.tsv.bgz")
        if os.path.exists(expected):
            combined_tsv = expected

    return PTMAtlasResult(
        output_ht=build_result.output_ht,
        combined_tsv=combined_tsv,
        n_sites=build_result.n_mapped,
        sources_used=list(sources),
    )


__all__ = [
    "DEFAULT_ATLAS_SOURCES",
    "PTMAtlasConfig",
    "PTMAtlasResult",
    "build_atlas",
]
