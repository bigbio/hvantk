"""Workflow orchestration for PTM build pipeline.

Wraps the pure algorithm core (:func:`hvantk.algorithms.ptm.pipeline.ptm_build_pipeline_core`)
with the skill-driven download and Hail Table build steps. This is the layer
that knows about :mod:`hvantk.skills.uniprot_ptm`; the algorithm layer stays pure.
"""
from __future__ import annotations

import logging

from hvantk.algorithms.ptm.pipeline import (
    PTMBuildConfig,
    PTMBuildResult,
    ptm_build_pipeline_core,
)

logger = logging.getLogger(__name__)


def download_uniprot_ptm(output_dir: str, overwrite: bool = False) -> str:
    """Download the latest UniProt PTM TSV via :class:`UniProtPTMDataset`.

    Skill-driven; lives in tools/ because it depends on
    :mod:`hvantk.skills.uniprot_ptm`.

    Parameters
    ----------
    output_dir : str
        Directory to save the TSV file.
    overwrite : bool
        If True, re-download even if a file exists.

    Returns
    -------
    str
        Path to the downloaded TSV file.
    """
    from hvantk.skills.uniprot_ptm.shared.datasets import UniProtPTMDataset

    dataset = UniProtPTMDataset.from_latest()
    return dataset.download(output_dir, overwrite=overwrite)


def ptm_build_pipeline(config: PTMBuildConfig) -> PTMBuildResult:
    """End-to-end PTM build pipeline: download + core mapping + Hail Table build.

    Steps:
        1. Download UniProt PTM TSV (if ``config.ptm_tsv`` is None)
        2. Delegate to :func:`ptm_build_pipeline_core` for GTF download,
           parsing, and coordinate mapping
        3. Build the Hail Table from the mapped TSV by invoking the
           ``uniprot_ptm:sites`` Phase B builder via
           :func:`hvantk.core.plugin.run_builder.run_builder_for_spec`.
           The build is stamped with platform Provenance.

    Parameters
    ----------
    config : PTMBuildConfig
        Pipeline configuration.  If ``config.ptm_tsv`` is None the UniProt
        TSV will be downloaded automatically.

    Returns
    -------
    PTMBuildResult
        Mapping statistics, ``mapped_tsv_path``, and ``output_ht``.

    Raises
    ------
    ValueError
        If configuration validation fails.
    """
    from pathlib import Path

    from hvantk.core.plugin import loader as plugin_loader
    from hvantk.core.plugin.run_builder import run_builder_for_spec

    # Step 1: Ensure UniProt PTM TSV is available before calling the core.
    # ptm_build_pipeline_core requires config.ptm_tsv to be set.
    if config.ptm_tsv is None:
        logger.info("Downloading UniProt PTM TSV to %s...", config.output_dir)
        config = PTMBuildConfig(
            output_dir=config.output_dir,
            output_ht=config.output_ht,
            gtf_path=config.gtf_path,
            ptm_tsv=download_uniprot_ptm(config.output_dir, config.overwrite),
            peptideatlas_tsv=config.peptideatlas_tsv,
            cptac_tsv=config.cptac_tsv,
            flanking_codons=config.flanking_codons,
            reference_genome=config.reference_genome,
            overwrite=config.overwrite,
        )

    # Step 2: Pure coordinate-mapping core (no skills).
    result = ptm_build_pipeline_core(config)

    # Step 3: Build Hail Table from the mapped TSV via the Phase B contract.
    if result.n_mapped > 0:
        logger.info("Building Hail Table at %s...", config.output_ht)
        reg = plugin_loader.get_registry()
        spec = reg.get_dataset("uniprot_ptm:sites")
        mapped_path = _resolve_mapped_path(config, result)
        run_builder_for_spec(
            spec,
            parsed_input=mapped_path,
            output_path=Path(config.output_ht),
            plugin_version=spec.plugin_version or "<unknown>",
            reference_genome=config.reference_genome,
            flanking_codons=config.flanking_codons,
        )
        result.output_ht = config.output_ht

    logger.info("PTM build pipeline complete: %d sites mapped", result.n_mapped)
    return result


def _resolve_mapped_path(config: PTMBuildConfig, result: PTMBuildResult) -> str:
    """Return the mapped TSV path produced by :func:`ptm_build_pipeline_core`.

    The core sets ``result.mapped_tsv_path`` to either:

    * ``<output_dir>/ptm_sites_mapped.tsv.bgz`` (single-source UniProt-only run), or
    * ``<output_dir>/ptm_sites_combined.tsv.bgz`` (multi-source run after concat).

    Both cases are captured in ``PTMBuildResult.mapped_tsv_path``, so we simply
    return that field.  If for some reason it is empty we fall back to the
    conventional single-source filename.
    """
    if result.mapped_tsv_path:
        return result.mapped_tsv_path

    import os

    # Fallback: replicate the naming convention from ptm_build_pipeline_core.
    return os.path.join(config.output_dir, "ptm_sites_mapped.tsv.bgz")


__all__ = [
    "download_uniprot_ptm",
    "ptm_build_pipeline",
]
