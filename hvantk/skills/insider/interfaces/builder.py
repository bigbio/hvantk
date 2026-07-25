"""Hail Table builder for the INSIDER per-protein interface reduction.

Companion to ``insider:variants``. That dataset is interval-keyed and answers "does this
variant fall in a PPI interface?"; this one is protein-keyed and answers "how much of
this gene's product is interface, and with how many partners?" -- a gene-level structural
property, usable as a feature axis after ``hvantk annotate prepare`` maps
``uniprot_id`` onto the spine via ``hgnc:lookup``.

Reads ``H_sapiens_interfacesALL.txt`` (~49 MB), not the 1.17 GB BED: the pair table
already carries both accessions and both interface-residue lists, so no genomic join is
needed for a gene-level reduction.
"""
from __future__ import annotations

import logging
import os

logger = logging.getLogger(__name__)

DEFAULT_FILENAME = "H_sapiens_interfacesALL.txt"


def _resolve_input(parsed_input) -> str:
    """Accept either the interface file itself or the raw dir containing it."""
    path = str(parsed_input)
    if os.path.isdir(path):
        return os.path.join(path, DEFAULT_FILENAME)
    return path


def build_insider_interfaces(parsed_input, ctx, **params):
    """Phase B builder -- returns an AnnotationTable keyed by ``uniprot_id``.

    Columns: ``n_partners``, ``n_partners_experimental``, ``n_partners_predicted``,
    ``n_interface_residues``. Interface-residue counts are a UNION across a protein's
    interactions, so they are bounded by protein length rather than by partner count.
    """
    import hail as hl

    from hvantk.core.models import AnnotationTable

    from .parse import parse_interfaces

    df = parse_interfaces(_resolve_input(parsed_input))
    if df.empty:
        raise ValueError(
            "INSIDER interface table produced 0 proteins; check the input file"
        )
    ht = hl.Table.from_pandas(df, key=["uniprot_id"])
    logger.info("build_insider_interfaces: %d proteins", df.shape[0])
    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="insider-interfaces-v1")
    )
