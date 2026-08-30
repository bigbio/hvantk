"""Legacy artifact loader for pre-plugin Hail Tables.

This module preserves the path-resolution behavior of the now-deleted
``hvantk.core.models.dataset.get_*_ht()`` functions. Each entry in
``_LEGACY_PATHS`` maps a dataset key to its conventional location
relative to ``source_dir``. Callers that have not yet been promoted to
the plugin system (see Phase K of the data-model platform refactor) use
``load_legacy_table(name)`` to read these tables.

The ``source_dir`` module-level variable is set by the user (typically
via CLI configuration or a notebook bootstrap step). When None,
``load_legacy_table`` raises a clear error.

Legacy datasets covered here will be retired as their canonical sources
become plugins (see Phase K and follow-ups). For now they continue to be
read by path with an unknown-provenance marker (see ``core/io/_legacy.py``)
attached when wrapped in an AnnotationTable.
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    import hail as hl  # noqa: F401 — referenced from string annotations

logger = logging.getLogger(__name__)


# Module-level config — set externally by the user before calling load_legacy_table.
source_dir: Optional[str] = None


# Dataset key → relative path under source_dir.
_LEGACY_PATHS = {
    "chd_denovo": "data/ht/DNM_Jin2017_Sifrim2016_GRCh38_lift.ht",
    "clinvar": "data/ht/clinvar.GRCh38.ht",
    "gene_expression": "data/ht/rnaseq.human.ht",
    "gene_ann": "data/ht/gene.ann.ensembl.ht",
    "gevir": "data/ht/gevir.metrics.ht",
    "ppi": "data/ht/interactome.GRCh38.ht",
    "dbnsfp_scores": "data/ht/dbNSFP4.1a_variant.ht",
    "gnomad_metrics": "data/ht/gnomad.metrics.ht",
    "gnomad_af": "data/ht/gnomad_3.0_sites_AF.ht",
    "deg": "data/ht/scell.heart.degs.ht",
    "hca": "data/ht/hca.heart.ht",
}


# Path to the gene-set TSV file (not a Hail Table).


def _require_source_dir(provided: Optional[str]) -> str:
    """Resolve source_dir from arg or module global; raise if neither set."""
    sd = provided if provided is not None else source_dir
    if not sd:
        raise ValueError(
            "hvantk.core.io.legacy_artifacts.source_dir is not set; "
            "either pass source_dir explicitly to load_legacy_table() or "
            "set the module global before calling."
        )
    return sd


def load_legacy_table(name: str, source_dir: Optional[str] = None) -> "hl.Table":
    """Read a legacy Hail Table by dataset key.

    Parameters
    ----------
    name : str
        One of the keys in _LEGACY_PATHS (e.g. "gevir", "dbnsfp_scores").
    source_dir : str, optional
        Explicit source dir. If None, falls back to the module-level
        ``source_dir`` global. Raises ValueError if neither is set.

    Returns
    -------
    hl.Table
        The Hail Table at <source_dir>/<relative-path>.
    """
    import hail as hl

    if name not in _LEGACY_PATHS:
        raise KeyError(
            f"Unknown legacy dataset: {name!r}. "
            f"Known keys: {sorted(_LEGACY_PATHS)}"
        )
    sd = _require_source_dir(source_dir)
    rel = _LEGACY_PATHS[name]
    path = f"{sd}/{rel}"
    logger.info("Reading legacy table %s from %s", name, path)
    return hl.read_table(path)


def load_legacy_gene_expression_table(
    organ: str = "Heart",
    tp_col: str = "mean_expr_time_point",
    source_dir: Optional[str] = None,
) -> "hl.Table":
    """Specialized loader for the gene expression table.

    The base table is keyed by gene with a struct-typed time-point map; this
    helper expands the time-point struct into per-time-point columns prefixed
    by organ name. This is the same logic the old
    ``hvantk.core.models.dataset.get_gene_expression_ht`` performed.
    """
    import hail as hl

    t = load_legacy_table("gene_expression", source_dir=source_dir)

    # Get available time points for the specified organ
    tps = (
        t[tp_col].key_set().filter(lambda x: x.organ == organ).time_point.collect()[0]
    )

    # Annotate expression values per time point
    t = t.annotate(
        **{
            f"{organ}.{tp}": t[tp_col].get(hl.struct(organ=organ, time_point=tp))
            for tp in tps
        }
    )
    t = t.drop(t["mean_expr_time_point"], t["mean_expr_dev_stage"]).key_by("Gene")
    return t
