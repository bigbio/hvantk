"""Builder for `onek-genomes:samples` — imports IGSR's canonical sample-metadata TSV."""
from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


def build_onek_genomes_samples(parsed_input, ctx, **params):
    """Phase B builder — returns an AnnotationTable of IGSR sample metadata.

    Parameters
    ----------
    parsed_input : str | Path
        Path to a directory containing ``igsr_samples.tsv`` (populated by the
        ``download_igsr_samples`` lifecycle stage).
    ctx : hvantk.core.models.BuildContext
        Platform-provided context.
    """
    import hail as hl
    from pathlib import Path
    from hvantk.core.models import AnnotationTable

    raw_dir = Path(parsed_input)
    samples_path = raw_dir / "igsr_samples.tsv"
    if not samples_path.exists():
        raise FileNotFoundError(
            f"Expected igsr_samples.tsv under {raw_dir}; run with --skip-download=false "
            "or place the file there manually."
        )
    logger.info("Importing IGSR samples from %s", samples_path)
    # The panel TSV has trailing tabs in the header; use impute=False with explicit types
    # to avoid Hail inferring extra columns from those tabs.
    ht = hl.import_table(
        str(samples_path),
        impute=False,
        types={"sample": "str", "pop": "str", "super_pop": "str", "gender": "str"},
        key="sample",
    )
    # Keep only the named columns (discard unnamed columns from trailing tabs)
    ht = ht.select(pop=ht.pop, super_pop=ht.super_pop, gender=ht.gender)
    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="onek-genomes-samples-v1"),
    )
