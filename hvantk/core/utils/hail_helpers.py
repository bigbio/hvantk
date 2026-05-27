"""Shared Hail Table infrastructure used by Phase B plugin builders.

``create_table_base`` is the common boilerplate for the small set of plugin
builders that still follow the import → optional transform → checkpoint →
optional TSV-export pattern. ``cleanup_temp_file`` provides best-effort
cleanup of local / Hadoop / S3 / GS temp files.

Both functions were previously in ``hvantk/core/builders/table.py`` alongside
the (now retired) Phase A ``create_<x>_tb`` builders; they moved here when
the Phase A surface was retired (issue #114).
"""

from __future__ import annotations

import logging
import os
from typing import Callable, List, Optional

import hail as hl

from hvantk.core.models.metadata import build_table_metadata

logger = logging.getLogger(__name__)
_FILE_URI_PREFIX = "file://"


def create_table_base(
    source_name: str,
    input_path: str,
    output_path: str,
    import_func: Callable[[], hl.Table],
    transform_func: Optional[Callable[[hl.Table], hl.Table]] = None,
    fields: Optional[List[str]] = None,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> hl.Table:
    """Common scaffold for Hail Table builders.

    1. Logs the import.
    2. Runs ``import_func()`` to obtain the table.
    3. Applies ``transform_func`` if provided.
    4. Optionally subsets to ``fields``.
    5. Stamps ``hvantk_metadata`` globals.
    6. Checkpoints to ``output_path``.
    7. Optionally exports a TSV alongside the checkpoint.
    """
    logger.info(f"Creating {source_name} table from {input_path}")
    ht = import_func()

    if transform_func is not None:
        ht = transform_func(ht)

    if fields is not None:
        logger.info(f"Selecting fields: {fields}")
        ht = ht.select(*fields)

    ht = ht.annotate_globals(
        hvantk_metadata=build_table_metadata(source_name, input_path, ht)
    )

    logger.info(f"Checkpointing table to {output_path}")
    ht = ht.checkpoint(output=output_path, overwrite=overwrite)

    if export_tsv:
        logger.info(f"Exporting table to {output_path}.tsv.bgz")
        ht.export(output_path + ".tsv.bgz")

    return ht


def cleanup_temp_file(tmp_path: Optional[str]) -> None:
    """Best-effort cleanup for local or Hadoop/S3/GS temp files."""
    if not tmp_path:
        return
    try:
        import hailtop.fs as hfs

        if hfs.exists(tmp_path):
            if hfs.is_dir(tmp_path):
                hfs.rmtree(tmp_path)
            else:
                hfs.remove(tmp_path)
        return
    except Exception:
        logger.debug(
            "Failed to remove temp path via hailtop.fs: %s", tmp_path, exc_info=True
        )

    try:
        local_path = tmp_path
        if local_path.startswith(_FILE_URI_PREFIX):
            local_path = local_path[len(_FILE_URI_PREFIX):]
        if os.path.exists(local_path):
            os.remove(local_path)
    except Exception:
        logger.debug(
            "Failed to remove temp path via os.remove: %s", tmp_path, exc_info=True
        )
