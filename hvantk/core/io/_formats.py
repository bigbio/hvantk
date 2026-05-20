"""Per-format handlers. Dispatch table lives in __init__.py."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.provenance import Provenance


def save_annotation_table_parquet(ann: AnnotationTable, path: Path) -> None:
    df = ann.to_pandas()
    df.to_parquet(path, index=False)


def load_annotation_table_parquet(path: Path, provenance: Provenance) -> AnnotationTable:
    df = pd.read_parquet(path)
    return AnnotationTable.from_pandas(df, provenance=provenance)


def save_annotation_table_ht(ann: AnnotationTable, path: Path) -> None:
    ht = ann.to_hail()
    ht.write(str(path), overwrite=True)


def load_annotation_table_ht(path: Path, provenance: Provenance) -> AnnotationTable:
    import hail as hl

    ht = hl.read_table(str(path))
    return AnnotationTable.from_hail(ht, provenance=provenance)
