# hvantk/algorithms/rerank/features.py
import functools
import logging

import pandas as pd

logger = logging.getLogger(__name__)


class FeatureAssembler:
    def assemble(self, config):
        frames = []
        col_owner = {}  # non-gene feature column -> the axis that first declared it
        for ax in config.features:
            fr = ax.load()
            if fr["gene"].duplicated().any():
                n_dup = int(fr["gene"].duplicated().sum())
                logger.warning(
                    "rerank feature axis '%s': %d duplicate gene key(s) dropped (kept first) "
                    "to avoid a cartesian blow-up on the outer join",
                    ax.name,
                    n_dup,
                )
                fr = fr.drop_duplicates("gene", keep="first")
            # Feature column names must be unique across axes: the outer join below would
            # otherwise suffix a shared name to `_x`/`_y`, which both feeds the model twice
            # and silently breaks the per-axis grouping in engine.py (it matches columns by
            # exact name). Fail loud instead. (PR #222 review.)
            clashes = [c for c in fr.columns if c != "gene" and c in col_owner]
            if clashes:
                raise ValueError(
                    f"rerank feature axis '{ax.name}' declares column(s) {clashes} already "
                    f"provided by axis '{col_owner[clashes[0]]}'; feature column names must be "
                    f"unique across axes. Rename the columns in the axis tables."
                )
            for c in fr.columns:
                if c != "gene":
                    col_owner[c] = ax.name
            frames.append(fr)
        matrix = functools.reduce(
            lambda l, r: l.merge(r, on="gene", how="outer"), frames
        )
        coverage = {}
        for ax, fr in zip(config.features, frames):
            present = matrix["gene"].isin(set(fr["gene"]))
            coverage[ax.name] = round(float(present.mean()), 3)
        return matrix, coverage
