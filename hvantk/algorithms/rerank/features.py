# hvantk/algorithms/rerank/features.py
import functools
import logging

import pandas as pd

logger = logging.getLogger(__name__)


class FeatureAssembler:
    def assemble(self, config):
        frames = []
        for ax in config.features:
            fr = ax.load()
            if fr["gene"].duplicated().any():
                n_dup = int(fr["gene"].duplicated().sum())
                logger.warning(
                    "rerank feature axis '%s': %d duplicate gene key(s) dropped (kept first) "
                    "to avoid a cartesian blow-up on the outer join", ax.name, n_dup)
                fr = fr.drop_duplicates("gene", keep="first")
            frames.append(fr)
        matrix = functools.reduce(lambda l, r: l.merge(r, on="gene", how="outer"), frames)
        coverage = {}
        for ax, fr in zip(config.features, frames):
            present = matrix["gene"].isin(set(fr["gene"]))
            coverage[ax.name] = round(float(present.mean()), 3)
        return matrix, coverage
