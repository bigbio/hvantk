# local/rerank_engine/features.py
import functools, pandas as pd

class FeatureAssembler:
    def assemble(self, config):
        frames = [ax.load() for ax in config.features]
        matrix = functools.reduce(lambda l, r: l.merge(r, on="gene", how="outer"), frames)
        coverage = {}
        for ax, fr in zip(config.features, frames):
            present = matrix["gene"].isin(set(fr["gene"]))
            coverage[ax.name] = round(float(present.mean()), 3)
        return matrix, coverage
