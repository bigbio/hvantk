# hvantk/algorithms/rerank/catalog/builders.py
"""Generic, path-free axis/label builders for the rerank CLI/library.
Domain-specific builders (gnomAD/Cardoso/GenCC) live in the user's workspace as recipes."""
import pandas as pd
from hvantk.algorithms.rerank.config import FeatureAxis, LabelSpec


def table_axis(name: str, path: str) -> FeatureAxis:
    """A FeatureAxis backed by a gene-keyed table file (.parquet/.tsv/.csv).
    The table must have a 'gene' column + one or more numeric feature columns."""
    def loader():
        if path.endswith(".parquet"):
            return pd.read_parquet(path)
        sep = "\t" if path.endswith((".tsv", ".bgz", ".gz")) else ","
        return pd.read_csv(path, sep=sep)
    return FeatureAxis(name, loader)


def genelist_labels(path: str) -> LabelSpec:
    """A LabelSpec from a file of gene symbols (one per line) or a table's first column."""
    def load():
        with open(path) as fh:
            return {ln.strip().split()[0] for ln in fh if ln.strip() and not ln.startswith("#")}
    return LabelSpec(load)
