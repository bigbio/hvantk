# local/rerank_engine/config.py
from dataclasses import dataclass, field
from typing import Callable, Optional, TYPE_CHECKING
import pandas as pd

if TYPE_CHECKING:
    from hvantk.algorithms.rerank.audit import Audit

@dataclass
class PriorSpec:
    path: str; unit_col: str; stat_col: str
    def load(self) -> pd.DataFrame:
        sep = "\t" if self.path.endswith((".tsv",".bgz",".gz")) else ","
        df = pd.read_csv(self.path, sep=sep)
        return df.rename(columns={self.unit_col:"unit", self.stat_col:"prior_stat"})[["unit","prior_stat"]]

@dataclass
class LabelSpec:
    loader: Callable[[], set]
    def load(self) -> set: return set(self.loader())

@dataclass
class FeatureAxis:
    name: str
    loader: Callable[[], pd.DataFrame]   # returns DataFrame with a 'gene' column + feature columns
    params: dict = field(default_factory=dict)
    # Cache assumes one logical load per FeatureAxis instance; configs construct fresh instances per run.
    _cached: Optional[pd.DataFrame] = field(default=None, init=False, repr=False, compare=False)
    def __post_init__(self):
        self._cached = None
    def load(self) -> pd.DataFrame:
        if self._cached is None:
            self._cached = self.loader()
        return self._cached

@dataclass
class CohortSpec:
    """DEPRECATED: superseded by hvantk.algorithms.cohort.spec.CohortManifest.

    This loads a fixed set of audit columns from a gene-keyed table and is consumed in
    exactly one place (engine.py, after scoring). It is retained because
    examples/rerank/config.yaml exposes it publicly and the committed reference output
    in hvantk/tests/rerank/test_example.py asserts flags derived from it. Migrating the
    rerank path onto CohortManifest is a separate increment with its own regression gate.
    """

    variant_table_path: str
    params: dict = field(default_factory=dict)
    def load(self) -> pd.DataFrame:
        cols = self.params.get("audit_cols", [])
        key = self.params.get("key", "gene")
        sep = "\t" if self.variant_table_path.endswith((".tsv", ".bgz", ".gz")) else ","
        df = pd.read_csv(self.variant_table_path, sep=sep)
        keep = [key] + [c for c in cols if c in df.columns]
        return df[keep].drop_duplicates(key).rename(columns={key: "gene"})

@dataclass
class Config:
    name: str
    prior: PriorSpec
    features: list
    labels: LabelSpec
    units: str = "gene"
    cohort: Optional[CohortSpec] = None
    audit: Optional["Audit"] = None
    calibration: str = "isotonic"
    folds: int = 5
    tiers: int = 5
    extra_flagged_genes: list = field(default_factory=list)
    """Genes to append to the output table after scoring, forced FLAGGED/unscored (not in model matrix)."""
    min_label_coverage: float = 0.5
    """Minimum fraction of label-positive units that must appear in the feature matrix.
    Set to 0.0 for intentional cross-disease transfer configs where labels come from a
    different gene universe (e.g. NDD labels scored against a CHD feature matrix)."""
    def __post_init__(self):
        if self.audit is None:
            from hvantk.algorithms.rerank.audit import NoAudit
            self.audit = NoAudit()

def validate(config: Config) -> None:
    from hvantk.algorithms.rerank.audit import NoAudit
    if config.units != "gene":
        raise NotImplementedError(f"units={config.units!r}: only 'gene' is supported in v1")
    if config.cohort is None and not isinstance(config.audit, NoAudit):
        raise ValueError("a non-NoAudit audit requires a cohort (variant-level data)")
    if not config.features:
        raise ValueError("at least one FeatureAxis is required")
    if not config.labels.load():
        raise ValueError("labels are empty: LabelSpec.load() returned no positive units")
