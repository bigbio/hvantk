# local/rerank_engine/config.py
from dataclasses import dataclass, field
from typing import Callable, Optional, TYPE_CHECKING
import pandas as pd

from hvantk.algorithms.cohort.frame import load_prior_frame
from hvantk.algorithms.cohort.spec import CohortManifest

if TYPE_CHECKING:
    from hvantk.algorithms.rerank.audit import Audit


@dataclass
class PriorSpec:
    path: str
    unit_col: str
    stat_col: str

    def load(self) -> pd.DataFrame:
        sep = "\t" if self.path.endswith((".tsv", ".bgz", ".gz")) else ","
        df = pd.read_csv(self.path, sep=sep)
        return df.rename(columns={self.unit_col: "unit", self.stat_col: "prior_stat"})[
            ["unit", "prior_stat"]
        ]


@dataclass
class LabelSpec:
    loader: Callable[[], set]

    def load(self) -> set:
        return set(self.loader())


@dataclass
class FeatureAxis:
    name: str
    loader: Callable[
        [], pd.DataFrame
    ]  # returns DataFrame with a 'gene' column + feature columns
    params: dict = field(default_factory=dict)
    # Cache assumes one logical load per FeatureAxis instance; configs construct fresh instances per run.
    _cached: Optional[pd.DataFrame] = field(
        default=None, init=False, repr=False, compare=False
    )

    def __post_init__(self):
        self._cached = None

    def load(self) -> pd.DataFrame:
        if self._cached is None:
            self._cached = self.loader()
        return self._cached


@dataclass
class _ManifestPrior:
    """PriorSpec-shaped view of a CohortManifest's prior column.

    ``Config.__post_init__`` constructs this automatically when a cohort manifest is
    set and no ``PriorSpec`` was supplied directly (M2: the manifest is now the single
    user-facing declaration -- callers stop authoring a ``PriorSpec`` by hand). Not
    part of the public API.

    Delegates to :func:`hvantk.algorithms.cohort.frame.load_prior_frame` instead of
    re-reading the file with ``PriorSpec.load()``'s own suffix-based delimiter
    heuristic, so this can never silently disagree with
    ``hvantk.algorithms.cohort.checks`` (the same delimiter/gzip logic ``hvantk cohort
    validate``/``attach`` use) about how the cohort table is parsed.
    """

    manifest: CohortManifest

    def load(self) -> pd.DataFrame:
        return load_prior_frame(self.manifest)


@dataclass
class Config:
    name: str
    features: list
    labels: LabelSpec
    units: str = "gene"
    cohort: Optional[CohortManifest] = None
    """External cohort declaration: the single source of the prior statistic and, when
    present, of the case/control architecture columns an Audit reads. Required for
    rerank()/validate() (M3) -- a PriorSpec-only config can still be constructed (e.g.
    by DiseaseProfile/build_config's existing callers) but cannot be scored."""
    prior: Optional[PriorSpec] = None
    """The unit -> prior_stat table engine.rerank merges in before scoring. Derived
    automatically from `cohort` in __post_init__ (via _ManifestPrior) when a cohort
    manifest is present and no prior was given explicitly; an explicitly-supplied
    prior is never overridden."""
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
        if self.cohort is not None and self.prior is None:
            self.prior = _ManifestPrior(self.cohort)


def validate(config: Config) -> None:
    from hvantk.algorithms.rerank.audit import NoAudit

    if config.units != "gene":
        raise NotImplementedError(
            f"units={config.units!r}: only 'gene' is supported in v1"
        )
    if config.cohort is None:
        raise ValueError(
            "rerank requires Config.cohort: a CohortManifest is the single supported "
            "source of the prior statistic (and, when an Audit needs them, the case/"
            "control architecture columns); a PriorSpec-only config is no longer "
            "accepted by rerank. Build a manifest with "
            "hvantk.algorithms.cohort.spec.load_cohort(path) and set it as Config.cohort."
        )
    # Narrower predecessor of the unconditional check above -- kept because it still
    # documents the original, more specific rule (a real Audit always needs a cohort);
    # it can no longer fire on its own since config.cohort is already required above.
    if config.cohort is None and not isinstance(config.audit, NoAudit):
        raise ValueError("a non-NoAudit audit requires a cohort (variant-level data)")
    if not config.features:
        raise ValueError("at least one FeatureAxis is required")
    if not config.labels.load():
        raise ValueError(
            "labels are empty: LabelSpec.load() returned no positive units"
        )
