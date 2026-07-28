# local/rerank_engine/config.py
import logging
from dataclasses import dataclass, field
from typing import Callable, Optional, TYPE_CHECKING
import pandas as pd

from hvantk.algorithms.cohort.frame import load_prior_frame
from hvantk.algorithms.cohort.spec import CohortManifest

if TYPE_CHECKING:
    from hvantk.algorithms.rerank.audit import Audit
    from hvantk.algorithms.rerank.leakage import LeakagePolicy
    from hvantk.algorithms.rerank.selection import SelectionPolicy

logger = logging.getLogger(__name__)


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
    prior is never overridden. If a cohort manifest is also present and its own prior
    source (table/key_column/prior column) disagrees with the explicit prior's
    (path/unit_col/stat_col), __post_init__ logs a warning naming both sources -- the
    explicit prior still wins, but the discrepancy is no longer silent."""
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
    selection: Optional["SelectionPolicy"] = None
    """Feature-selection policy. None (default) disables selection entirely and reproduces
    the pre-selection code path exactly."""
    leakage: Optional["LeakagePolicy"] = None
    """Presence-leakage control: bar columns whose MISSINGNESS predicts the label.

    Independent of `selection` and of `feature_provenance`. Provenance asks what a
    predictor was TRAINED on; this asks which units it was RUN on. A predictor can pass the
    first and fail the second -- EVE is unsupervised on alignments, so it is correctly
    clean on provenance, while the bare flag "was EVE computed for this gene" scored AUC
    0.716 against a ClinGen/GenCC label in the cohort that motivated this, above the whole
    constraint axis.

    None (default) disables the control and reproduces the previous code path exactly."""
    feature_provenance: Optional[dict] = None
    """column -> frozenset of sources the predictor was trained on, or None if undeclared.
    None for the whole dict means provenance is unavailable: a single 'all' arm is run."""
    label_provenance: Optional[frozenset] = None
    """Sources the LABELS were derived from. Circularity is a property of the pair --
    REVEL against a ClinGen-derived label is badly circular; against a purely
    burden-derived one it is far less so -- so neither half means anything alone.

    None means UNDECLARED and is rejected by `rerank_arms`; an explicit `frozenset()`
    means "these labels derive from nothing curated" and is accepted. The distinction is
    the same one `feature_provenance` draws, and it exists because the failure is silent:
    an empty label source conflicts with nothing, so a forgotten declaration produces a
    `clean` arm that admits every circular predictor and looks entirely healthy."""
    provenance_equivalence: Optional[dict] = None
    """Source-name classes for the circularity check; None uses DEFAULT_EQUIVALENCE.
    `selection.yaml` can override the vocabulary and `load_policy` returns it, so it needs
    somewhere to live -- without this field the override is silently discarded and the
    clean/all split is computed against the defaults, which is a wrong answer rather than
    an error."""

    def __post_init__(self):
        if self.audit is None:
            from hvantk.algorithms.rerank.audit import NoAudit

            self.audit = NoAudit()
        if self.leakage is not None:
            from hvantk.algorithms.rerank.leakage import LeakagePolicy

            # Checked rather than duck-typed: the engine reads `.q` and `.min_auc`, so a
            # bare float or dict here would raise deep inside a per-fold selector, or worse,
            # be swallowed and leave the control silently off while the caller believes it
            # is on. A disabled safety control that looks enabled is the failure to avoid.
            if not isinstance(self.leakage, LeakagePolicy):
                raise TypeError(
                    f"Config.leakage must be a LeakagePolicy or None; got "
                    f"{type(self.leakage).__name__}. To use defaults, pass "
                    f"LeakagePolicy()."
                )
        if self.cohort is not None:
            if self.prior is None:
                self.prior = _ManifestPrior(self.cohort)
            elif isinstance(self.prior, PriorSpec):
                # "Disagree" is defined structurally, not by data content: the explicit
                # PriorSpec's (path, unit_col, stat_col) vs. the (table, key_column,
                # prior.column) the manifest would have derived via _ManifestPrior. This
                # is a cheap comparison of source descriptors -- it never reads either
                # table -- and it is intentionally silent when the two sources happen to
                # describe the same file/columns (the common, correct case in
                # registry.build_config, which always passes the prior explicitly
                # alongside a cohort). The explicit prior always wins either way; this
                # only decides whether that precedence is worth flagging.
                manifest_source = (
                    self.cohort.table,
                    self.cohort.key_column,
                    self.cohort.prior.column,
                )
                explicit_source = (
                    self.prior.path,
                    self.prior.unit_col,
                    self.prior.stat_col,
                )
                if explicit_source != manifest_source:
                    logger.warning(
                        "Config %r: an explicitly-supplied Config.prior (path=%r, "
                        "unit_col=%r, stat_col=%r) disagrees with cohort %r's own prior "
                        "(table=%r, key_column=%r, column=%r); the explicit "
                        "Config.prior wins and the cohort's prior is ignored.",
                        self.name,
                        self.prior.path,
                        self.prior.unit_col,
                        self.prior.stat_col,
                        self.cohort.name,
                        self.cohort.table,
                        self.cohort.key_column,
                        self.cohort.prior.column,
                    )


def validate(config: Config) -> None:
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
    if not config.features:
        raise ValueError("at least one FeatureAxis is required")
    if not config.labels.load():
        raise ValueError(
            "labels are empty: LabelSpec.load() returned no positive units"
        )
