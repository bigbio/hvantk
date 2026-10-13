from hvantk.algorithms.rerank.config import Config
from hvantk.algorithms.rerank.audit import (
    NoAudit,
    CaseControlArchitectureAudit,
    has_architecture_columns,
)
from hvantk.algorithms.rerank.seeds import DEFAULT_SEED

AXES = {}


def axis(name):
    def deco(fn):
        AXES[name] = fn
        return fn

    return deco


def _constraint_first(feats):
    con = [f for f in feats if f.name == "constraint"]
    rest = [f for f in feats if f.name != "constraint"]
    return con + rest


def _default_audit(profile):
    # Must agree with hvantk.tools.rerank.rerank_cli's audit-selection: a cohort is
    # only eligible for CaseControlArchitectureAudit when it actually declares the
    # three columns the audit needs (has_architecture_columns), not merely by being
    # present. Since Config.cohort is mandatory for every rerank config, "cohort
    # is not None" is true for every config -- so it can no longer stand in for "has
    # the architecture columns" the way it once did.
    #
    # The eligibility check is against axis_columns(), not declared_columns():
    # engine.rerank() merges the cohort frame with include_prior=False, so the prior
    # column never reaches the audit table. Testing declared_columns() (which still
    # carries the prior column) can green-light a manifest whose prior happens to be
    # named e.g. "n_case_var" even though the engine will never actually merge that
    # column in -- the audit would then raise instead of running, on a manifest that
    # looked eligible.
    if profile.cohort is None:
        return NoAudit()
    if has_architecture_columns(profile.cohort.axis_columns()):
        return CaseControlArchitectureAudit()
    return NoAudit()


def build_config(
    profile,
    axes,
    *,
    audit=None,
    calibration="isotonic",
    folds=5,
    tiers=5,
    seed=DEFAULT_SEED,
    seed_sweep=1,
    nulls=None,
    blocks=None,
):
    """Assemble a ``Config`` from a ``DiseaseProfile`` and a list of axis names/instances.

    ``audit``/``calibration``/``folds``/``tiers`` were the original surface; ``seed``,
    ``seed_sweep``, ``nulls`` and ``blocks`` (#390/#408) are later additions threaded
    straight through to the matching ``Config`` field with no transformation. Every one of
    them defaults to that field's OWN default (``seed=DEFAULT_SEED``, ``seed_sweep=1``,
    ``nulls=None``, ``blocks=None``) rather than a value re-typed here, so a caller who
    builds a ``Config`` through this registry and one who constructs ``Config`` directly
    start from the identical, single source of truth for what "default" means.

    seed
        The one CV/estimator/bootstrap/null seed (``Config.seed``).
    seed_sweep
        How many CV seeds the ablation interval is swept over (``Config.seed_sweep``).
    nulls
        Permutation-null settings, or ``None`` to skip the multiplicity correction
        (``Config.nulls``, a ``NullConfig`` when set).
    blocks
        Paralogue-blocked cross-validation policy, or ``None`` for unblocked folds
        (``Config.blocks``, a ``BlockPolicy`` when set).

    None of the four is validated here: ``Config.__post_init__`` already rejects a
    wrong-typed ``nulls``/``blocks``/``seed``/``seed_sweep``, and duplicating that check in
    two places would only give it two chances to drift apart.
    """
    if profile.prior is None:
        raise ValueError("DiseaseProfile.prior is required for rerank (build_config)")
    feats = [AXES[a](profile) if isinstance(a, str) else a for a in axes]
    feats = _constraint_first(feats)
    return Config(
        name=profile.name,
        prior=profile.prior,
        features=feats,
        labels=AXES["labels"](profile),
        cohort=profile.cohort,
        audit=audit if audit is not None else _default_audit(profile),
        extra_flagged_genes=list(profile.extra_flagged_genes),
        min_label_coverage=profile.min_label_coverage,
        calibration=calibration,
        folds=folds,
        tiers=tiers,
        seed=seed,
        seed_sweep=seed_sweep,
        nulls=nulls,
        blocks=blocks,
    )
