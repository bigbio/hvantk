# local/rerank_engine/catalog/registry.py
from hvantk.algorithms.rerank.config import Config
from hvantk.algorithms.rerank.veto import NoOpVeto, CaseControlArchitectureVeto

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

def _default_veto(profile):
    return CaseControlArchitectureVeto() if profile.cohort is not None else NoOpVeto()

def build_config(profile, axes, *, veto=None, calibration="isotonic", folds=5, tiers=5):
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
        veto=veto if veto is not None else _default_veto(profile),
        extra_vetoed_genes=list(profile.extra_vetoed_genes),
        min_label_coverage=profile.min_label_coverage,
        calibration=calibration, folds=folds, tiers=tiers,
    )
