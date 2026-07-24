# local/rerank_engine/catalog/registry.py
from hvantk.algorithms.rerank.config import Config
from hvantk.algorithms.rerank.audit import NoAudit, CaseControlArchitectureAudit

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
    return CaseControlArchitectureAudit() if profile.cohort is not None else NoAudit()

def build_config(profile, axes, *, audit=None, calibration="isotonic", folds=5, tiers=5):
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
        calibration=calibration, folds=folds, tiers=tiers,
    )
