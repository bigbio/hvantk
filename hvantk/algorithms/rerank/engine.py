# local/rerank_engine/engine.py
from dataclasses import dataclass
import numpy as np, pandas as pd
from hvantk.algorithms.cohort.frame import load_cohort_frame
from hvantk.algorithms.rerank.config import validate
from hvantk.algorithms.rerank.features import FeatureAssembler
from hvantk.algorithms.rerank.reranker import ReRanker
from hvantk.algorithms.rerank.tiers import TierAssigner
from hvantk.algorithms.rerank.evaluator import Evaluator, EvalResult


@dataclass
class SelectionSummary:
    """What feature selection did, for one arm of one run.

    ``auc_nested`` is the headline. ``auc_global`` is reported next to it, never instead of
    it: the global pass selects once on all the data, so its estimate is contaminated by
    label exposure.

    Read the gap as a diagnostic, NOT as a bias estimate. It is the difference between two
    estimators that differ in two ways at once -- the global one has seen every label
    (which inflates it), and it also commits every fold to a single feature set instead of
    one fitted per fold (which can help or hurt). Those pull in opposite directions and the
    sum has no guaranteed sign. Measured on the real CHD cohort (1362 genes, 51 features):
    +0.0098 on the `all` arm but -0.0100 on `clean`, from the same data and the same
    policy. A large gap in either direction says the selection is unstable on this cohort;
    a small one says little.
    """

    arm: str
    frequency: dict          # axis -> column -> number of folds that selected it
    global_features: dict    # axis -> columns selected by the global (all-data) pass
    auc_nested: float
    auc_global: float
    n_conflicted: int
    n_unknown: int


@dataclass
class RerankResult:
    table: pd.DataFrame
    metrics: EvalResult
    coverage: dict
    selection: "SelectionSummary | None" = None


def _axis_selector(policy, groups, frequency=None):
    """Build the per-fold selector: selection runs INDEPENDENTLY WITHIN each axis.

    Within-axis is what makes a per-axis delta-AUC well defined -- "axis X's contribution"
    means "X's best subset over the baseline's best subset", regardless of which other
    axes happen to be present. It also stops a 50-column axis from being penalised
    against a 1-column axis purely for carrying redundant columns.

    ``columns`` is whatever slice the caller is scoring (the ablation path passes
    baseline+one axis), so each axis is intersected with it rather than assumed present.
    """
    from hvantk.algorithms.rerank.selection import select_axis

    def selector(X_train, y_train, columns):
        wanted = set(columns)
        kept, claimed = [], set()
        for axis, axis_cols in groups.items():
            present = [c for c in axis_cols if c in wanted]
            if not present:
                continue
            claimed.update(present)
            chosen = select_axis(X_train, y_train, present, policy).kept
            kept.extend(chosen)
            if frequency is not None:
                counts = frequency.setdefault(axis, {})
                for c in chosen:
                    counts[c] = counts.get(c, 0) + 1
        # Columns owned by no axis are passed through rather than dropped: selection is a
        # within-axis operation and has nothing to say about them.
        kept.extend(c for c in columns if c not in claimed)
        return kept

    return selector


def rerank(config, _allowed_columns=None, _arm="all", _n_conflicted=0,
           _n_unknown=0) -> RerankResult:
    validate(config)
    matrix, coverage = FeatureAssembler().assemble(config)
    prior = config.prior.load().rename(columns={"unit": "gene"})
    pos = config.labels.load()
    df = matrix.merge(prior, on="gene", how="left")
    df["y"] = df["gene"].isin(pos).astype(int)
    feat_cols = [c for c in matrix.columns if c != "gene"]
    if not feat_cols:
        raise ValueError(
            "no feature columns in the assembled matrix: every feature axis contributed only "
            "a 'gene' column. Check the feature specs and that each axis table carries numeric "
            "columns."
        )
    # A column the provenance map never mentions is UNDECLARED, which the contract treats
    # exactly like an explicit None: usable, never clean. Counting it here (before the arm
    # restriction) is what stops an axis nobody declared from silently vanishing from the
    # run -- a dropped column would quietly move the headline with nothing to show for it.
    if config.feature_provenance is not None:
        _n_unknown += sum(1 for c in feat_cols if c not in config.feature_provenance)
    if _allowed_columns is not None:
        feat_cols = [c for c in feat_cols if c in _allowed_columns]
        if not feat_cols:
            raise ValueError(
                f"rerank arm {_arm!r}: no feature column survives the provenance filter. "
                f"Every column conflicts with the label provenance "
                f"{sorted(config.label_provenance)!r} or is undeclared, so this arm has "
                "nothing to score. Declare 'trained_on' for the columns that are not "
                "actually derived from those sources, or re-derive the labels."
            )
    for c in feat_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    y = df["y"].values
    covered = sum(1 for g in pos if g in set(matrix["gene"]))
    if pos and covered / len(pos) < config.min_label_coverage:
        raise ValueError(
            f"label/feature join coverage too low: only {covered}/{len(pos)} "
            f"({covered/len(pos):.0%}) positive units are in the feature matrix — likely a gene-symbol/build mismatch"
        )
    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)
    if min(n_pos, n_neg) < config.folds:
        raise ValueError(
            f"too few examples for {config.folds}-fold OOF: {n_pos} positive / {n_neg} negative units "
            f"(need >= {config.folds} of each class). Provide more labels or lower Config.folds."
        )
    # Axis grouping is needed before scoring, not just for the ablation table: selection
    # is a within-axis operation, so the selector has to know which column belongs where.
    groups = {
        ax.name: [c for c in ax.load().columns if c != "gene" and c in feat_cols]
        for ax in config.features
    }
    groups = {k: v for k, v in groups.items() if v}
    baseline = next(iter(groups))
    reranker = ReRanker(config.calibration, config.folds)
    selector = summary = None
    if config.selection is not None:
        frequency: dict = {}
        selector = _axis_selector(config.selection, groups)
        recording = _axis_selector(config.selection, groups, frequency)
        scores = reranker.score(df, feat_cols, y, selector=recording)
        summary = _selection_summary(
            config, df, y, groups, scores, frequency, _arm, _n_conflicted, _n_unknown
        )
    else:
        scores = reranker.score(df, feat_cols, y)
    audit_table = df
    if config.cohort is not None:
        # The prior column was already consumed above (as `prior_stat`), so it is
        # excluded here: re-merging it under its raw name would add nothing new and
        # would make it a spurious collision candidate against a feature axis that
        # legitimately reuses the same column name (e.g. a cohort whose prior is the
        # same p-value a "burden" axis also carries as a model feature).
        cohort_cols = load_cohort_frame(config.cohort, include_prior=False)
        collisions = [
            c for c in cohort_cols.columns if c != "gene" and c in audit_table.columns
        ]
        if collisions:

            def _source(col):
                if col == "prior_stat":
                    return "the prior"
                if col == "y":
                    return "the label column"
                if col in feat_cols:
                    return "a feature axis"
                return "the scored table"

            detail = "; ".join(
                f"{c!r} (already supplied by {_source(c)})" for c in sorted(collisions)
            )
            raise ValueError(
                f"cohort {config.cohort.name!r} declares column(s) that collide with "
                f"the scored table: {detail}. A cohort column must never silently "
                "override, or be silently dropped in favour of, a same-named model/"
                "prior column -- rename the colliding column(s) in the cohort manifest."
            )
        audit_table = df.merge(cohort_cols, on="gene", how="left")
    flag_reason = config.audit.apply(audit_table).reset_index(drop=True)
    flag = flag_reason != ""
    tiers = TierAssigner(config.tiers).assign(scores)  # pure credibility, no flag input
    metrics = Evaluator().evaluate(df, feat_cols, y, scores, groups, baseline, selector)
    table = pd.DataFrame(
        {
            "gene": df.gene,
            "prior_stat": df.prior_stat,
            "score": scores,
            "flag": flag.values,
            "flag_reason": flag_reason.values,
            "y": y,
        }
    )
    table = pd.concat(
        [table.reset_index(drop=True), tiers.reset_index(drop=True)], axis=1
    )
    table["score_percentile"] = (table["score"].rank(pct=True) * 100).round(1)
    # Append extra genes excluded from model scoring (e.g. HCAR1: too few case variants).
    # They are unscored and carry an advisory flag; they are NOT forced onto the tier ladder.
    if config.extra_flagged_genes:
        extra_y = {g: int(g in pos) for g in config.extra_flagged_genes}
        extra_rows = pd.DataFrame(
            [
                {
                    "gene": g,
                    "prior_stat": float("nan"),
                    "score": float("nan"),
                    "score_percentile": float("nan"),
                    "flag": True,
                    "flag_reason": "insufficient_data",
                    "y": extra_y.get(g, 0),
                    "tier": "unscored",
                    "verdict": "unscored",
                }
                for g in config.extra_flagged_genes
                if g not in set(df.gene)
            ]
        )
        if len(extra_rows):
            table = pd.concat([table, extra_rows], ignore_index=True)
    table = table[
        [
            "gene",
            "prior_stat",
            "score",
            "score_percentile",
            "tier",
            "verdict",
            "flag",
            "flag_reason",
            "y",
        ]
    ]
    return RerankResult(
        table=table, metrics=metrics, coverage=coverage, selection=summary
    )


def _selection_summary(config, df, y, groups, scores, frequency, arm,
                       n_conflicted, n_unknown) -> SelectionSummary:
    """Run the global pass and package it with the nested result.

    The global pass exists ONLY to produce a human-readable "these are the features"
    list -- one selection over all the data, which is what a reader can actually inspect
    and argue with. Its AUC is reported alongside as ``auc_global`` so the selection bias
    it carries is visible rather than hidden; nothing downstream ranks on it.
    """
    from sklearn.metrics import roc_auc_score

    from hvantk.algorithms.rerank.selection import select_axis

    global_features = {
        axis: select_axis(df, y, cols, config.selection).kept
        for axis, cols in groups.items()
    }
    picked = [c for cols in global_features.values() for c in cols]
    auc_global = float("nan")
    if picked:
        global_scores = ReRanker(config.calibration, config.folds).score(df, picked, y)
        auc_global = float(roc_auc_score(y, global_scores))
    return SelectionSummary(
        arm=arm,
        frequency=frequency,
        global_features=global_features,
        auc_nested=float(roc_auc_score(y, scores)),
        auc_global=auc_global,
        n_conflicted=n_conflicted,
        n_unknown=n_unknown,
    )


def rerank_arms(config) -> dict:
    """Run rerank once per provenance arm and return ``{arm_name: RerankResult}``.

    Two arms when selection and provenance are both configured:
      clean -- columns with no provenance conflict against this label source. ALWAYS the
               headline: statistical filtering cannot detect circularity, it REWARDS it
               (REVEL correlates with the label partly because it was trained on genes
               like these).
      all   -- clean + conflicted + unknown. Exists only so the circularity channel is a
               measured number instead of an assumption.

    Both arms use identical folds, so the delta between them is paired. An undeclared
    column is usable but never contributes to the headline, so ``all - clean`` bundles the
    circularity channel with whatever undeclared provenance is worth; the summary keeps
    ``n_conflicted`` and ``n_unknown`` separate so the two are not confused.
    """
    from hvantk.algorithms.rerank.provenance import DEFAULT_EQUIVALENCE, resolve_arms

    if config.selection is None or config.feature_provenance is None:
        return {"all": rerank(config)}

    if config.label_provenance is None:
        raise ValueError(
            "feature_provenance is set but label_provenance is not. Circularity is a "
            "property of the feature/label PAIR, so an undeclared label source makes the "
            "clean arm meaningless: nothing conflicts with nothing, every circular "
            "predictor is admitted, and the run looks healthy. Declare what the labels "
            "were derived from (e.g. frozenset({'ClinGen', 'GenCC'})), or pass an "
            "explicit frozenset() to assert they derive from nothing curated."
        )

    assignment = resolve_arms(
        config.feature_provenance,
        config.label_provenance,
        config.provenance_equivalence or DEFAULT_EQUIVALENCE,
    )
    # `all` is left unrestricted rather than set to `assignment.all_columns`: the matrix may
    # carry columns the provenance map never mentions, and those belong in `all` (they are
    # merely undeclared, not disqualified). `rerank` counts them into n_unknown.
    arms = {"clean": set(assignment.clean), "all": None}
    return {
        name: rerank(
            config,
            _allowed_columns=allowed,
            _arm=name,
            _n_conflicted=len(assignment.conflicted),
            _n_unknown=len(assignment.unknown),
        )
        for name, allowed in arms.items()
    }
