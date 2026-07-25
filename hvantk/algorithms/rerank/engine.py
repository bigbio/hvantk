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
class RerankResult:
    table: pd.DataFrame
    metrics: EvalResult
    coverage: dict


def rerank(config) -> RerankResult:
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
    scores = ReRanker(config.calibration, config.folds).score(df, feat_cols, y)
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
    groups = {
        ax.name: [c for c in ax.load().columns if c != "gene" and c in feat_cols]
        for ax in config.features
    }
    groups = {k: v for k, v in groups.items() if v}
    baseline = next(iter(groups))
    metrics = Evaluator().evaluate(df, feat_cols, y, scores, groups, baseline)
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
    return RerankResult(table=table, metrics=metrics, coverage=coverage)
