"""Evaluate a `hvantk rerank` ranked table against its own `y` labels.

No Hail: pandas + scikit-learn only, so this runs inside the image (through
run-in-image.sh for a quick check, or as a job -- it needs no JVM either way).

    python evaluate-ranking.py --ranked ranked.tsv --receipt eval.json \\
        [--top-decile-pct 90]

Computes exactly what `hvantk rerank` itself does not print: there is no overall-AUC
line on its console (only the per-axis ablation table, which is a different,
uncalibrated quantity -- skills/rerank-evaluate/references/feature-axes.md). Prints
`auc: <value>` as its own line, meant to be copied verbatim into a run's
`run-notes.md`. Exits 1 if the table is missing a required column, or has no
label-positive or no label-negative scored row -- there is nothing to compute an
AUC/average-precision over.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

REQUIRED_COLUMNS = (
    "gene",
    "prior_stat",
    "score",
    "score_percentile",
    "tier",
    "verdict",
    "flag",
    "flag_reason",
    "y",
)


def _to_bool(series):
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin({"true", "1"})


def evaluate(ranked_path: Path, top_decile_pct: float) -> dict:
    import pandas as pd
    from sklearn.metrics import average_precision_score, roc_auc_score

    df = pd.read_csv(ranked_path, sep="\t")
    missing = sorted(set(REQUIRED_COLUMNS) - set(df.columns))
    if missing:
        raise ValueError(f"{ranked_path}: missing column(s) {missing}")
    df["flag"] = _to_bool(df["flag"])
    # An empty flag_reason ("" for a non-flagged gene, exactly what rerank writes)
    # round-trips through pandas' default NA tokens as NaN; restore it.
    df["flag_reason"] = df["flag_reason"].fillna("")

    scored = df[df["score"].notna()].copy()
    y = scored["y"].astype(int)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos == 0:
        raise ValueError(
            f"{ranked_path}: no label-positive (y=1) row among {len(scored)} scored "
            "gene(s)"
        )
    if n_neg == 0:
        raise ValueError(
            f"{ranked_path}: no label-negative (y=0) row among {len(scored)} scored "
            "gene(s)"
        )

    score = scored["score"].astype(float)
    auc = float(roc_auc_score(y, score))
    ap = float(average_precision_score(y, score))
    top_decile = scored[scored["score_percentile"].astype(float) >= top_decile_pct]
    n_top_decile_positive = int((top_decile["y"] == 1).sum())

    tier_counts = {str(k): int(v) for k, v in scored["tier"].value_counts().items()}
    verdict_counts = {
        str(k): int(v) for k, v in scored["verdict"].value_counts().items()
    }
    flagged = scored[scored["flag"]]
    flag_counts = {
        str(k): int(v) for k, v in flagged["flag_reason"].value_counts().items()
    }

    return {
        "ranked": str(ranked_path),
        "n_genes": int(len(scored)),
        "n_positive": n_pos,
        "n_negative": n_neg,
        "auc": auc,
        "average_precision": ap,
        "top_decile_pct": top_decile_pct,
        "n_top_decile_positive": n_top_decile_positive,
        "tier_counts": tier_counts,
        "verdict_counts": verdict_counts,
        "flag_counts": flag_counts,
        "ranked_sha256": hashlib.sha256(ranked_path.read_bytes()).hexdigest(),
    }


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--ranked", required=True)
    ap.add_argument("--receipt", required=True)
    ap.add_argument(
        "--top-decile-pct",
        type=float,
        default=90.0,
        help="score_percentile threshold counted as the top decile [default: 90]",
    )
    args = ap.parse_args()

    try:
        result = evaluate(Path(args.ranked), args.top_decile_pct)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    out_path = Path(args.receipt)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2))
    print(
        f"n_genes={result['n_genes']} n_positive={result['n_positive']} "
        f"n_negative={result['n_negative']} auc={result['auc']:.4f} "
        f"average_precision={result['average_precision']:.4f} "
        f"top_decile_positive={result['n_top_decile_positive']}"
    )
    print(f"auc: {result['auc']:.4f}")
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
