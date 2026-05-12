"""Stratified constraint analysis at PTM codons.

Compares gnomAD allele-frequency distributions between PTM-proximal and
non-PTM variants, stratified by tissue, cell type, or any categorical
metadata field derived from an expression dataset.

This module is an **orchestrator**, not a per-variant scorer. For per-site
PTM flags use :func:`hvantk.ptm.annotate.annotate_variants_with_ptm`.
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats

from hvantk.ptm.constraint_expression import load_gene_by_group_matrix
from hvantk.utils.tissue_specificity import compute_specificity

logger = logging.getLogger(__name__)

__all__ = [
    "PTMConstraintConfig",
    "PTMConstraintResult",
    "run_ptm_constraint",
]

_DEFAULT_PTM_CATEGORIES = (
    "phosphorylation",
    "acetylation",
    "methylation",
    "ubiquitination",
    "sumoylation",
    "glycosylation",
    "other",
)


@dataclass
class PTMConstraintConfig:
    """Configuration for :func:`run_ptm_constraint`."""

    variants_ht_path: str
    expression_source: str
    expression_path: str
    grouping: str
    output_dir: str

    label_filter: str = "TN"
    label_field: str = "rf_label"
    gene_field: str = "gene_symbol"
    af_field: str = "gnomad_af_genomes"
    loeuf_field: str = "loeuf"
    ptm_category_field: str = "ptm_types"
    gene_id_mapping: Optional[str] = None
    expression_metric: str = "median"
    min_cells_per_group: int = 50
    min_variants_per_group: int = 20
    flanking_codons: int = 7
    expressed_threshold: float = 1.0
    overwrite: bool = False

    def validate(self) -> List[str]:
        errors: List[str] = []
        if self.expression_source not in {"hail-mt", "anndata", "tabular"}:
            errors.append(
                f"--expression-source must be one of 'hail-mt', 'anndata', 'tabular' "
                f"(got '{self.expression_source}')."
            )
        if self.label_filter not in {"TN", "TP", "all"}:
            errors.append(
                f"--label-filter must be one of 'TN', 'TP', 'all' "
                f"(got '{self.label_filter}')."
            )
        if not os.path.exists(self.variants_ht_path):
            errors.append(f"--variants-ht path does not exist: {self.variants_ht_path}")
        if not os.path.exists(self.expression_path):
            errors.append(
                f"--expression-path does not exist: {self.expression_path}"
            )
        if self.gene_id_mapping and not os.path.exists(self.gene_id_mapping):
            errors.append(
                f"--gene-id-mapping does not exist: {self.gene_id_mapping}"
            )
        if self.min_cells_per_group < 0:
            errors.append("--min-cells-per-group must be >= 0.")
        if self.min_variants_per_group < 1:
            errors.append("--min-variants-per-group must be >= 1.")
        return errors


@dataclass
class PTMConstraintResult:
    """Structured outputs of :func:`run_ptm_constraint`."""

    n_variants: int
    n_variants_ptm: int
    n_variants_non_ptm: int
    n_groups: int
    n_genes: int
    top_groups: List[Dict[str, Any]] = field(default_factory=list)
    tau_quartile: List[Dict[str, Any]] = field(default_factory=list)
    loeuf_factorial: List[Dict[str, Any]] = field(default_factory=list)
    within_gene_paired: Dict[str, Any] = field(default_factory=dict)
    category_heatmap: List[Dict[str, Any]] = field(default_factory=list)
    output_dir: str = ""
    config: Dict[str, Any] = field(default_factory=dict)

    def summary(self) -> str:
        top = (
            f" ({self.top_groups[0]['group']}: {self.top_groups[0]['ratio']:.1f}x)"
            if self.top_groups
            else ""
        )
        return (
            f"PTM constraint analysis: {self.n_variants} variants, "
            f"{self.n_groups} groups, {self.n_genes} genes{top}"
        )


def run_ptm_constraint(config: PTMConstraintConfig) -> PTMConstraintResult:
    """Execute the five stratified constraint tests end-to-end.

    Parameters
    ----------
    config
        Fully-validated :class:`PTMConstraintConfig`.

    Returns
    -------
    PTMConstraintResult
        Top-line numbers; full per-test tables written to
        ``config.output_dir``.
    """
    errors = config.validate()
    if errors:
        raise ValueError("Invalid config: " + "; ".join(errors))

    os.makedirs(config.output_dir, exist_ok=True)
    os.makedirs(os.path.join(config.output_dir, "plots"), exist_ok=True)

    variants_df = _load_variants(config)
    logger.info(
        "Loaded %d variants (ptm_site=%d, proximal=%d)",
        len(variants_df),
        int(variants_df["is_ptm_site"].sum()),
        int(variants_df["is_ptm_proximal"].sum()),
    )

    expr_wide = load_gene_by_group_matrix(
        source=config.expression_source,
        path=config.expression_path,
        grouping=config.grouping,
        aggfunc=config.expression_metric,
        min_cells_per_group=config.min_cells_per_group,
    )

    gene_id_map = _load_gene_id_map(config.gene_id_mapping)
    gene_features = _compute_gene_features(
        expr_wide,
        expressed_threshold=config.expressed_threshold,
        gene_id_map=gene_id_map,
    )

    merged = _join_and_filter(variants_df, gene_features, config)

    logger.info(
        "Post-join: %d variants across %d genes, %d groups",
        len(merged),
        merged["gene_key"].nunique(),
        merged["primary_group"].nunique(),
    )

    per_group = _test_per_group_ranking(merged, config)
    tau_q = _test_tau_stratification(merged, config)
    loeuf_fx = _test_loeuf_group_factorial(merged, config)
    category_hm = _test_category_group_heatmap(merged, config)
    within_gene = _test_within_gene_paired(merged)

    tables = {
        "per_group_ranking": per_group,
        "tau_quartile": tau_q,
        "loeuf_factorial": loeuf_fx,
        "category_group_heatmap": category_hm,
        "within_gene_paired": within_gene["per_gene"],
    }
    _write_tables(tables, config.output_dir)

    result = PTMConstraintResult(
        n_variants=len(merged),
        n_variants_ptm=int(merged["ptm_any"].sum()),
        n_variants_non_ptm=int((~merged["ptm_any"]).sum()),
        n_groups=merged["primary_group"].nunique(),
        n_genes=merged["gene_key"].nunique(),
        top_groups=per_group.head(10).to_dict(orient="records"),
        tau_quartile=tau_q.to_dict(orient="records"),
        loeuf_factorial=loeuf_fx.to_dict(orient="records"),
        within_gene_paired=within_gene["summary"],
        category_heatmap=category_hm.to_dict(orient="records"),
        output_dir=config.output_dir,
        config=asdict(config),
    )

    summary_path = os.path.join(config.output_dir, "summary.json")
    with open(summary_path, "w") as f:
        json.dump(asdict(result), f, indent=2, default=_json_default)
    logger.info("Wrote summary to %s", summary_path)

    try:
        from hvantk.ptm.constraint_plots import render_panels
        from hvantk.ptm.constraint_report import render_html

        render_panels(tables, config.output_dir)
        render_html(result, config, config.output_dir)
    except Exception:
        logger.exception("Plot/report generation failed; continuing with tables only.")

    return result


def _json_default(obj: Any) -> Any:
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj) if np.isfinite(obj) else None
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (pd.Timestamp,)):
        return obj.isoformat()
    raise TypeError(f"Unserialisable type: {type(obj)}")


def _load_variants(config: PTMConstraintConfig) -> pd.DataFrame:
    """Load the PTM-annotated variant HT into pandas, keeping the needed columns."""
    from hvantk.core.hail_context import hl, init_hail

    init_hail()
    ht = hl.read_table(config.variants_ht_path)
    row_fields = set(ht.row)

    required = {"is_ptm_site", "is_ptm_proximal"}
    missing = required - row_fields
    if missing:
        raise KeyError(
            f"Variants HT missing required PTM fields: {missing}. "
            "Run `hvantk ptm annotate` first."
        )

    select: Dict[str, Any] = {
        "is_ptm_site": ht.is_ptm_site,
        "is_ptm_proximal": ht.is_ptm_proximal,
    }
    if config.gene_field in row_fields:
        select["gene_key"] = ht[config.gene_field]
    else:
        raise KeyError(
            f"Gene field '{config.gene_field}' not in variants HT. "
            f"Available row fields: {sorted(row_fields)[:30]}"
        )

    if config.af_field in row_fields:
        select["af"] = ht[config.af_field]
    else:
        logger.warning(
            "AF field '%s' not found; filling with zeros.", config.af_field
        )
        select["af"] = hl.float64(0.0)

    if config.loeuf_field in row_fields:
        select["loeuf"] = ht[config.loeuf_field]
    else:
        select["loeuf"] = hl.missing(hl.tfloat64)

    if config.label_field in row_fields:
        select["label"] = ht[config.label_field]
    else:
        select["label"] = hl.missing(hl.tstr)

    if config.ptm_category_field in row_fields:
        select["ptm_categories"] = ht[config.ptm_category_field]
    else:
        select["ptm_categories"] = hl.missing(hl.tset(hl.tstr))

    ht = ht.select(**select)

    if config.label_filter != "all":
        if config.label_field in row_fields:
            ht = ht.filter(ht.label == config.label_filter)
        else:
            raise ValueError(
                f"--label-filter '{config.label_filter}' requested but label "
                f"field '{config.label_field}' is absent from the variants HT. "
                "Pass --label-filter all to disable filtering or supply the "
                "correct --label-field."
            )

    df = ht.to_pandas()

    df["is_ptm_site"] = df["is_ptm_site"].fillna(False).astype(bool)
    df["is_ptm_proximal"] = df["is_ptm_proximal"].fillna(False).astype(bool)
    df["ptm_any"] = df["is_ptm_site"] | df["is_ptm_proximal"]
    df["af"] = pd.to_numeric(df["af"], errors="coerce").fillna(0.0)
    df["loeuf"] = pd.to_numeric(df["loeuf"], errors="coerce")
    df["gene_key"] = df["gene_key"].astype(str).str.replace(r"\.\d+$", "", regex=True)

    return df


def _load_gene_id_map(path: Optional[str]) -> Optional[Dict[str, str]]:
    if path is None:
        return None
    df = pd.read_csv(path, sep="\t")
    if df.shape[1] < 2:
        raise ValueError(
            "Gene ID mapping TSV needs at least two columns (source, symbol)."
        )
    src, dst = df.columns[0], df.columns[1]
    return dict(zip(df[src].astype(str), df[dst].astype(str)))


def _compute_gene_features(
    expr_wide: pd.DataFrame,
    expressed_threshold: float,
    gene_id_map: Optional[Dict[str, str]],
) -> pd.DataFrame:
    """Produce per-gene τ, primary_group, and max-expression flag."""
    index = expr_wide.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    expr_wide = expr_wide.copy()
    expr_wide.index = index

    if gene_id_map is not None:
        mapped = expr_wide.index.map(lambda x: gene_id_map.get(x, x))
        expr_wide.index = pd.Index(mapped, name="gene_id")

    expr_wide = expr_wide.groupby(level=0).max()

    tau = compute_specificity(expr_wide, method="tau")

    primary_group = expr_wide.idxmax(axis=1)
    max_expr = expr_wide.max(axis=1)

    features = pd.DataFrame(
        {
            "tau": tau.values,
            "primary_group": primary_group.values,
            "max_expr": max_expr.values,
            "expressed": (max_expr >= expressed_threshold).values,
        },
        index=expr_wide.index,
    )
    features.index.name = "gene_key"
    return features


def _join_and_filter(
    variants_df: pd.DataFrame,
    gene_features: pd.DataFrame,
    config: PTMConstraintConfig,
) -> pd.DataFrame:
    """Inner-join variants with gene features; keep expressed genes only."""
    merged = variants_df.merge(
        gene_features, how="inner", left_on="gene_key", right_index=True
    )
    if len(merged) == 0:
        raise ValueError(
            "No variants matched expression data. Check that the gene identifiers "
            "align between variants HT and expression matrix (symbols vs Ensembl IDs)."
        )

    merged = merged[merged["expressed"]].copy()
    if len(merged) == 0:
        raise ValueError(
            f"No variants remain after filtering to expressed genes "
            f"(threshold={config.expressed_threshold})."
        )

    return merged


def _mw_log2_ratio(
    afs_ptm: np.ndarray,
    afs_non: np.ndarray,
    pseudocount: float = 1e-8,
) -> Dict[str, float]:
    """Mann-Whitney U + effect sizes used by every per-strata test."""
    n_ptm = int(np.size(afs_ptm))
    n_non = int(np.size(afs_non))

    med_ptm = float(np.median(afs_ptm)) if n_ptm else np.nan
    med_non = float(np.median(afs_non)) if n_non else np.nan
    mean_ptm = float(np.mean(afs_ptm)) if n_ptm else np.nan
    mean_non = float(np.mean(afs_non)) if n_non else np.nan

    ratio = (
        (mean_non + pseudocount) / (mean_ptm + pseudocount)
        if n_ptm and n_non
        else np.nan
    )
    log2_ratio = float(np.log2(ratio)) if np.isfinite(ratio) and ratio > 0 else np.nan

    pvalue = np.nan
    if n_ptm >= 3 and n_non >= 3:
        try:
            _, pvalue = stats.mannwhitneyu(
                afs_non, afs_ptm, alternative="greater"
            )
            pvalue = float(pvalue)
        except ValueError:
            pvalue = np.nan

    return {
        "n_ptm": n_ptm,
        "n_non": n_non,
        "med_ptm": med_ptm,
        "med_non": med_non,
        "mean_ptm": mean_ptm,
        "mean_non": mean_non,
        "ratio": float(ratio) if np.isfinite(ratio) else np.nan,
        "log2_ratio": log2_ratio,
        "p_value": pvalue,
    }


def _test_per_group_ranking(
    df: pd.DataFrame, config: PTMConstraintConfig
) -> pd.DataFrame:
    """Mann-Whitney per group, sorted by log2 depletion ratio."""
    rows: List[Dict[str, Any]] = []
    for group, block in df.groupby("primary_group", dropna=True):
        ptm = block.loc[block["ptm_any"], "af"].to_numpy()
        non = block.loc[~block["ptm_any"], "af"].to_numpy()

        if len(ptm) + len(non) < config.min_variants_per_group:
            continue

        row = {"group": str(group), **_mw_log2_ratio(ptm, non)}
        rows.append(row)

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("log2_ratio", ascending=False, na_position="last")


def _test_tau_stratification(
    df: pd.DataFrame, config: PTMConstraintConfig
) -> pd.DataFrame:
    """τ quartile bins × PTM depletion."""
    if df["tau"].isna().all():
        return pd.DataFrame()

    labels = ["Q1_housekeeping", "Q2", "Q3", "Q4_tissue_specific"]
    try:
        df = df.copy()
        df["tau_quartile"] = pd.qcut(
            df["tau"].to_numpy(), q=4, labels=labels, duplicates="drop"
        )
    except ValueError as exc:
        logger.warning("Could not bin τ into quartiles: %s", exc)
        return pd.DataFrame()

    rows: List[Dict[str, Any]] = []
    for q, block in df.groupby("tau_quartile", observed=True):
        ptm = block.loc[block["ptm_any"], "af"].to_numpy()
        non = block.loc[~block["ptm_any"], "af"].to_numpy()
        if len(ptm) + len(non) < config.min_variants_per_group:
            continue
        rows.append({"tau_quartile": str(q), **_mw_log2_ratio(ptm, non)})

    return pd.DataFrame(rows)


def _test_loeuf_group_factorial(
    df: pd.DataFrame, config: PTMConstraintConfig
) -> pd.DataFrame:
    """Two-way factorial: τ bin (tissue-specific vs broad) × LOEUF bin (constrained vs unconstrained)."""
    if df["loeuf"].isna().all() or df["tau"].isna().all():
        return pd.DataFrame()

    df = df.copy()
    df = df.dropna(subset=["loeuf", "tau"]).copy()
    if df.empty:
        return pd.DataFrame()

    try:
        df["tau_bin"] = pd.qcut(
            df["tau"].to_numpy(),
            q=2,
            labels=["broad", "tissue_specific"],
            duplicates="drop",
        )
        df["loeuf_bin"] = pd.qcut(
            df["loeuf"].to_numpy(),
            q=2,
            labels=["constrained", "unconstrained"],
            duplicates="drop",
        )
    except ValueError as exc:
        logger.warning("Could not bin factorial: %s", exc)
        return pd.DataFrame()

    rows: List[Dict[str, Any]] = []
    for (tau_bin, loeuf_bin), block in df.groupby(
        ["tau_bin", "loeuf_bin"], observed=True
    ):
        ptm = block.loc[block["ptm_any"], "af"].to_numpy()
        non = block.loc[~block["ptm_any"], "af"].to_numpy()
        if len(ptm) + len(non) < config.min_variants_per_group:
            continue
        rows.append(
            {
                "tau_bin": str(tau_bin),
                "loeuf_bin": str(loeuf_bin),
                **_mw_log2_ratio(ptm, non),
            }
        )

    return pd.DataFrame(rows)


def _normalise_category(raw: Any) -> List[str]:
    """Explode a set/list/string of PTM categories into lowercase tokens."""
    if raw is None:
        return []
    if isinstance(raw, (set, frozenset, list, tuple, np.ndarray)):
        cats = [str(x).lower() for x in raw if x is not None]
    else:
        cats = [str(raw).lower()]
    return cats or ["other"]


def _test_category_group_heatmap(
    df: pd.DataFrame, config: PTMConstraintConfig
) -> pd.DataFrame:
    """PTM category × group log2 ratios.

    PTM variants are exploded by category; the non-PTM comparator is the pool
    of non-PTM variants inside the same group.
    """
    rows: List[Dict[str, Any]] = []
    non_ptm = df[~df["ptm_any"]]

    for group, block in df.groupby("primary_group"):
        non_afs = non_ptm.loc[non_ptm["primary_group"] == group, "af"].to_numpy()
        if len(non_afs) < max(3, config.min_variants_per_group // 2):
            continue

        ptm_block = block[block["ptm_any"]].copy()
        if ptm_block.empty:
            continue

        ptm_block["_cats"] = ptm_block["ptm_categories"].apply(_normalise_category)
        exploded = ptm_block.explode("_cats")

        for category, cat_block in exploded.groupby("_cats"):
            if category is None or category == "":
                continue
            ptm_afs = cat_block["af"].to_numpy()
            if len(ptm_afs) < 3:
                continue
            rows.append(
                {
                    "group": str(group),
                    "category": str(category),
                    **_mw_log2_ratio(ptm_afs, non_afs),
                }
            )

    return pd.DataFrame(rows)


def _test_within_gene_paired(df: pd.DataFrame) -> Dict[str, Any]:
    """Per-gene paired comparison: median AF PTM vs non-PTM, Wilcoxon signed-rank."""
    rows: List[Dict[str, Any]] = []
    for gene, block in df.groupby("gene_key"):
        ptm = block.loc[block["ptm_any"], "af"].to_numpy()
        non = block.loc[~block["ptm_any"], "af"].to_numpy()
        if len(ptm) == 0 or len(non) == 0:
            continue
        med_ptm = float(np.median(ptm))
        med_non = float(np.median(non))
        rows.append(
            {
                "gene": str(gene),
                "n_ptm": int(len(ptm)),
                "n_non": int(len(non)),
                "med_ptm": med_ptm,
                "med_non": med_non,
                "delta": med_non - med_ptm,
            }
        )

    per_gene = pd.DataFrame(rows)

    summary: Dict[str, Any] = {
        "n_genes": len(per_gene),
        "p_value": np.nan,
        "statistic": np.nan,
        "pct_genes_depleted": np.nan,
        "median_delta": np.nan,
    }

    if len(per_gene) >= 3:
        deltas = per_gene["delta"].to_numpy()
        depleted = float((deltas > 0).mean())
        try:
            if np.all(deltas == 0):
                stat, pvalue = np.nan, np.nan
            else:
                result = stats.wilcoxon(deltas, alternative="greater")
                stat, pvalue = float(result.statistic), float(result.pvalue)
        except ValueError:
            stat, pvalue = np.nan, np.nan
        summary.update(
            {
                "n_genes": int(len(per_gene)),
                "p_value": pvalue,
                "statistic": stat,
                "pct_genes_depleted": depleted,
                "median_delta": float(np.median(deltas)),
            }
        )

    return {"per_gene": per_gene, "summary": summary}


def _write_tables(tables: Dict[str, pd.DataFrame], output_dir: str) -> None:
    for name, df in tables.items():
        path = os.path.join(output_dir, f"{name}.tsv")
        if df is None or df.empty:
            logger.warning("Test '%s' produced no rows; writing empty TSV.", name)
            pd.DataFrame().to_csv(path, sep="\t", index=False)
            continue
        df.to_csv(path, sep="\t", index=False)
        logger.info("Wrote %s (%d rows)", path, len(df))
