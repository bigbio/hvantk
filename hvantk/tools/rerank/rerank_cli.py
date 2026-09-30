# hvantk/tools/rerank/rerank_cli.py
import logging

import click
import jsonschema
import yaml

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)

# The axis label a cohort author would name an architecture axis after -- used only to
# detect a near-miss (declares this axis but not all three columns
# hvantk.algorithms.rerank.audit.ARCHITECTURE_AUDIT_COLUMNS needs) and warn about it.
# Any other axis combination that happens to supply the three columns is unaffected.
_ARCHITECTURE_AXIS_NAME = "architecture"

# Top-level keys this CLI understands. A leftover 'prior:' block from a config that
# predates the CohortManifest migration -- or any other stray key -- is rejected
# loudly rather than silently ignored: unlike the pre-migration PriorSpec-only config,
# the prior now lives entirely inside the cohort manifest this CLI loads via
# 'cohort:', so a 'prior:' block here can never be honoured and silently keeping it
# around invites exactly the confusion this check exists to prevent.
_KNOWN_TOP_LEVEL_KEYS = {"name", "cohort", "features", "labels", "min_label_coverage"}


@click.command(
    name="rerank",
    context_settings=CONTEXT_SETTINGS,
    help="Re-rank genes by multi-omic credibility from a declarative YAML config.",
)
@click.option(
    "-c",
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True),
    help="YAML config: cohort (a path to a cohort manifest YAML), features "
    "(gene-keyed tables), labels.",
)
@click.option(
    "-o", "--output", required=True, type=click.Path(), help="Output ranked TSV."
)
@click.option(
    "--seed",
    type=int,
    default=None,
    help="Random seed for the CV partition, the estimator, the bootstrap and the "
    "permutation null. Default: 42, the historic value.",
)
@click.option(
    "--seed-sweep",
    type=int,
    default=1,
    show_default=True,
    help="Recompute each axis's interval over this many CV seeds (seed, seed+1, ...) and "
    "report the envelope (the union of the gene-bootstrap interval and the across-seed "
    "range; never narrower) beside the gene-resampling interval. 1 leaves the reported "
    "interval blind to which genes landed in which fold.",
)
@click.option(
    "--blocks",
    "blocks_path",
    type=click.Path(exists=True),
    default=None,
    help="HGNC complete-set TSV (symbol/gene_group/status). Switches cross-validation to "
    "paralogue-blocked folds so a gene family cannot straddle a fold.",
)
@click.option(
    "--max-block-frac",
    type=float,
    default=None,
    help="Abort if the largest paralogue block exceeds this fraction of the universe "
    "[default: 0.1]. Only meaningful with --blocks.",
)
@click.option(
    "--n-perm",
    type=int,
    default=0,
    show_default=True,
    help="Label permutations for the per-axis and selected-maximum nulls. 0 disables the "
    "multiplicity correction entirely (the historic behaviour). Each permutation refits "
    "the baseline and every candidate axis (about n_perm x (1 + axes) five-fold fits; a "
    "selection policy, RFECV especially, multiplies that further) -- expect tens of "
    "minutes for a few hundred permutations on ~1,000 genes. Run `hvantk -v rerank ...` "
    "to see progress.",
)
@click.option(
    "--null-out",
    type=click.Path(),
    default=None,
    help="Write the per-axis null summary (TSV) here (the control setting the p-values "
    "belong to is printed to the console, not written to this file). Requires --n-perm. "
    "The API also supports chunked nulls (NullConfig(chunk=, n_chunks=) + "
    "NullDistribution.merge) for cluster array jobs; there is no --chunk/--n-chunks here.",
)
def rerank_cmd(config_path, output, seed, seed_sweep, blocks_path, max_block_frac,
               n_perm, null_out):
    # Heavy ML imports are deferred to invocation time so that importing the hvantk CLI
    # (and every other subcommand) does NOT require scikit-learn, which is an OPTIONAL
    # dependency. Mirrors the psroc/ancestry deferral pattern.
    from hvantk.algorithms.rerank import rerank, Config
    from hvantk.algorithms.rerank.blocks import BlockPolicy, DEFAULT_MAX_BLOCK_FRAC
    from hvantk.algorithms.rerank.catalog.builders import table_axis, genelist_labels
    from hvantk.algorithms.rerank.audit import (
        ARCHITECTURE_AUDIT_COLUMNS,
        CaseControlArchitectureAudit,
        NoAudit,
        has_architecture_columns,
    )
    from hvantk.algorithms.rerank.nulls import NullConfig
    from hvantk.algorithms.rerank.seeds import DEFAULT_SEED
    from hvantk.algorithms.cohort.spec import load_cohort

    # Flag-combination errors are checked before any file is read or any work starts --
    # a bad combination here is a typo, not a run worth paying for.
    if max_block_frac is not None and blocks_path is None:
        raise click.ClickException(
            "--max-block-frac only applies to paralogue-blocked folds; pass --blocks too "
            "(or drop it)."
        )
    if n_perm < 0:
        raise click.ClickException(f"--n-perm must be >= 0; got {n_perm}.")
    if n_perm == 0 and null_out:
        raise click.ClickException(
            "--null-out describes a permutation null; pass --n-perm N too."
        )

    with open(config_path) as fh:
        spec = yaml.safe_load(fh)
    if not isinstance(spec, dict):
        raise click.ClickException(
            f"config {config_path}: expected a YAML mapping, got {type(spec).__name__}"
        )
    for key in ("cohort", "features", "labels"):
        if key not in spec:
            raise click.ClickException(
                f"config {config_path}: missing required key '{key}'"
            )
    unknown = sorted(set(spec) - _KNOWN_TOP_LEVEL_KEYS)
    if unknown:
        prior_hint = (
            " 'prior:' is a pre-migration key -- the prior now lives in the cohort "
            "manifest's own 'prior:' block, referenced here via 'cohort:'; delete it."
            if "prior" in unknown
            else ""
        )
        raise click.ClickException(
            f"config {config_path}: unknown key(s) {unknown}."
            + prior_hint
            + f" Known keys: {sorted(_KNOWN_TOP_LEVEL_KEYS)}."
        )

    cohort_path = spec["cohort"]
    if not isinstance(cohort_path, str):
        raise click.ClickException(
            f"config {config_path}: 'cohort' must be a path to a cohort manifest YAML "
            f"(e.g. cohort: path/to/cohort.yaml), got {type(cohort_path).__name__}"
        )
    try:
        cohort = load_cohort(cohort_path)
    except OSError as exc:
        raise click.ClickException(
            f"config {config_path}: cohort manifest '{cohort_path}' could not be read: {exc}"
        )
    except (jsonschema.ValidationError, ValueError, yaml.YAMLError) as exc:
        raise click.ClickException(
            f"config {config_path}: invalid cohort manifest '{cohort_path}': {exc}"
        )

    try:
        feats = [table_axis(f["name"], f["path"]) for f in spec["features"]]
        labels = genelist_labels(spec["labels"]["path"])
    except (KeyError, TypeError) as exc:
        raise click.ClickException(
            f"config {config_path}: malformed features/labels block ({exc})"
        )

    # axis_columns(), not declared_columns(): engine.rerank() merges the cohort frame
    # with include_prior=False, so eligibility must be judged on exactly the columns
    # that merge actually contributes (see the matching comment in
    # hvantk.algorithms.rerank.catalog.registry._default_audit, which this CLI path
    # must never disagree with).
    has_architecture = has_architecture_columns(cohort.axis_columns())
    audit = CaseControlArchitectureAudit() if has_architecture else NoAudit()

    # Only worth warning about a near-miss (an "architecture"-named axis that doesn't
    # cover all three required columns) when the audit did NOT end up wired -- e.g. the
    # three columns are split across an "architecture" axis and a differently-named
    # axis (both contributing to declared_columns()), which has_architecture_columns()
    # already accounts for. Warning here regardless of has_architecture would fire on
    # working configs where CaseControlArchitectureAudit genuinely ran.
    if not has_architecture:
        for entry in cohort.cohort_axes:
            if entry.axis == _ARCHITECTURE_AXIS_NAME:
                missing = sorted(ARCHITECTURE_AUDIT_COLUMNS - set(entry.columns))
                if missing:
                    logger.warning(
                        "cohort %r: cohort_axes entry %r has columns %r, which is not "
                        "a superset of the columns CaseControlArchitectureAudit "
                        "requires %r -- missing %r; falling back to NoAudit unless "
                        "those columns are declared elsewhere in the manifest.",
                        cohort.name,
                        entry.axis,
                        sorted(entry.columns),
                        sorted(ARCHITECTURE_AUDIT_COLUMNS),
                        missing,
                    )
                break

    try:
        cfg = Config(
            name=spec.get("name", "rerank"),
            features=feats,
            labels=labels,
            cohort=cohort,
            audit=audit,
            min_label_coverage=float(spec.get("min_label_coverage", 0.5)),
            seed=DEFAULT_SEED if seed is None else seed,
            seed_sweep=seed_sweep,
            blocks=(
                None
                if blocks_path is None
                else BlockPolicy(
                    table=blocks_path,
                    max_block_frac=(
                        DEFAULT_MAX_BLOCK_FRAC if max_block_frac is None else max_block_frac
                    ),
                )
            ),
            nulls=(
                None
                if n_perm <= 0
                else NullConfig(n_perm=n_perm, seed=DEFAULT_SEED if seed is None else seed)
            ),
        )
    except (TypeError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc

    try:
        res = rerank(cfg)
    except ValueError as exc:
        # A hard abort, surfaced as a clean CLI failure rather than a traceback: the run
        # must not write a scored table whose pooled out-of-fold AUC estimates something
        # other than what the caller will compare it against. Every ValueError the engine
        # raises is a user-facing message (e.g. a gene-group table that matches no gene,
        # too few paralogue blocks for the fold count), so the CLI reports it cleanly
        # instead of a traceback.
        raise click.ClickException(str(exc)) from exc
    if n_perm > 0 and res.nulls is None:
        # cfg.nulls is a NullConfig here (n_perm > 0), so the only way _run_nulls returns
        # None is "no candidate axis besides the baseline" -- checked before any output is
        # written, the same as every other flag-combination error above.
        axis_names = [f.name for f in cfg.features]
        raise click.ClickException(
            "--n-perm needs at least one feature axis besides the baseline; this config "
            f"has only {axis_names}."
        )
    res.table.to_csv(output, sep="\t", index=False)
    click.echo(f"Wrote {len(res.table)} genes -> {output}")
    scored = res.table[res.table["score"].notna()]
    n_robust = int((scored["verdict"] == "ROBUST").sum())
    n_flagged = int(res.table["flag"].sum())
    click.echo(
        f"Credibility: ROBUST={n_robust} / INTERMEDIATE={len(scored) - n_robust} "
        f"(over {len(scored)} scored genes)"
    )
    click.echo(
        f"Audit: {n_flagged} gene(s) FLAGGED for review (advisory; ranking not overridden). "
        f"Reasons: {res.table.loc[res.table['flag'], 'flag_reason'].value_counts().to_dict()}"
    )
    click.echo("Per-axis ablation (delta-AUC over constraint):")
    click.echo(res.metrics.ablation.to_string(index=False))
    if res.blocks is not None:
        universe = len(res.blocks.blocks)
        click.echo(
            f"Blocks: {res.blocks.n_blocks} paralogue block(s); largest "
            f"{res.blocks.largest} ({res.blocks.largest_frac:.1%}); "
            f"{res.blocks.n_in_multi} gene(s) in a block of size > 1; "
            f"{res.blocks.n_matched:,}/{universe:,} genes matched the gene-group table; "
            f"digest {res.blocks.digest[:12]}"
        )
    if res.nulls is not None:
        # res.nulls.summary() (not the ablation's d_md -- a bootstrap median rounded to 3
        # dp, a different statistic from the one the null itself is made of): the null's
        # own scorer already computed `observed`, `p_per_axis` and `p_selected_max`.
        summary = res.nulls.summary()
        axis_word = "axis" if len(summary) == 1 else "axes"
        click.echo(
            f"Multiplicity: {len(summary)} {axis_word} searched over "
            f"{res.nulls.n_perm:,} permutations; selected-maximum null median "
            f"{summary['selmax_median'].iloc[0]:+.4f}"
        )
        click.echo(f"Setting: {res.nulls.setting.describe()}")
        click.echo(
            summary[["axis", "observed", "selmax_median", "p_selected_max"]].to_string(
                index=False
            )
        )
        if null_out:
            summary.to_csv(null_out, sep="\t", index=False)
            click.echo(f"Wrote the null summary -> {null_out}")
