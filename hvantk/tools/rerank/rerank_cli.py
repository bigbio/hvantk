# hvantk/tools/rerank/rerank_cli.py
import logging

import click
import jsonschema
import yaml

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)

# The columns CaseControlArchitectureAudit needs (audit.py). When the cohort manifest
# declares all three -- via any cohort_axes entry -- the audit is wired up automatically;
# otherwise the run falls back to NoAudit.
_ARCHITECTURE_AUDIT_COLUMNS = {"n_case_var", "conc", "driver_af"}
# The axis label a cohort author would name an architecture axis after -- used only to
# detect a near-miss (declares this axis but not all three columns above) and warn about
# it. Any other axis combination that happens to supply the three columns is unaffected.
_ARCHITECTURE_AXIS_NAME = "architecture"


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
def rerank_cmd(config_path, output):
    # Heavy ML imports are deferred to invocation time so that importing the hvantk CLI
    # (and every other subcommand) does NOT require scikit-learn, which is an OPTIONAL
    # dependency. Mirrors the psroc/ancestry deferral pattern.
    from hvantk.algorithms.rerank import rerank, Config
    from hvantk.algorithms.rerank.catalog.builders import table_axis, genelist_labels
    from hvantk.algorithms.rerank.audit import CaseControlArchitectureAudit, NoAudit
    from hvantk.algorithms.cohort.spec import load_cohort

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

    for entry in cohort.cohort_axes:
        if entry.axis == _ARCHITECTURE_AXIS_NAME:
            missing = sorted(_ARCHITECTURE_AUDIT_COLUMNS - set(entry.columns))
            if missing:
                logger.warning(
                    "cohort %r: cohort_axes entry %r has columns %r, which is not a "
                    "superset of the columns CaseControlArchitectureAudit requires "
                    "%r -- missing %r; falling back to NoAudit unless those columns "
                    "are declared elsewhere in the manifest.",
                    cohort.name,
                    entry.axis,
                    sorted(entry.columns),
                    sorted(_ARCHITECTURE_AUDIT_COLUMNS),
                    missing,
                )
            break

    has_architecture = _ARCHITECTURE_AUDIT_COLUMNS <= set(cohort.declared_columns())
    audit = CaseControlArchitectureAudit() if has_architecture else NoAudit()

    cfg = Config(
        name=spec.get("name", "rerank"),
        features=feats,
        labels=labels,
        cohort=cohort,
        audit=audit,
        min_label_coverage=float(spec.get("min_label_coverage", 0.5)),
    )
    res = rerank(cfg)
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
