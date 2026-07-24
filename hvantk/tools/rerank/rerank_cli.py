# hvantk/tools/rerank/rerank_cli.py
import click
import yaml

from hvantk.core.config import CONTEXT_SETTINGS


@click.command(name="rerank", context_settings=CONTEXT_SETTINGS,
               help="Re-rank genes by multi-omic credibility from a declarative YAML config.")
@click.option("-c", "--config", "config_path", required=True, type=click.Path(exists=True),
              help="YAML config: prior, features (gene-keyed tables), labels, optional cohort.")
@click.option("-o", "--output", required=True, type=click.Path(), help="Output ranked TSV.")
def rerank_cmd(config_path, output):
    # Heavy ML imports are deferred to invocation time so that importing the hvantk CLI
    # (and every other subcommand) does NOT require scikit-learn, which is an OPTIONAL
    # dependency. Mirrors the psroc/ancestry deferral pattern.
    from hvantk.algorithms.rerank import rerank, Config, PriorSpec, CohortSpec
    from hvantk.algorithms.rerank.catalog.builders import table_axis, genelist_labels
    from hvantk.algorithms.rerank.audit import CaseControlArchitectureAudit, NoAudit

    with open(config_path) as fh:
        spec = yaml.safe_load(fh)
    if not isinstance(spec, dict):
        raise click.ClickException(
            f"config {config_path}: expected a YAML mapping, got {type(spec).__name__}")
    for key in ("prior", "features", "labels"):
        if key not in spec:
            raise click.ClickException(f"config {config_path}: missing required key '{key}'")
    try:
        prior = PriorSpec(spec["prior"]["path"], spec["prior"]["unit_col"], spec["prior"]["stat_col"])
        feats = [table_axis(f["name"], f["path"]) for f in spec["features"]]
        labels = genelist_labels(spec["labels"]["path"])
    except (KeyError, TypeError) as exc:
        raise click.ClickException(
            f"config {config_path}: malformed prior/features/labels block ({exc})")

    cohort = None
    if spec.get("cohort"):
        cohort = CohortSpec(spec["cohort"]["path"],
                            params={"audit_cols": spec["cohort"].get("audit_cols", []), "key": "gene"})
    cfg = Config(name=spec.get("name", "rerank"), prior=prior, features=feats, labels=labels,
                 cohort=cohort, audit=CaseControlArchitectureAudit() if cohort else NoAudit(),
                 min_label_coverage=float(spec.get("min_label_coverage", 0.5)))
    res = rerank(cfg)
    res.table.to_csv(output, sep="\t", index=False)
    click.echo(f"Wrote {len(res.table)} genes -> {output}")
    scored = res.table[res.table["score"].notna()]
    n_robust = int((scored["verdict"] == "ROBUST").sum())
    n_flagged = int(res.table["flag"].sum())
    click.echo(f"Credibility: ROBUST={n_robust} / INTERMEDIATE={len(scored) - n_robust} "
               f"(over {len(scored)} scored genes)")
    click.echo(f"Audit: {n_flagged} gene(s) FLAGGED for review (advisory; ranking not overridden). "
               f"Reasons: {res.table.loc[res.table['flag'], 'flag_reason'].value_counts().to_dict()}")
    click.echo("Per-axis ablation (delta-AUC over constraint):")
    click.echo(res.metrics.ablation.to_string(index=False))
