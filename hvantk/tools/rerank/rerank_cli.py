# hvantk/tools/rerank/rerank_cli.py
import click, yaml, pandas as pd
from hvantk.core.config import CONTEXT_SETTINGS
from hvantk.algorithms.rerank import rerank, Config, PriorSpec, CohortSpec, LabelSpec
from hvantk.algorithms.rerank.catalog.builders import table_axis, genelist_labels


@click.command(name="rerank", context_settings=CONTEXT_SETTINGS,
               help="Re-rank genes by multi-omic credibility from a declarative YAML config.")
@click.option("-c", "--config", "config_path", required=True, type=click.Path(exists=True),
              help="YAML config: prior, features (gene-keyed tables), labels, optional cohort.")
@click.option("-o", "--output", required=True, type=click.Path(), help="Output ranked TSV.")
def rerank_cmd(config_path, output):
    with open(config_path) as fh:
        spec = yaml.safe_load(fh)
    prior = PriorSpec(spec["prior"]["path"], spec["prior"]["unit_col"], spec["prior"]["stat_col"])
    feats = [table_axis(f["name"], f["path"]) for f in spec["features"]]
    labels = genelist_labels(spec["labels"]["path"])
    cohort = None
    if spec.get("cohort"):
        cohort = CohortSpec(spec["cohort"]["path"],
                            params={"veto_cols": spec["cohort"].get("veto_cols", []), "key": "gene"})
    from hvantk.algorithms.rerank.veto import CaseControlArchitectureVeto, NoOpVeto
    cfg = Config(name=spec.get("name", "rerank"), prior=prior, features=feats, labels=labels,
                 cohort=cohort, veto=CaseControlArchitectureVeto() if cohort else NoOpVeto(),
                 min_label_coverage=float(spec.get("min_label_coverage", 0.5)))
    res = rerank(cfg)
    res.table.to_csv(output, sep="\t", index=False)
    click.echo(f"Wrote {len(res.table)} genes -> {output}")
    click.echo("Per-axis ablation (delta-AUC over constraint):")
    click.echo(res.metrics.ablation.to_string(index=False))
