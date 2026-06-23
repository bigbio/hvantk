# hvantk/tests/rerank/test_cli.py
import numpy as np, pandas as pd, yaml
from click.testing import CliRunner
from hvantk.tools.rerank import rerank_cmd


def test_rerank_cli_end_to_end(tmp_path):
    rng = np.random.default_rng(0); n = 150; genes = [f"g{i}" for i in range(n)]
    y = (rng.random(n) < 0.3).astype(int)
    pd.DataFrame({"gene": genes, "p": rng.random(n)}).to_csv(tmp_path/"prior.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": genes, "x": y + rng.normal(0, .5, n)}).to_parquet(tmp_path/"feat.parquet")
    (tmp_path/"labels.txt").write_text("\n".join(g for g, yy in zip(genes, y) if yy == 1))
    cfg = {"name": "toy",
           "prior": {"path": str(tmp_path/"prior.tsv"), "unit_col": "gene", "stat_col": "p"},
           "features": [{"name": "constraint", "path": str(tmp_path/"feat.parquet")}],
           "labels": {"path": str(tmp_path/"labels.txt")}, "min_label_coverage": 0.0}
    (tmp_path/"c.yaml").write_text(yaml.safe_dump(cfg))
    out = tmp_path/"out.tsv"
    r = CliRunner().invoke(rerank_cmd, ["-c", str(tmp_path/"c.yaml"), "-o", str(out)])
    assert r.exit_code == 0, r.output
    t = pd.read_csv(out, sep="\t"); assert len(t) == n and {"gene", "score", "verdict"} <= set(t.columns)
