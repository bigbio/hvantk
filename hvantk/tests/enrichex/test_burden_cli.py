"""
Tests for the EnrichEx `burden` / `burden-pipeline` CLI commands.

These exercise pure pandas-level helpers and config construction only --
no Hail session is started.
"""

import pandas as pd
import pytest
from click.testing import CliRunner

import hvantk.algorithms.enrichex.pipeline as pipeline_module
from hvantk.algorithms.enrichex.burden import apply_correction_and_flag_significance
from hvantk.tools.enrichex.burden_cli import burden_pipeline_cmd


def test_apply_correction_and_flag_significance_handles_non_convergent_pvalue():
    """A non-convergent fit leaves p_value missing (Hail's to_pandas() gives
    a pandas NA there, per the comment in BurdenPipeline._run_single).
    Feeding that straight into apply_correction() raised TypeError
    ("boolean value of NA is ambiguous"); the CLI must instead skip that
    row in the correction and mark it not significant, like the pipeline
    path already does.
    """
    df = pd.DataFrame({"gene_set_name": ["A", "B", "C"], "p_value": [0.01, pd.NA, 0.2]})

    result = apply_correction_and_flag_significance(
        df, correction="benjamini-hochberg", alpha=0.05, n_total=3
    )

    assert result.loc[1, "p_value"] != result.loc[1, "p_value"]  # NaN
    assert result.loc[1, "p_adjusted"] != result.loc[1, "p_adjusted"]  # NaN
    assert bool(result.loc[1, "significant"]) is False
    # The valid rows are still corrected against n_total=3.
    assert result.loc[0, "p_adjusted"] == pytest.approx(0.03)
    assert result.loc[2, "p_adjusted"] == pytest.approx(0.3)


def test_burden_pipeline_cli_applies_filters_without_variant_classes(
    tmp_path, monkeypatch
):
    """burden-pipeline must honour --max-af/--pass-only/field options even
    when --variant-classes is not given. Previously, without
    --variant-classes, vc_dict stayed {} and the pipeline ran the single
    unstratified pass with variant_filter=None, silently dropping
    --max-af/--min-score/--pass-only/the field options.

    BurdenPipeline is stubbed so this never touches Hail: it only needs to
    inspect the BurdenConfig the real command builds.
    """
    cohort_mt = tmp_path / "cohort.mt"
    cohort_mt.write_text("")
    gene_sets_path = tmp_path / "genes.json"
    gene_sets_path.write_text("{}")
    output_dir = tmp_path / "out"

    captured = {}

    class _StubPipeline:
        def __init__(self, config):
            captured["config"] = config

        def show_plan(self):
            return ""

    monkeypatch.setattr(pipeline_module, "BurdenPipeline", _StubPipeline)

    runner = CliRunner()
    result = runner.invoke(
        burden_pipeline_cmd,
        [
            "-m",
            str(cohort_mt),
            "--gene-sets",
            f"heart:{gene_sets_path}",
            "--max-af",
            "0.001",
            "--pass-only",
            "--af-field",
            "internal_af",
            "-o",
            str(output_dir),
            "--dry-run",
        ],
    )

    assert result.exit_code == 0, result.output
    config = captured["config"]
    assert set(config.variant_classes) == {"all"}
    vf = config.variant_classes["all"]
    assert vf.max_af == 0.001
    assert vf.pass_only is True
    assert vf.af_field == "internal_af"


def test_burden_pipeline_dry_run_creates_no_directories(tmp_path):
    """--dry-run must not touch disk.

    BurdenPipeline.__init__ used to create output_dir/per_run/plots eagerly,
    and the CLI builds the pipeline before checking dry_run, so even a dry
    run left directories behind. Uses the real BurdenPipeline (unstubbed):
    show_plan() needs no Hail, so this stays a fast test.
    """
    cohort_mt = tmp_path / "cohort.mt"
    cohort_mt.write_text("")
    gene_sets_path = tmp_path / "genes.json"
    gene_sets_path.write_text("{}")
    output_dir = tmp_path / "out"

    runner = CliRunner()
    result = runner.invoke(
        burden_pipeline_cmd,
        [
            "-m",
            str(cohort_mt),
            "--gene-sets",
            f"heart:{gene_sets_path}",
            "-o",
            str(output_dir),
            "--dry-run",
        ],
    )

    assert result.exit_code == 0, result.output
    assert not output_dir.exists()
