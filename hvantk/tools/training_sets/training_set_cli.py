"""CLI for ML training datasets.

    hvantk training-set build --spec S.yaml --input clinvar=clinvar.ht ... \\
        --gene-matrix matrix.ht --gene-list is_chd_gene=chd.txt --output out/name

The gene matrix is an ``hvantk annotate compose`` output built from the same spec. The
command does its work in-process; scheduling is external (sbatch wrappers submit it).
"""

from __future__ import annotations

import logging

import click

from hvantk.tools._cli_utils import parse_named_paths

logger = logging.getLogger(__name__)


@click.group("training-set")
def training_set_group():
    """Build ML training datasets (labelled, tiered feature tables)."""


def _check_names(option: str, expected, given) -> None:
    """Raise a ``click.UsageError`` when ``given`` names differ from the spec's."""
    missing = sorted(set(expected) - set(given))
    unknown = sorted(set(given) - set(expected))
    if missing or unknown:
        parts = []
        if missing:
            parts.append("missing for: " + ", ".join(missing))
        if unknown:
            parts.append("not declared in the spec: " + ", ".join(unknown))
        raise click.UsageError(f"{option} does not match the spec ({'; '.join(parts)})")


def _read_gene_list(path: str) -> list[str]:
    from pathlib import Path

    symbols = Path(path).read_text().split()
    if not symbols:
        raise click.UsageError(f"gene list {path!r} is empty")
    return symbols


def _read_legacy_names(path: str) -> dict[str, str]:
    import csv

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or not {"column", "legacy_name"} <= set(
            reader.fieldnames
        ):
            raise click.UsageError(
                f"--legacy-names {path!r} needs a header with 'column' and 'legacy_name'"
            )
        return {row["column"]: row["legacy_name"] for row in reader}


def _read_gene_axes(path: str) -> dict[str, list[str]]:
    import json
    from pathlib import Path

    if not Path(path).exists():
        raise click.UsageError(
            f"gene-matrix manifest {path!r} not found; pass --gene-matrix-manifest "
            "(annotate compose writes <output>.manifest.json by default)"
        )
    axes = json.loads(Path(path).read_text()).get("axes")
    if not isinstance(axes, dict):
        raise click.UsageError(
            f"{path!r} has no 'axes' block; is it a compose manifest?"
        )
    return {axis: list(cols) for axis, cols in axes.items()}


@training_set_group.command("build")
@click.option(
    "--spec",
    "spec_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Feature-spec YAML with label, variant, gene and layer1 sections.",
)
@click.option(
    "--input",
    "input_opts",
    required=True,
    multiple=True,
    help="Repeatable NAME=PATH (.ht), one per input the spec names: the label "
    "source, each variant source and the gene section's HGNC table.",
)
@click.option(
    "--gene-matrix",
    "gene_matrix_path",
    required=True,
    help="The annotate compose output (.ht) built from the same spec.",
)
@click.option(
    "--gene-matrix-manifest",
    default=None,
    help="Its compose manifest (default: <gene-matrix>.manifest.json).",
)
@click.option(
    "--gene-list",
    "gene_list_opts",
    multiple=True,
    help="Repeatable COLUMN=FILE, one per gene_lists column: whitespace-separated "
    "gene symbols.",
)
@click.option(
    "--legacy-names",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Optional TSV (header: column, legacy_name) that fills the dictionary's "
    "legacy_name column.",
)
@click.option(
    "--output",
    required=True,
    help="Output prefix: writes <prefix>.ht, .tsv.gz, .dictionary.tsv and "
    ".manifest.json.",
)
def build_cmd(
    spec_path,
    input_opts,
    gene_matrix_path,
    gene_matrix_manifest,
    gene_list_opts,
    legacy_names,
    output,
):
    """Build a labelled, tiered ML training dataset from a feature spec."""
    import hashlib
    from pathlib import Path

    from hvantk.algorithms.annotation.spec import load_spec

    spec = load_spec(spec_path)
    if spec.label is None or spec.gene is None:
        raise click.UsageError(
            f"spec {spec.name!r} needs a label: and a gene: section for training-set build"
        )
    input_paths = parse_named_paths(input_opts, "--input")
    _check_names("--input", spec.input_names(), input_paths)
    gene_list_paths = parse_named_paths(gene_list_opts, "--gene-list")
    _check_names("--gene-list", [g.column for g in spec.gene_lists], gene_list_paths)
    gene_lists = {col: _read_gene_list(p) for col, p in gene_list_paths.items()}
    gene_axes = _read_gene_axes(
        gene_matrix_manifest or f"{gene_matrix_path}.manifest.json"
    )
    legacy = _read_legacy_names(legacy_names) if legacy_names else None
    spec_sha256 = hashlib.sha256(Path(spec_path).read_bytes()).hexdigest()

    import hail as hl

    from hvantk.core.io import read_manifest
    from hvantk.core.plugin.loader import get_registry
    from hvantk.core.utils.hail_context import init_hail
    from hvantk.skills.hgnc.streamers import HGNCGeneCatalogStreamer
    from hvantk.tools.training_sets.build import build_training_set, write_training_set

    # Declared training provenance of every source whose columns are scores: the variant
    # sources and the layer1 axes that aggregate per-variant scores (key: variant).
    registry = get_registry()
    scored = [v.source for v in spec.variant]
    scored += [e.source for e in spec.layer1 if e.aggregate is not None]
    score_provenance = {}
    for source in dict.fromkeys(scored):
        try:
            score_provenance[source] = registry.get_manifest(source).score_provenance
        except KeyError:
            logger.warning(
                "source %r is not a registered dataset; its columns need a tier in "
                "the spec",
                source,
            )
            score_provenance[source] = {}
    input_provenance = {
        name: read_manifest(Path(path)) for name, path in input_paths.items()
    }

    init_hail()
    result = build_training_set(
        spec,
        {name: hl.read_table(path) for name, path in input_paths.items()},
        hl.read_table(gene_matrix_path),
        gene_axes,
        gene_lists,
        HGNCGeneCatalogStreamer.from_path(input_paths[spec.gene.hgnc_input]),
        score_provenance,
        input_provenance,
        checkpoint=f"{output}.work.ht",
        spec_sha256=spec_sha256,
        legacy_names=legacy,
    )
    result.manifest["inputs"] = dict(input_paths)
    result.manifest["gene_matrix"] = gene_matrix_path
    result.manifest["gene_list_files"] = dict(gene_list_paths)
    paths = write_training_set(result, output)

    counts = result.manifest["counts"]
    tiers = result.manifest["tiers"]
    click.echo(
        f"training-set {spec.name}: {counts['rows']} rows x "
        f"{len(result.dictionary)} columns; label {counts['label']}; "
        f"tiers {tiers}; excluded {len(result.manifest['excluded'])}"
    )
    for kind, path in paths.items():
        click.echo(f"  {kind}: {path}")
