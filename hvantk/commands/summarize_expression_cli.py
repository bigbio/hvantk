"""
CLI commands for expression MatrixTable analysis.

Grouped under ``hvantk expression``:

- ``hvantk expression describe``: inspect metadata fields in an expression MT
- ``hvantk expression summarize``: collapse an expression MT into a gene-level
  summary Table grouped by metadata fields
- ``hvantk expression markers``: extract marker genes from a summary Table
"""

import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


@click.group("expression", context_settings=CONTEXT_SETTINGS)
def expression_group():
    """Expression MatrixTable analysis commands."""
    pass


@expression_group.command("describe")
@click.option(
    "-m",
    "--matrix-table",
    type=click.Path(),
    required=True,
    help="Path to an expression MatrixTable (.mt).",
)
def describe_expression_cmd(matrix_table):
    """Inspect column metadata fields in an expression MatrixTable.

    Prints available grouping variables, their types, and level counts.
    Uses pre-computed column_summary when available (instant, no Spark job).

    \b
    Example:
      hvantk expression describe -m data/heart_sc.mt
    """
    from hvantk.utils.matrix_utils import describe_expression_mt

    describe_expression_mt(matrix_table)


@expression_group.command("summarize")
@click.option(
    "-m",
    "--matrix-table",
    type=click.Path(),
    required=True,
    help="Path to an expression MatrixTable (.mt).",
)
@click.option(
    "--group-by",
    multiple=True,
    required=True,
    help="Column metadata field(s) to group by. Repeat for multi-field grouping.",
)
@click.option(
    "--filter-by",
    multiple=True,
    default=None,
    help="Pre-filter columns: FIELD=VALUE (repeatable). "
    "Example: --filter-by time_point=9wpc --filter-by region=LV",
)
@click.option(
    "--expr-field",
    default="x",
    show_default=True,
    help="Entry field containing expression values.",
)
@click.option(
    "--gene-id-field",
    default="GeneID",
    show_default=True,
    help="Row field for gene IDs.",
)
@click.option(
    "--gene-name-field",
    default="Gene Name",
    show_default=True,
    help="Row field for gene names (use '' to omit).",
)
@click.option(
    "--min-cells",
    type=int,
    default=50,
    show_default=True,
    help="Minimum cells per group.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    required=True,
    help="Output path for the summary Hail Table (.ht).",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite existing output.",
)
def summarize_expression_cmd(
    matrix_table,
    group_by,
    filter_by,
    expr_field,
    gene_id_field,
    gene_name_field,
    min_cells,
    output,
    overwrite,
):
    """Collapse an expression MatrixTable into a gene-level summary Table.

    \b
    Examples:

      # Single-field grouping
      hvantk expression summarize \\
          -m data/heart_sc.mt \\
          --group-by cell_type \\
          -o data/heart_celltype_summary.ht

      # Multi-field grouping
      hvantk expression summarize \\
          -m data/heart_sc.mt \\
          --group-by cell_type --group-by region \\
          -o data/heart_celltype_region_summary.ht

      # With pre-filtering
      hvantk expression summarize \\
          -m data/heart_sc.mt \\
          --group-by cell_type \\
          --filter-by time_point=9wpc --filter-by region=LV \\
          --min-cells 50 \\
          -o data/heart_9wpc_LV_summary.ht
    """
    import hail as hl
    from hvantk.utils.matrix_utils import summarize_expression

    # Parse filter_by from "field=value" strings
    filters = None
    if filter_by:
        filters = {}
        for item in filter_by:
            if "=" not in item:
                raise click.BadParameter(
                    f"Expected FIELD=VALUE format, got: '{item}'",
                    param_hint="--filter-by",
                )
            key, value = item.split("=", 1)
            filters[key.strip()] = value.strip()

    gene_name = gene_name_field if gene_name_field else None

    mt = hl.read_matrix_table(matrix_table)

    tb = summarize_expression(
        mt,
        group_by=list(group_by),
        filter_by=filters,
        expr_field=expr_field,
        gene_id_field=gene_id_field,
        gene_name_field=gene_name,
        min_cells_per_group=min_cells,
        output_path=output,
        overwrite=overwrite,
    )

    n_genes = tb.count()
    # Collect group labels and sample sizes from one row
    sample_stats = hl.eval(tb.take(1)[0].stats) if n_genes > 0 else {}
    group_labels = sorted(sample_stats.keys())

    click.echo(f"\nSummary table written to: {output}")
    click.echo(f"  Genes:  {n_genes:,}")
    click.echo(f"  Groups: {len(group_labels)} (from {list(group_by)})")
    click.echo("")
    tb.describe()
    click.echo("")
    if group_labels:
        click.echo(f"Group labels ({len(group_labels)}):")
        # Show first few with cell counts
        for label in group_labels[:10]:
            info = sample_stats[label]
            click.echo(f"  {label:<40s}  {info.n_cells:>6,} cells")
        if len(group_labels) > 10:
            click.echo(f"  ... and {len(group_labels) - 10} more")


@expression_group.command("markers")
@click.option(
    "-s",
    "--summary",
    type=click.Path(),
    default=None,
    help="Path to an expression summary Hail Table (.ht) from 'expression summarize'. "
    "Required for fold_change/specificity methods; optional for wilcoxon (pre-filter).",
)
@click.option(
    "-m",
    "--matrix-table",
    type=click.Path(),
    default=None,
    help="Path to expression MatrixTable (.mt). Required for wilcoxon method.",
)
@click.option(
    "--method",
    type=click.Choice(["fold_change", "specificity", "wilcoxon"]),
    default="fold_change",
    show_default=True,
    help="Marker scoring method.",
)
@click.option(
    "--group-by",
    multiple=True,
    default=None,
    help="Metadata field(s) for grouping. Required for wilcoxon method.",
)
@click.option(
    "--filter-by",
    multiple=True,
    default=None,
    help="Pre-filter columns: FIELD=VALUE (repeatable). For wilcoxon method.",
)
@click.option(
    "--top-n",
    type=int,
    default=200,
    show_default=True,
    help="Maximum markers per group.",
)
@click.option(
    "--min-fold-change",
    type=float,
    default=1.5,
    show_default=True,
    help="Minimum fold change to qualify as marker.",
)
@click.option(
    "--min-fraction-expressed",
    type=float,
    default=0.1,
    show_default=True,
    help="Minimum fraction of cells expressing a gene in the group.",
)
@click.option(
    "--expr-field",
    default="x",
    show_default=True,
    help="Entry field containing expression values (wilcoxon).",
)
@click.option(
    "--gene-id-field",
    default="GeneID",
    show_default=True,
    help="Row field for gene IDs (wilcoxon).",
)
@click.option(
    "--gene-name-field",
    default="Gene Name",
    show_default=True,
    help="Row field for gene names (wilcoxon; use '' to omit).",
)
@click.option(
    "--min-cells",
    type=int,
    default=3,
    show_default=True,
    help="Minimum cells per group (wilcoxon).",
)
@click.option(
    "--max-candidates",
    type=int,
    default=2000,
    show_default=True,
    help="Maximum candidate genes to test after pre-filtering (wilcoxon).",
)
@click.option(
    "--correction",
    type=click.Choice(["benjamini-hochberg", "bonferroni", "none"]),
    default="benjamini-hochberg",
    show_default=True,
    help="Multiple testing correction method (wilcoxon).",
)
@click.option(
    "--alpha",
    type=float,
    default=0.05,
    show_default=True,
    help="Adjusted p-value threshold (wilcoxon).",
)
@click.option(
    "--results-tsv",
    type=click.Path(),
    default=None,
    help="Save full Wilcoxon results table to TSV.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    required=True,
    help="Output path (.json or .gmt).",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite existing output.",
)
def markers_cmd(
    summary,
    matrix_table,
    method,
    group_by,
    filter_by,
    top_n,
    min_fold_change,
    min_fraction_expressed,
    expr_field,
    gene_id_field,
    gene_name_field,
    min_cells,
    max_candidates,
    correction,
    alpha,
    results_tsv,
    output,
    overwrite,
):
    """Extract top marker genes per group from expression data.

    Three scoring methods are available, operating on different inputs:

    \b
    fold_change / specificity (ratio-based, from summary table):
      Requires --summary (-s) pointing to a pre-computed summary Hail Table
      produced by 'hvantk expression summarize'. Fast, no statistical test.

    \b
    wilcoxon (statistical test, from MatrixTable):
      Requires --matrix-table (-m) and --group-by. Runs a one-vs-rest
      Wilcoxon rank-sum test (Mann-Whitney U) with tie correction and
      Benjamini-Hochberg p-value adjustment. Optionally accepts --summary
      for Phase 1 candidate pre-filtering to reduce compute.

    \b
    Examples:

      # Fold-change markers from summary table
      hvantk expression markers \\
          -s data/heart_celltype_summary.ht \\
          --method fold_change \\
          --top-n 200 \\
          -o gene_sets/heart_cell_types.json

      # Wilcoxon rank-sum markers from MatrixTable
      hvantk expression markers \\
          -m data/heart_sc.mt \\
          --method wilcoxon \\
          --group-by cell_type \\
          --top-n 200 --alpha 0.05 \\
          --results-tsv results/wilcoxon_full.tsv \\
          -o gene_sets/heart_wilcoxon.json

      # Wilcoxon with summary pre-filter (faster on large datasets)
      hvantk expression markers \\
          -m data/heart_sc.mt \\
          -s data/heart_celltype_summary.ht \\
          --method wilcoxon \\
          --group-by cell_type \\
          -o gene_sets/heart_wilcoxon.json
    """
    from pathlib import Path

    # --- Validate method-specific requirements ---
    if method in ("fold_change", "specificity"):
        if summary is None:
            raise click.UsageError(f"--summary is required for method '{method}'.")
    elif method == "wilcoxon":
        if matrix_table is None:
            raise click.UsageError("--matrix-table is required for method 'wilcoxon'.")
        if not group_by:
            raise click.UsageError("--group-by is required for method 'wilcoxon'.")

    output_path = Path(output)
    if output_path.exists() and not overwrite:
        click.echo(
            f"Error: Output file already exists: {output_path}\n"
            "Use --overwrite to replace it.",
            err=True,
        )
        raise SystemExit(1)

    # --- Dispatch ---
    if method in ("fold_change", "specificity"):
        from hvantk.utils.gene_sets import extract_marker_gene_sets

        collection = extract_marker_gene_sets(
            summary=summary,
            n_markers=top_n,
            min_fold_change=min_fold_change,
            min_fraction_expressed=min_fraction_expressed,
            method=method,
        )

    else:  # wilcoxon
        import hail as hl
        from hvantk.utils.wilcoxon import WilcoxonParams
        from hvantk.utils.wilcoxon_hail import wilcoxon_markers_from_mt

        # Parse filter_by
        filters = None
        if filter_by:
            filters = {}
            for item in filter_by:
                if "=" not in item:
                    raise click.BadParameter(
                        f"Expected FIELD=VALUE format, got: '{item}'",
                        param_hint="--filter-by",
                    )
                key, value = item.split("=", 1)
                filters[key.strip()] = value.strip()

        gene_name = gene_name_field if gene_name_field else None

        params = WilcoxonParams(
            min_fold_change=min_fold_change,
            min_fraction_expressed=min_fraction_expressed,
            max_candidates=max_candidates,
            top_n=top_n,
            correction_method=correction,
            alpha=alpha,
        )

        mt = hl.read_matrix_table(matrix_table)
        results_df, collection = wilcoxon_markers_from_mt(
            mt,
            group_by=list(group_by),
            filter_by=filters,
            summary=summary,
            params=params,
            expr_field=expr_field,
            gene_id_field=gene_id_field,
            gene_name_field=gene_name,
            min_cells_per_group=min_cells,
        )

        # Save full results TSV if requested
        if results_tsv:
            tsv_path = Path(results_tsv)
            tsv_path.parent.mkdir(parents=True, exist_ok=True)
            results_df.to_csv(str(tsv_path), sep="\t", index=False)
            click.echo(f"Full results table: {tsv_path} ({len(results_df):,} rows)")

    if not collection.gene_sets:
        click.echo("Error: No marker gene sets produced.", err=True)
        raise SystemExit(1)

    # Report
    click.echo(f"Extracted markers for {len(collection)} groups (method={method})")
    for gs in sorted(collection, key=lambda g: -g.n_genes):
        click.echo(f"  {gs.name}: {gs.n_genes} markers")

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.suffix.lower() == ".gmt":
        collection.save_gmt(str(output_path))
    else:
        collection.save(str(output_path))

    click.echo(f"\nSaved to: {output_path}")
