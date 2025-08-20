"""
Visualization module for hvantk.

This module provides functions and classes for visualizing multiomics data,
including heatmaps, volcano plots, PCA plots, and other common visualizations
used in multiomics analysis.
"""

from .base import (
    set_default_style,
    save_figure
)

from .omics import (
    plot_heatmap,
    plot_volcano,
    plot_pca,
    plot_umap,
    plot_expression_distribution,
    plot_sample_correlation,
    plot_variant_lollipop,
    plot_manhattan
)

from .clinical import (
    plot_survival_curve,
    plot_clinical_association
)

from .enrichment import (
    plot_enrichment_barplot,
    plot_enrichment_dotplot,
    plot_pathway_network
)

__all__ = [
    'set_default_style',
    'save_figure',
    'plot_heatmap',
    'plot_volcano',
    'plot_pca',
    'plot_umap',
    'plot_expression_distribution',
    'plot_sample_correlation',
    'plot_variant_lollipop',
    'plot_manhattan',
    'plot_survival_curve',
    'plot_clinical_association',
    'plot_enrichment_barplot',
    'plot_enrichment_dotplot',
    'plot_pathway_network'
]
