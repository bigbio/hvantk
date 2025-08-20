"""
Visualization module for hvantk.

This module provides functions and classes for visualizing multiomics data.
"""

from .base import (
    set_default_style,
    save_figure,
    get_colors,
    add_figure_labels,
)

# Lazy facade to avoid importing heavy backends at module import time

def visualize_expression_distribution(*args, **kwargs):
    from .expression.hail import visualize_expression_distribution as _impl
    return _impl(*args, **kwargs)

__all__ = [
    'set_default_style',
    'save_figure',
    'get_colors',
    'add_figure_labels',
    'visualize_expression_distribution',
]
