"""
Base visualization utilities for hvantk.

This module provides foundational functions for visualization settings,
styling, and common operations used across different visualization types.
"""

import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from typing import Optional, Union, Tuple, Dict, Any, List


def set_default_style(
    style: str = "default",
    context: str = "notebook",
    palette: str = "deep",
    font_scale: float = 1.0,
) -> None:
    """
    Set the default matplotlib style for visualizations.

    Parameters
    ----------
    style : str
        The style to use. Options include 'default', 'whitegrid', 'darkgrid', etc.
        If 'publication', a style suitable for publications is used.
    context : str
        The context setting. Options: 'paper', 'notebook', 'talk', 'poster'
    palette : str
        Color palette to use
    font_scale : float
        Scaling factor for font sizes

    Notes
    -----
    This function will attempt to use seaborn if available, otherwise falls back to
    matplotlib's built-in styles.
    """
    try:
        import seaborn as sns

        orig_style = style
        # Map styles to valid seaborn styles
        style_mapping = {"default": "whitegrid", "publication": "whitegrid"}
        seaborn_style = style_mapping.get(style, style)
        sns.set_theme(
            style=seaborn_style, context=context, palette=palette, font_scale=font_scale
        )
        if orig_style == "publication":
            # Publication-ready style settings
            plt.rcParams["font.family"] = "sans-serif"
            plt.rcParams["font.sans-serif"] = [
                "Arial",
                "DejaVu Sans",
                "Liberation Sans",
            ]
            plt.rcParams["axes.linewidth"] = 0.8
            plt.rcParams["axes.labelsize"] = 8
            plt.rcParams["xtick.labelsize"] = 7
            plt.rcParams["ytick.labelsize"] = 7
            plt.rcParams["legend.fontsize"] = 7
            plt.rcParams["figure.titlesize"] = 10
    except ImportError:
        # Fall back to matplotlib styles if seaborn is not available
        available_styles = plt.style.available
        if style in available_styles:
            plt.style.use(style)
        else:
            plt.style.use("default")

        # Basic font size adjustments
        plt.rcParams["font.size"] = 10 * font_scale
        plt.rcParams["axes.labelsize"] = 11 * font_scale
        plt.rcParams["axes.titlesize"] = 12 * font_scale
        plt.rcParams["xtick.labelsize"] = 9 * font_scale
        plt.rcParams["ytick.labelsize"] = 9 * font_scale


def save_figure(
    fig: plt.Figure,
    filename: str,
    dpi: int = 300,
    formats: Optional[List[str]] = None,
    transparent: bool = False,
    bbox_inches: str = "tight",
    output_dir: Optional[str] = None,
    **kwargs,
) -> None:
    """
    Save a matplotlib figure in multiple formats.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save
    filename : str
        Base name of the file to save (without extension)
    dpi : int
        Resolution in dots per inch
    formats : list of str, optional
        List of file formats to save (e.g., ['png', 'pdf', 'svg']).
        Defaults to ['png'] if not specified.
    transparent : bool
        Whether to save with a transparent background
    bbox_inches : str
        Bounding box in inches ('tight' typically works best)
    output_dir : str, optional
        Directory to save the figure. If None, saves in current directory.
    **kwargs
        Additional keyword arguments passed to plt.savefig()
    """
    # Initialize formats with default value if None
    if formats is None:
        formats = ["png"]

    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)

    for fmt in formats:
        # Handle case where filename already has an extension and convert Path to str
        filename_str = str(filename)
        base_name = (
            filename_str.replace(".png", "").replace(".pdf", "").replace(".svg", "")
        )
        save_path = (
            os.path.join(output_dir, f"{base_name}.{fmt}")
            if output_dir
            else f"{base_name}.{fmt}"
        )
        fig.savefig(
            save_path,
            format=fmt,
            dpi=dpi,
            transparent=transparent,
            bbox_inches=bbox_inches,
            **kwargs,
        )
        print(f"Figure saved: {save_path}")


def get_colors(
    n_colors: int, palette: str = "deep", as_cmap: bool = False
) -> Union[list, mpl.colors.Colormap]:
    """
    Generate a list of colors or a colormap.

    Parameters
    ----------
    n_colors : int
        Number of colors needed
    palette : str
        Name of the colormap/palette to use
    as_cmap : bool
        If True, returns a colormap object; otherwise returns a list of colors

    Returns
    -------
    Union[list, matplotlib.colors.Colormap]
        Either a list of colors as hex strings or a colormap object
    """
    try:
        import seaborn as sns

        if as_cmap:
            return sns.color_palette(palette, as_cmap=True)
        else:
            return sns.color_palette(palette, n_colors=n_colors).as_hex()
    except ImportError:
        # Fall back to matplotlib colormaps if seaborn is not available
        if as_cmap:
            return plt.get_cmap(palette)
        else:
            cmap = plt.get_cmap(palette)
            return [
                mpl.colors.rgb2hex(cmap(i / (n_colors - 1 if n_colors > 1 else 1)))
                for i in range(n_colors)
            ]


def add_figure_labels(
    fig: plt.Figure,
    labels: Dict[str, Tuple[float, float, str]],
    fontsize: int = 12,
    fontweight: str = "bold",
    **kwargs,
) -> None:
    """
    Add subplot labels (A, B, C, etc.) to figure.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to add labels to
    labels : dict
        Dictionary mapping label text to (x, y, alignment) coordinates
        Example: {'A': (0.05, 0.95, 'left top'), 'B': (0.5, 0.95, 'top center')}
        The order of alignment tokens is free-form (e.g., 'left top' or 'top left').
    fontsize : int
        Font size for labels
    fontweight : str
        Font weight for labels ('normal', 'bold', etc.)
    **kwargs
        Additional keyword arguments passed to fig.text()
    """
    valid_ha = {"left", "center", "right"}
    valid_va = {"top", "center", "bottom"}
    for label, (x, y, alignment) in labels.items():
        ha = va = None
        tokens = alignment.lower().split()
        for token in tokens:
            if ha is None and token in valid_ha:
                ha = token
                continue
            if va is None and token in valid_va:
                va = token
                continue
        # Sensible defaults
        if ha is None:
            ha = "center"
        if va is None:
            va = "center"
        # Validate
        if ha not in valid_ha:
            raise ValueError(
                f"Invalid horizontal alignment '{ha}' for label '{label}'. Must be one of {valid_ha}."
            )
        if va not in valid_va:
            raise ValueError(
                f"Invalid vertical alignment '{va}' for label '{label}'. Must be one of {valid_va}."
            )
        fig.text(
            x,
            y,
            label,
            ha=ha,
            va=va,
            fontsize=fontsize,
            fontweight=fontweight,
            **kwargs,
        )
