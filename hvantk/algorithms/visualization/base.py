"""
Base visualization utilities for hvantk.

This module provides foundational functions for visualization settings,
styling, and common operations used across different visualization types.
"""

import base64
import io
import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib as mpl
from typing import Optional, Union, Tuple, Dict, List

logger = logging.getLogger(__name__)


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


_FORMAT_MIME_TYPES: dict[str, str] = {
    "png": "image/png",
    "jpg": "image/jpeg",
    "jpeg": "image/jpeg",
    "svg": "image/svg+xml",
    "pdf": "application/pdf",
    "tiff": "image/tiff",
    "tif": "image/tiff",
    "webp": "image/webp",
}


def encode_figure_to_base64(
    fig: plt.Figure,
    format: str = "png",
    dpi: int = 200,
    *,
    as_data_uri: bool = False,
) -> str:
    """Convert a matplotlib figure to a base64-encoded image string.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to encode.
    format : str
        Image format ('png', 'svg', etc.). Leading dots are stripped and the
        value is lowercased.
    dpi : int
        Resolution for raster formats.
    as_data_uri : bool, keyword-only
        When True, wrap the payload as ``data:image/{format};base64,{payload}``
        for direct embedding in an HTML ``<img src=...>`` attribute. When False
        (default), return the raw base64 payload.

    Returns
    -------
    str
        The base64 payload, optionally wrapped as a data URI.
    """
    format = format.lstrip(".").lower()
    buffer = io.BytesIO()
    fig.savefig(buffer, format=format, dpi=dpi, bbox_inches="tight")
    buffer.seek(0)
    payload = base64.b64encode(buffer.read()).decode("utf-8")
    buffer.close()
    if as_data_uri:
        mime = _FORMAT_MIME_TYPES.get(format, f"image/{format}")
        return f"data:{mime};base64,{payload}"
    return payload


def save_figure_to_path(
    fig: plt.Figure,
    output_path: Optional[str],
    format: Optional[str] = None,
    dpi: int = 300,
) -> None:
    """Save a figure to a single output path.

    No-op when ``output_path`` is falsy. When ``format`` is provided, the path's
    suffix is normalized to match it and the format is passed explicitly to
    :meth:`savefig`; otherwise the path's own suffix determines the format.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    output_path : str or None
        Destination path. If falsy, nothing is written.
    format : str, optional
        Explicit export format. If given, the path suffix is normalized to it.
    dpi : int
        Resolution in dots per inch.
    """
    if not output_path:
        return
    if format is not None:
        format = format.lstrip(".").lower()
    path = Path(output_path)
    if format is not None and path.suffix.lower() != f".{format}":
        path = path.with_suffix(f".{format}")
    path.parent.mkdir(parents=True, exist_ok=True)
    if format is not None:
        fig.savefig(path, dpi=dpi, bbox_inches="tight", format=format)
    else:
        fig.savefig(str(path), dpi=dpi, bbox_inches="tight")
    logger.info("Saved figure to %s", path)


def empty_figure(
    output_path: Optional[str] = None,
    format: str = "png",
    dpi: int = 300,
    title: str = "No Data",
    message: str = "No data available",
    figsize: Tuple[int, int] = (8, 4),
) -> plt.Figure:
    """Create a placeholder figure when there is no data to plot."""
    fig, ax = plt.subplots(figsize=figsize)
    ax.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        fontsize=14,
        color="#888888",
        transform=ax.transAxes,
    )
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    if output_path is not None:
        save_figure_to_path(fig, output_path, format=format, dpi=dpi)
    return fig
