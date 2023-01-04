#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Literal, Sequence
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex
from ._type_alias import Color, Point, Axes


def get_cmap_colors(
    cmap_name: str, format_: Literal["hex", "rgb"] = "hex"
) -> list[Color]:
    """
    https://gist.github.com/jdbcode/33d37999f950a36b43e058d15280b536

    >>> get_cmap_colors("Set2")
    ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3']
    """
    cmap = plt.get_cmap(cmap_name)
    colors = list({cmap(i)[:3]: None for i in range(cmap.N)})
    if format_ == "rgb":
        pass
    elif format_ == "hex":
        colors = [rgb2hex(c) for c in colors]
    else:
        raise ValueError("`format_` must be one of {'rgb', 'hex'}")
    return colors


def get_random_colors(n: int, *, seed=0, cmap="hsv") -> list[Color]:
    """
    Returns `n` random colors selected uniformly from `cmap`.

    >>> get_random_colors(3)
    [(0.0, 0.22463448382566054, 1.0, 1.0), (0.4018366371307548, 1.0, 0.0, 1.0), (1.0, 0.2316178786767022, 0.0, 1.0)]
    """
    rng = np.random.default_rng(seed=0)
    cmap = plt.get_cmap("hsv")
    colors = [cmap(rng.random()) for __ in range(n)]
    return colors


def draw_rigid_polygon(
    ax: Axes,
    shape: Sequence[Point],
    position: Point,
    position_transform: mpl.transforms.Transform,
    **kw,
):
    """
    Draw a polygon patch with shape defined in inches and position defined in data, Axes or physical units.
    """
    if (
        len(position) != 2
        or not isinstance(position[0], float)
        or not isinstance(position[1], float)
    ):
        raise ValueError(f"position={position}")
    transform = ax.figure.dpi_scale_trans + mpl.transforms.ScaledTranslation(
        *position, position_transform
    )
    polygon = mpl.patches.Polygon(
        shape,
        transform=transform,
        **kw,
    )
    ax.add_patch(polygon)


def get_ax_size(ax: Axes):
    """
    Return the size of a given Axes in inches.

    >>> fig, ax = plt.subplots(figsize=(8, 6))
    >>> get_ax_size(ax)
    (6.2, 4.62)
    """
    fig = ax.figure
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    return (width, height)


def scientific_notation(
    x: float, significant_figures: int = 3, *, quote: str = "$"
) -> str:
    r"""
    Returns scientific notation in Matplotlib mathtext format.

    >>> scientific_notation(0.000000013923, 4)
    '$1.392\\times 10^{-8}$'
    """
    if significant_figures < 1:
        raise ValueError()
    if not isinstance(significant_figures, int):
        raise ValueError()
    float_digits = significant_figures - 1
    # Python native scientific notation
    coefficient, exponent = format(x, f".{float_digits}e").split("e")
    # Remove leading zeros
    exponent = str(int(exponent))
    # Format with mathtext
    s = quote + coefficient + r"\times 10^" + "{" + exponent + "}" + quote
    return s



class BasePairFormatter(mpl.ticker.FuncFormatter):
    """
    A Matplotlib `tick formatter <https://matplotlib.org/stable/api/ticker_api.html#tick-formatting>`_ for common base pair units (bp, kb, Mb, Gb, Tb).

    >>> formatter = BasePairFormatter('kb')
    >>> formatter(24310)
    '24.310 kb'
    >>> formatter = BasePairFormatter('Mb', 2, show_suffix=False)
    >>> formatter(844_293_192)
    '844.29'
    """

    def __init__(
        self,
        unit: Literal["bp", "kb", "Mb", "Gb", "Tb"],
        decimals: int | None = None,
        *,
        show_suffix: bool = True,
    ):
        """
        :param unit: Base pair unit.
        :param decimals: The number of decimal places to show for the coefficient. The default is 1 for 'bp', or 3 for other values of `unit`.
        :param show_suffix: Whether to show the unit as a suffix.
        """
        unit_divisor = self._get_unit_divisor(unit)
        
        if decimals is None:
            if unit == 'bp':
                decimals = 1
            else:
                decimals = 3

        def formatter_function(x: float, pos) -> str:
            tick_label = format(x / unit_divisor, f".{decimals}f")
            #print(unit_divisor)
            if show_suffix:
                tick_label += " " + unit
            return tick_label

        super().__init__(formatter_function)

    @staticmethod
    def _get_unit_divisor(unit: Literal["bp", "kb", "Mb", "Gb", "Tb"]) -> int:
        UNIT_DIVISOR_DICT: dict[str, int] = dict(
            bp=1, kb=int(1e3), Mb=int(1e6), Gb=int(1e9), Tb=int(1e12)
        )
        if unit in UNIT_DIVISOR_DICT:
            unit_divisor: int = UNIT_DIVISOR_DICT[unit]
            return unit_divisor
        else:
            raise ValueError(
                f"Invalid value for `unit`: {unit!r}. Supported values: {tuple(UNIT_DIVISOR_DICT)!r}."
            )
