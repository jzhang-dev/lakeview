#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Literal, Callable, Any, Optional
from collections.abc import Iterable, Sequence
import os
import pathlib
import tempfile
from numbers import Real
import matplotlib as mpl
from matplotlib.colors import rgb2hex
import pysam
from .custom_types import Color, NativeHashable


def get_cmap_colors(
    cmap_name: str, format_: Literal["hex", "rgb"] = "hex"
) -> list[Color]:
    """
    https://gist.github.com/jdbcode/33d37999f950a36b43e058d15280b536

    >>> get_cmap_colors("Set2")
    ['#66c2a5',
     '#fc8d62',
     '#8da0cb',
     '#e78ac3',
     '#a6d854',
     '#ffd92f',
     '#e5c494',
     '#b3b3b3']
    """
    cmap = mpl.colormaps.get_cmap(cmap_name, 256)
    colors = list({cmap(i)[:3]: None for i in range(cmap.N)})
    if format_ == "rgb":
        pass
    elif format_ == "hex":
        colors = [rgb2hex(c) for c in colors]
    else:
        raise ValueError("`format_` must be one of {'rgb', 'hex'}")
    return colors


def sort_by(
    *iterables: Sequence, by: Sequence[NativeHashable], reverse: bool = False
) -> list[list]:
    """
    Sort multiple equal-length lists by the value of another list.

    >>> a = ['a', 'b', 'c']
    >>> b = ['x', 'y', 'z']
    >>> r = [3, 1, 2]
    >>> sort_by(a, b, by=r, reverse=False)
    [('b', 'c', 'a'), ('y', 'z', 'x')]
    """
    # Check iterable lengths are equal
    lengths = [len(x) for x in iterables]
    if len(set(lengths)) > 1:
        raise ValueError("Iterables must have the same length.")
    # Return early if iterables are empty
    if any(l == 0 for l in lengths):
        return [[] for x in iterables]
    
    zipped_lists = list(zip(*iterables))
    sorted_zipped_lists = [
        x
        for (_, x) in sorted(
            zip(by, zipped_lists), key=lambda pair: pair[0], reverse=reverse
        )
    ]
    sorted_lists = list(zip(*sorted_zipped_lists))
    return sorted_lists


def filter_by(*iterables: Sequence, by: Sequence[bool]) -> list[list]:
    """
    Filter multiple equal-length lists by the value of another list.

    >>> a = ['a', 'b', 'c']
    >>> b = ['x', 'y', 'z']
    >>> f = [True, False, True]
    >>> filter_by(a, b, by=f)
    [('a', 'c'), ('x', 'z')]
    """
    # Check iterable lengths are equal
    lengths = [len(x) for x in iterables]
    if len(set(lengths)) > 1:
        raise ValueError("Iterables must have the same length.")
    # Return early if iterables are empty
    if any(l == 0 for l in lengths):
        return [[] for x in iterables]

    zipped_lists = list(zip(*iterables))
    filtered_zipped_lists = [x for (x, b) in zip(zipped_lists, by) if b]
    filtered_lists = list(zip(*filtered_zipped_lists))
    return filtered_lists


def draw_rigid_polygon(
    ax,
    shape: Iterable[tuple[Real, Real]],
    position: tuple[Real, Real],
    position_transform: mpl.transforms.Transform,
    **kw,
):
    """
    Draw a polygon patch with shape defined in inches and position defined in data, Axes or physical units.
    """
    if (
        len(position) != 2
        or not isinstance(position[0], Real)
        or not isinstance(position[1], Real)
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


def download_bam(
    bam_url: str,
    bai_url: str,
    region: str,
    output_bam_path: str,
    output_bai_path: Optional[str] = None,
    *,
    index: bool = True,
    override: bool = False,
):
    if not os.path.isfile(output_bam_path) or override:
        workdir = os.getcwd()
        with tempfile.TemporaryDirectory() as d:
            try:
                os.chdir(d)
                bam_data = pysam.view("-X", "-b", bam_url, bai_url, region)  # type: ignore
            finally:
                os.chdir(workdir)
        target_dir = os.path.split(output_bam_path)[0]
        pathlib.Path(target_dir).mkdir(parents=True, exist_ok=True)
        with open(output_bam_path, "wb") as f:
            f.write(bam_data)
    if output_bai_path is None:
        output_bai_path = output_bam_path + ".bai"
    if index and (not os.path.isfile(output_bai_path) or override):
        target_dir = os.path.split(output_bai_path)[0]
        pathlib.Path(target_dir).mkdir(parents=True, exist_ok=True)
        pysam.index(output_bam_path, output_bai_path)  # type: ignore


def pack_intervals(intervals: Iterable[tuple[float, float]]) -> list[int]:
    """
    Assign an non-negative offset to each input interval so that intervals sharing the same offset will not overlap with each other, while minimising offset values.
    Intervals are treated as being closed.

    >>> pack_intervals([(1, 2), (3, 4), (1, 3), (0, 5)])
    [0, 0, 1, 2]
    """
    occupied_intervals: list[list[tuple[float, float]]] = [[]]
    offsets = []
    for (start, end) in intervals:
        for offset, intervals in enumerate(occupied_intervals):
            for occupied_start, occupied_end in intervals:
                if (
                    occupied_start <= start <= occupied_end
                    or occupied_start <= end <= occupied_end
                    or start <= occupied_start <= end
                ):
                    break
            else:
                # Found place at the current offset.
                offsets.append(offset)
                occupied_intervals[offset].append((start, end))
                break
        else:
            # No places available.
            # Increment the offset.
            occupied_intervals.append([(start, end)])
            offsets.append(offset + 1)
    return offsets


def get_ax_size(ax):
    """
    Return the size of a given Axes in inches.
    """
    fig = ax.figure
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    return (width, height)


def scientific_notation(
    x: float, significant_figures: int = 3, *, quote: str = "$"
) -> str:
    """
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
