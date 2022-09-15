from typing import Tuple, Iterable
from numbers import Real
import matplotlib as mpl
from matplotlib.colors import rgb2hex
from matplotlib import cm


def get_cmap_colors(cmap_name, format_="hex"):
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
    cmap = cm.get_cmap(cmap_name, 256)
    colors = list({cmap(i)[:3]: None for i in range(cmap.N)})
    if format_ == "rgb":
        pass
    elif format_ == "hex":
        colors = [rgb2hex(c) for c in colors]
    else:
        raise ValueError("`format_` must be one of {'rgb', 'hex'}")
    return colors


def sort_by(*iterables, by, reverse=False):
    """
    Sort multiple equal-length lists by the value of another list.

    >>> a = ['a', 'b', 'c']
    >>> b = ['x', 'y', 'z']
    >>> r = [3, 1, 2]
    >>> sort_by(a, b, by=r, reverse=False)
    [('b', 'c', 'a'), ('y', 'z', 'x')]
    """
    zipped_lists = list(zip(*iterables))
    sorted_zipped_lists = [
        x
        for (_, x) in sorted(
            zip(by, zipped_lists), key=lambda pair: pair[0], reverse=reverse
        )
    ]
    sorted_lists = list(zip(*sorted_zipped_lists))
    return sorted_lists


def filter_by(*iterables, by):
    """
    Filter multiple equal-length lists by the value of another list.

    >>> a = ['a', 'b', 'c']
    >>> b = ['x', 'y', 'z']
    >>> f = [True, False, True]
    >>> filter_by(a, b, by=f)
    [('a', 'c'), ('x', 'z')]
    """
    zipped_lists = list(zip(*iterables))
    filtered_zipped_lists = [x for (x, b) in zip(zipped_lists, by) if b]
    filtered_lists = list(zip(*filtered_zipped_lists))
    return filtered_lists


# def bind(instance, func, as_name=None):
#     """
#     Bind the function *func* to *instance*, with either provided name *as_name*
#     or the existing name of *func*. The provided *func* should accept the
#     instance as the first argument, i.e. "self".
#     https://stackoverflow.com/questions/1015307/python-bind-an-unbound-method/1015405#1015405
#     """
#     if as_name is None:
#         as_name = func.__name__
#     bound_method = func.__get__(instance, instance.__class__)
#     setattr(instance, as_name, bound_method)
#     return bound_method


def draw_rigid_polygon(
    ax,
    shape: Iterable[Tuple[Real, Real]],
    position: Tuple[Real, Real],
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
    polygon = mpl.patches.Polygon(shape, transform=transform, **kw,)
    ax.add_patch(polygon)
