from typing import Tuple, Iterable
from numbers import Real
import matplotlib as mpl


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
        raise ValueError(f"{position=}")
    transform = ax.figure.dpi_scale_trans + mpl.transforms.ScaledTranslation(
        *position, position_transform
    )
    polygon = mpl.patches.Polygon(shape, transform=transform, **kw,)
    ax.add_patch(polygon)
