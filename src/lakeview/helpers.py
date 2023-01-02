#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Literal, Callable, Any, Optional, TypeVar
from collections.abc import Iterable, Sequence
import os
import pathlib
import tempfile
from numbers import Real
import matplotlib as mpl
from matplotlib.colors import rgb2hex
import pysam
from ._type_alias import Color, Identifier

T = TypeVar("T")

def key_sort(sequence: Sequence[T], keys: Sequence[Identifier], reverse: bool = False) -> Sequence[T]:
    """
    Sort `sequence` by the values of `keys`.

    >>> x = ['a', 'b', 'c']
    >>> keys = [3, 1, 2]
    >>> key_sort(x, keys=keys)
    ['b', 'c', 'a']
    """
    if len(sequence) != len(keys):
        raise ValueError(f"`sequence` and `keys` have different lengths: {len(sequence)=}; {len(keys)=}.")
    if len(sequence) == 0:
        return []
    zipped_sequence = list(zip(keys, sequence))
    sorted_zipped_sequence = sorted(zipped_sequence, key=lambda pair: pair[0], reverse=reverse)
    sorted_sequence = [pair[1] for pair in sorted_zipped_sequence]
    return sorted_sequence


def key_filter(sequence: Sequence[T], keys: Sequence[bool]) -> Sequence[T]:
    """
    Filter `sequence` by the values of `keys`.

    >>> x = ['a', 'b', 'c']
    >>> keys = [True, False, True]
    >>> key_filter(x, keys=keys)
    ['a', 'c']
    """
    if len(sequence) != len(keys):
        raise ValueError(f"`sequence` and `keys` have different lengths: {len(sequence)=}; {len(keys)=}.")
    return [element for element, key in zip(sequence, keys) if key]


def sort_by_keys(
    *iterables: Sequence, keys: Sequence[Identifier], reverse: bool = False
) -> list[list]:
    """
    Sort multiple equal-length sequences by the value of another sequences.

    >>> a = ['a', 'b', 'c']
    >>> b = ['x', 'y', 'z']
    >>> r = [3, 1, 2]
    >>> sort_by_keys(a, b, keys=r, reverse=False)
    [('b', 'c', 'a'), ('y', 'z', 'x')]
    """
    # Check iterable lengths are equal
    lengths = [len(x) for x in iterables] + [len(keys)]
    if len(set(lengths)) > 1:
        raise ValueError("All iterables and `keys` must have the same length.")
    # Return early if iterables are empty
    if any(l == 0 for l in lengths):
        return [[] for x in iterables]

    zipped_lists = list(zip(*iterables))
    sorted_zipped_lists = [
        x
        for (_, x) in sorted(
            zip(keys, zipped_lists), key=lambda pair: pair[0], reverse=reverse
        )
    ]
    sorted_lists = list(zip(*sorted_zipped_lists))
    return sorted_lists


def filter_by_keys(*iterables: Sequence, keys: Sequence[bool]) -> list[list]:
    """
    Filter multiple equal-length lists by the value of another list.

    >>> a = ['a', 'b', 'c']
    >>> b = ['x', 'y', 'z']
    >>> f = [True, False, True]
    >>> filter_by_keys(a, b, keys=f)
    [('a', 'c'), ('x', 'z')]
    """
    # Check iterable lengths are equal
    lengths = [len(x) for x in iterables] + [len(keys)]
    if len(set(lengths)) > 1:
        raise ValueError("All iterables and `keys` must have the same length.")
    # Return early if iterables are empty
    if any(l == 0 for l in lengths):
        return [[] for x in iterables]

    zipped_lists = list(zip(*iterables))
    filtered_zipped_lists = [x for (x, b) in zip(zipped_lists, keys) if b]
    filtered_lists = list(zip(*filtered_zipped_lists))
    return filtered_lists





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
