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
    

def pack_intervals(intervals: Iterable[tuple[float, float]], max_offset:float=float('inf')) -> Sequence[int]:
    """
    Assign an non-negative offset to each input interval so that intervals sharing the same offset will not overlap with each other, while minimising offset values.
    Intervals are treated as being closed.

    >>> pack_intervals([(1, 2), (3, 4), (1, 3), (0, 5)])
    [0, 0, 1, 2]
    >>> pack_intervals([(1, 2), (1, 2), (1, 2), (3, 4)], max_offset=1)
    [0, 1, -1, 0]
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
            new_offset = offset + 1
            if new_offset > max_offset:
                # Reached max offset; try next interval
                offsets.append(-1) # -1 -> failed to pack within max offset
                continue
            occupied_intervals.append([(start, end)])
            offsets.append(offset + 1)
    return offsets
