#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Union, Literal, Tuple, Sequence
from collections.abc import Hashable
import matplotlib as mpl

NativeHashable = Union[int, float, str, Tuple[Hashable, ...], frozenset]
GroupIdentifier = NativeHashable
LinkIdentifier = NativeHashable
Color = Union[Tuple[float, float, float], Tuple[float, float, float, float], str]
Position = Union[int, float]
Point = Tuple[float, float]
Line = Sequence[Point]
Figure = mpl.figure.Figure
Axes = mpl.axes.Axes
Base = Literal["A", "T", "C", "G"]
