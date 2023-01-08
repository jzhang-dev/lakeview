#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Union, Literal, Tuple, Sequence
from collections.abc import Hashable
import matplotlib as mpl

Identifier = Union[int , float , str , Tuple[Hashable, ...]]
"Identifier"
LinkIdentifier = Identifier
"Link identifier"
GroupIdentifier = Identifier
"Group identifier"
Color = Union[Tuple[float, float, float], Tuple[float, float, float, float], str]
"See https://matplotlib.org/stable/tutorials/colors/colors.html"
Position = Union[int, float]
Point = Tuple[float, float]
Line = Sequence[Point]
Figure = mpl.figure.Figure
Axes = mpl.axes.Axes
Base = Literal["A", "T", "C", "G"]
