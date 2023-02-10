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
"A color represented as one of the formats supported by Matplotlib. See https://matplotlib.org/stable/tutorials/colors/colors.html"
Point = Tuple[float, float]
"A point represented as tuple(x, y)."
Line = Sequence[Point]
Figure = mpl.figure.Figure
Axes = mpl.axes.Axes
Base = Literal["A", "T", "C", "G"]
"a DNA base represented as a string"
