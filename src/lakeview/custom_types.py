from typing import Tuple, Hashable, FrozenSet, Union, TypeVar, Literal, Sequence
from matplotlib import axes

NativeHashable = Union[int, float, str, Tuple[Hashable], FrozenSet]
GroupIdentifier = NativeHashable
LinkIdentifier = NativeHashable
Color = Union[Tuple[float, float, float], Tuple[float, float, float, float], str]
Position = Union[int, float]
Point = Tuple[float, float]
Line = Sequence[Point]
Axes = axes.Axes
Base = Literal["A", "T", "C", "G"]
