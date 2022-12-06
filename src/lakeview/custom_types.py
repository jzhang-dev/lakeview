from typing import Tuple, Hashable, FrozenSet, Union, TypeVar
from matplotlib import axes

NativeHashable = Union[int, float, str, Tuple[Hashable], FrozenSet]
GroupIdentifier = NativeHashable
LinkIdentifier = NativeHashable
Color = Union[Tuple[float, float, float], Tuple[float, float, float, float], str]
Position = Union[int, float]
Axes = axes.Axes
