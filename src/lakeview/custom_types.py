from typing import Tuple, Hashable, FrozenSet, Union, TypeVar

NativeHashable = Union[int, float, str, Tuple[Hashable], FrozenSet]
GroupIdentifier = TypeVar("GroupIdentifier", bound=NativeHashable)
LinkIdentifier = TypeVar("LinkIdentifier", bound=NativeHashable)
Color = TypeVar("Color", bound=NativeHashable) # TODO: more specific color type
Position = Union[int, float]