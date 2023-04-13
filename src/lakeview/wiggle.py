from __future__ import annotations
import pyBigWig
from typing import Sequence
from dataclasses import dataclass, field
from functools import cached_property
from .region_notation import (
    parse_region_notation,
    get_region_notation,
    InvalidRegionNotationError,
)
from ._type_alias import Region, Axes, Color


@dataclass(init=False)
class Wiggle:
    intervals: Sequence[tuple[int, int]]
    values: Sequence[float]

    def __init__(self, intervals: Sequence[tuple[int, int]], values: Sequence[float]):
        if len(intervals) != len(values):
            raise ValueError()
        self.intervals = intervals
        self.values = values

    @staticmethod
    def _parse_bigwig(
        bigwig: pyBigWig.bigWigFile, region: Region
    ) -> tuple[Sequence[tuple[int, int]], Sequence[float]]:
        # TODO: support bins for performance
        chromosome: str
        interval: tuple[int, int] | None
        if isinstance(region, str):
            chromosome, interval = parse_region_notation(region)
        else:
            chromosome, interval = parse_region_notation(get_region_notation(*region))
        bigwig_intervals_args: list = (
            [chromosome] if interval is None else [chromosome, *interval]
        )
        intervals: list[tuple[int, int]] = []
        values: list[float] = []
        data: Sequence[tuple[int, int, float]] | None = bigwig.intervals(
            *bigwig_intervals_args
        )
        if data is None: # No intervals found in the region
            return intervals, values
        for interval_start, interval_end, value in data:
            intervals.append((interval_start, interval_end))
            values.append(value)
        return intervals, values

    @classmethod
    def from_bigwig(cls, file_path_or_url: str, region: Region):
        with pyBigWig.open(file_path_or_url, "r") as bw:
            if bw is None:
                raise IOError(f"Failed to open {file_path_or_url}")
            intervals, values = cls._parse_bigwig(bw, region=region)
            return cls(intervals=intervals, values=values)


    @cached_property
    def _data(self) -> tuple[Sequence[float], Sequence[float]]:
        if len(self.intervals) == 0:
            return [], []
        xs:list[float] = [self.intervals[0][0] - 0.5]
        ys:list[float] = [0]
        previous_end: int = self.intervals[0][0]
        for (start, end), value in zip(self.intervals, self.values):
            if start != previous_end:
                xs += [xs[-1], start-0.5]
                ys += [0, 0]
            xs += [start-0.5, end-0.5]
            ys += [value, value]
            previous_end = end
        xs.append(end-0.5)
        ys.append(0)
        return xs, ys


    def draw(self, ax: Axes, *, plot_kw={}) -> None:
        xs, ys = self._data
        if len(xs) == 0:
            return
        ax.plot(xs, ys, **plot_kw)
        ax.set_ylim(bottom=0)
