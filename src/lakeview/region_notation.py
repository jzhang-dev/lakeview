#!/usr/bin/env python
# coding: utf-8

"""
In Lakeview, the genomic region of interest can be specified in one of the following two formats:

- a region notation string (e.g. ``"chr1:15,000-20,000"``), with optional commas as thousands separators
- a nested tuple in the form of ``(sequence_name, (start, end))`` (e.g. ``('chr1', (15_000, 20_000))``), with optional underscores as thousands separators following `the Python standard syntax <https://peps.python.org/pep-0515/#underscore-placement-rules>`_.

In both cases, the start and end coordinates can be omitted together to include the entire sequence (e.g. ``"chr1"`` or ``('chr1', (None, None))``).

This module contains functions for the parsing and conversion of valid region notations.

"""

from __future__ import annotations
from typing import Optional, Union, Any
from math import floor, ceil


def get_region_notation(
    sequence_name, interval: Optional[tuple[float, float]] = None
) -> str:
    """
    Returns a samtools-compatible region string from ``sequence_name``, ``start``, and ``end``.
    Note that ``start`` and ``end`` coordinates are 0-based half-open intervals, while the samtools-compatible notation instead represent a 1-based closed interval.
    See `pysam documentation <https://pysam.readthedocs.io/en/latest/glossary.html#term-region>`_.

    >>> get_region_notation("chr1", (15000, 20000))
    'chr1:15001-20000'
    >>> get_region_notation("chr1")
    'chr1'
    >>> get_region_notation("chr1", (234.2, 480.8))
    'chr1:235-481'
    """
    if interval is None:
        return sequence_name
    else:
        start, end = interval
        start = floor(start)
        end = ceil(end)
        return f"{sequence_name}:{start+1}-{end}"


def parse_region_notation(
    region_notation: str,
) -> tuple[str, Optional[tuple[int, int]]]:
    """
    Parse ``region_notation`` into ``(sequence_name, (start, end))``.
    If ``region_notation`` only contains the sequence name, both ``start`` and ``end`` will be None.

    :param region_notation: region notation string to be parsed
    :returns: a two-element tuple ``(sequence_name, (start, end))``
    :raises ValueError: if ``region_notation`` is not correctly formatted



    >>> parse_region_notation("chr1:15001-20000")
    ('chr1', (15000, 20000))
    >>> parse_region_notation("chr14:104,586,347-107,043,718")
    ('chr14', (104586346, 107043718))
    >>> parse_region_notation("chr14")
    ('chr14', None)
    >>> parse_region_notation("chr14:200")
    Traceback (most recent call last):
        ...
    lakeview.region_notation.InvalidRegionNotationError: Unable to parse region: 'chr14:200'.
    """
    if ":" not in region_notation:
        return region_notation, None
    try:
        sequence_name, coordinate_str = region_notation.split(":")
        start_str, end_str = coordinate_str.split("-")
    except ValueError:
        raise InvalidRegionNotationError(region_notation)
    if not sequence_name:
        raise InvalidRegionNotationError(region_notation)
    start_str = "".join(x for x in start_str if x != ",")
    end_str = "".join(x for x in end_str if x != ",")
    if not start_str.isnumeric() or not end_str.isnumeric():
        raise InvalidRegionNotationError(region_notation)
    start, end = int(start_str) - 1, int(end_str)
    return sequence_name, (start, end)


def normalize_region_notation(region_notation: str) -> str:
    """
    Normalize ``region_notation`` to be samtools-compatible. Spaces are removed. Commas are removed from the coordinates.

    :raises RegionNotationError: if `region_notation` is not correctly formatted

    >>> normalize_region_notation("chr14:104,586,347-107,043,718")
    'chr14:104586347-107043718'
    >>> normalize_region_notation("chrX: 26,23,922 - 26,235,150")
    'chrX:2623922-26235150'
    >>> normalize_region_notation("chr14")
    'chr14'
    """
    sequence_name, interval = parse_region_notation(
        region_notation.replace(" ", "").replace("\t", "")
    )
    return get_region_notation(sequence_name, interval)


class InvalidRegionNotationError(TypeError):
    """Exception raised for invalid region of interest format."""

    def __init__(self, region: Any):
        message = f"Unable to parse region: {region!r}."
        super().__init__(message)
