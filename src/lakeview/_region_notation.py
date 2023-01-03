#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Optional, Union

# Optional: make this module public


def get_region_notation(sequence_name, interval: Optional[tuple[int, int]] = None) -> str:
    """
    Get a samtools-compatible region string from `sequence_name`, `start`, and `end`.
    Note that `start` and `end` coordinates are 0-based half-open intervals, while the samtools-compatible notation instead represent a 1-based closed interval.
    See https://pysam.readthedocs.io/en/latest/glossary.html#term-region

    >>> get_region_notation("chr1", (15000, 20000))
    'chr1:15001-20000'
    >>> get_region_notation("chr1")
    'chr1'
    """
    if interval is None:
        return sequence_name
    else:
        start, end = interval
        return f"{sequence_name}:{start+1}-{end}"


def parse_region_notation(
    region_notation: str,
) -> tuple[str, Optional[tuple[int, int]]]:
    """
    Parse `region_notation` into (sequence_name, interval).
    If `region_notation` only contains sequence name, both `start` and `end` will be None.
    Raises ValueError if `region_notation` is not correctly formatted.

    >>> parse_region_notation("chr1:15001-20000")
    ('chr1', (15000, 20000))
    >>> parse_region_notation("chr14:104,586,347-107,043,718")
    ('chr14', (104586346, 107043718))
    >>> parse_region_notation("chr14")
    ('chr14', None)
    >>> parse_region_notation("chr14:200")
    Traceback (most recent call last):
        ...
    ValueError: Invalid `region_notation`: 'chr14:200'
    """
    if ":" not in region_notation:
        return region_notation, None
    error_message = f"Invalid `region_notation`: {region_notation!r}"
    try:
        sequence_name, coordinate_str = region_notation.split(":")
        start_str, end_str = coordinate_str.split("-")
    except ValueError:
        raise ValueError(error_message)
    if not sequence_name:
        raise ValueError(error_message)
    start_str = "".join(x for x in start_str if x != ",")
    end_str = "".join(x for x in end_str if x != ",")
    if not start_str.isnumeric() or not end_str.isnumeric():
        raise ValueError(error_message)
    start, end = int(start_str) - 1, int(end_str)
    return sequence_name, (start, end)


def normalize_region_notation(region_notation: str) -> str:
    """
    Normalize `region_notation` to be samtools-compatible. Commas are removed from the coordinates.
    Raises ValueError if `region_notation` is not correctly formatted.

    >>> normalize_region_notation("chr14:104,586,347-107,043,718")
    'chr14:104586347-107043718'
    >>> normalize_region_notation("chr14")
    'chr14'
    """
    sequence_name, interval = parse_region_notation(region_notation)
    return get_region_notation(sequence_name, interval)
