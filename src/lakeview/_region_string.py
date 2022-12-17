#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Optional


def get_region_string(
    sequence_name: str, start: Optional[int] = None, end: Optional[int] = None
) -> str:
    """
    Get a samtools-compatible region string from `sequence_name`, `start`, and `end`.
    Note that `start` and `end` coordinates are 0-based half-open intervals, while the samtools-compatible notation instead represent a 1-based closed interval.
    See https://pysam.readthedocs.io/en/latest/glossary.html#term-region

    >>> get_region_string("chr1", 15000, 20000)
    'chr1:15001-20000'
    >>> get_region_string("chr1")
    'chr1'
    """
    if start is None and end is None:
        return sequence_name
    elif start is not None and end is not None:
        return f"{sequence_name}:{start+1}-{end}"
    else:
        raise TypeError("`start` and `end` must be both integers or None")


def parse_region_string(region_string: str) -> tuple[str, Optional[int], Optional[int]]:
    """
    Parse `region_string` into (sequence_name, start, end).
    If `region_string` only contains sequence name, both `start` and `end` will be None.
    Raises ValueError if `region_string` is not correctly formatted.

    >>> parse_region_string("chr1:15001-20000")
    ('chr1', 15000, 20000)
    >>> parse_region_string("chr14:104,586,347-107,043,718")
    ('chr14', 104586346, 107043718)
    >>> parse_region_string("chr14")
    ('chr14', None, None)
    >>> parse_region_string("chr14:200")
    Traceback (most recent call last):
        ...
    ValueError: Invalid `region_string`: 'chr14:200'
    """
    if ":" not in region_string:
        return region_string, None, None
    error_message = f"Invalid `region_string`: {region_string!r}"
    try:
        sequence_name, coordinate_str = region_string.split(":")
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
    return sequence_name, start, end


def normalize_region_string(region_string: str) -> str:
    """
    Normalize `region_string` to be samtools-compatible. Commas are removed from the coordinates.
    Raises ValueError if `region_string` is not correctly formatted.

    >>> normalize_region_string("chr14:104,586,347-107,043,718")
    'chr14:104586347-107043718'
    >>> normalize_region_string("chr14")
    'chr14'
    """
    sequence_name, start, end = parse_region_string(region_string)
    return get_region_string(sequence_name, start, end)
