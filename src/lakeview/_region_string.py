#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations


def get_region_string(sequence_name: str, start: int, end: int) -> str:
    """
    Get a samtools-compatible region string from `sequence_name`, `start`, and `end`.
    Note that `start` and `end` coordinates are 0-based half-open intervals, while the samtools-compatible notation instead represent a 1-based closed interval.
    See https://pysam.readthedocs.io/en/latest/glossary.html#term-region

    >>> get_region_string("chr1", 15000, 20000)
    'chr1:15001-20000'
    """
    return f"{sequence_name}:{start+1}-{end}"


def parse_region_string(region_string: str) -> tuple[str, int, int]:
    """
    Parse `region_string` into (sequence_name, start, end).
    Raises ValueError if `region_string` is not correctly formatted.

    >>> parse_region_string("chr1:15001-20000")
    ('chr1', 15000, 20000)
    >>> parse_region_string("chr14:104,586,347-107,043,718")
    ('chr14', 104586346, 107043718)
    """
    sequence_name, coordinate_str = region_string.split(":")
    error_message = f"Invalid `region_string`: {region_string!r}. Expecting '<sequence_name>:<start>-<end>'."
    if not sequence_name:
        raise ValueError(error_message)
    start_str, end_str = coordinate_str.split("-")
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
    """
    sequence_name, start, end = parse_region_string(region_string)
    return get_region_string(sequence_name, start, end)
