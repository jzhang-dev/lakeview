#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations



def get_region_string(sequence_name: str, start: float, end: float) -> str:
    """
    Get samtools-compatible region string from `sequence_name`, `start`, and `end`.
    Float `start` and `end` coordinates are rounded to nearest integers. 

    >>> get_region_string("chr14", 104586347, 107043718)
    "chr14:104586347-107043718"
    """
    return f"{sequence_name}:{round(start)}-{round(end)}"

def parse_region_string(region_string: str) -> tuple[str, float, float]:
    """
    Parse `region_string` into (sequence_name, start, end).

    >>> parse_region_string("chr14:104,586,347-107,043,718")
    ("chr14", 104586347, 107043718)
    """
    sequence_name, coordinate_str = region_string.split(":")
    error_message = f"Invalid `region_string`: {region_string!r}. Expecting '<sequence_name>:<start>-<end>'."
    if not sequence_name:
        raise ValueError(error_message)
    start_str, end_str = coordinate_str.split('-')
    start_str = "".join(x for x in start_str if x != ',')
    end_str = "".join(x for x in end_str if x != ',')
    if not start_str.isnumeric() or not end_str.isnumeric():
        raise ValueError(error_message)
    start, end = float(start_str), float(end_str)
    return sequence_name, start, end


def normalize_region_string(region_string: str) -> str:
    """
    Normalize `region_string` to be samtools-compatible. 
    Commas are removed from the coordinates. Float `start` and `end` coordinates are rounded to nearest integers. 

    >>> normalize_region_string("chr14:104,586,347-107,043,718")
    "chr14:104586347-107043718"
    >>> normalize_region_string("ref:2.1-5.5")
    "ref:2-6"
    """
    sequence_name, start, end = parse_region_string(region_string)
    return get_region_string(sequence_name, start, end)

