#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
import tempfile
import os
import pathlib
from typing import Literal, Callable, Optional, Union
import pysam
from .region_notation import (
    parse_region_notation,
    normalize_region_notation,
    get_region_notation,
)


def download_bam(
    bam_url: str,
    output_bam_path: str,
    region: str | tuple[str, tuple[int, int]] | tuple[str, None] | None = None,
    *,
    index_url: Optional[str] = None,
    output_index_path: Optional[str] = None,
    index: bool = True,
    override: bool = False,
) -> None:
    if index_url is None:
        index_url = bam_url + ".bai"
    if output_index_path is None:
        output_index_path = output_bam_path + ".bai"
    region_notation: Optional[str]
    if region is None:
        region_notation = None
    elif isinstance(region, str):
        region_notation = normalize_region_notation(region)
    elif isinstance(region, tuple):
        reference_name, interval = region
        region_notation = get_region_notation(reference_name, interval)
    else:
        raise TypeError(
            f"Invalid type for `region`: {type(region)!r}. Expecting str | tuple[str, tuple[int, int]] | tuple[str, None] | None."
        )

    if not os.path.isfile(output_bam_path) or override:
        workdir = os.getcwd()
        with tempfile.TemporaryDirectory() as d:
            try:
                os.chdir(d)
                params = ["-X", "-b", bam_url, index_url]
                if region_notation:
                    params.append(region_notation)
                bam_data = pysam.view(*params)  # type: ignore # pysam.view is not correctly typed.
            finally:
                os.chdir(workdir)
        target_dir = os.path.split(output_bam_path)[0]
        pathlib.Path(target_dir).mkdir(parents=True, exist_ok=True)
        with open(output_bam_path, "wb") as f:
            f.write(bam_data)
    if index and (not os.path.isfile(output_index_path) or override):
        target_dir = os.path.split(output_index_path)[0]
        pathlib.Path(target_dir).mkdir(parents=True, exist_ok=True)
        pysam.index(output_bam_path, output_index_path)  # type: ignore # pysam.index is not correctly typed.
