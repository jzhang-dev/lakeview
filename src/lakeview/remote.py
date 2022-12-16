#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
import tempfile
import os
import pathlib
from typing import Literal, Callable, Optional, Union
import pysam
from ._region_string import (
    parse_region_string,
    normalize_region_string,
    get_region_string,
)


def download_bam(
    bam_url: str,
    output_bam_path: str,
    region: Union[str, tuple[str, float, float], None] = None,
    *,
    index_url: Optional[str] = None,
    output_index_path: Optional[str] = None,
    index: bool = True,
    override: bool = False,
) -> None:
    if index_url is None:
        index_url = bam_url + ".bai"
    if output_index_path is None:
        output_index_path = output_bam_path + '.bai'
    region_string: Optional[str]
    if region is None:
        region_string = None
    elif isinstance(region, str):
        region_string = normalize_region_string(region)
    else:
        region_string = get_region_string(*region)
    
    if not os.path.isfile(output_bam_path) or override:
        workdir = os.getcwd()
        with tempfile.TemporaryDirectory() as d:
            try:
                os.chdir(d)
                params = ["-X", "-b", bam_url, index_url]
                if region_string:
                    params.append(region_string)
                bam_data = pysam.view(*params) # type: ignore # pysam.view is not correctly typed. 
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