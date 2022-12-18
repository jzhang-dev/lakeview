#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
import collections
import tempfile
import os
import pytest
import lakeview as lv


# TODO: load user supplied reference sequence.


def test_load_local_bam() -> None:
    CHROMOSOME = "17"
    p = lv.SequenceAlignment.from_file(
        "tests/data/SKBR3_Illumina_550bp_pcrFREE.bam", CHROMOSOME
    )
    assert len(p.segments) == 1800

    p = lv.SequenceAlignment.from_file(
        "tests/data/SKBR3_Illumina_550bp_pcrFREE.bam",
        region=("17", (64041800, 64043800)),
    )
    assert len(p.segments) == 912

    with pytest.raises(ValueError):
        p = lv.SequenceAlignment.from_file(
            "tests/data/SKBR3_Illumina_550bp_pcrFREE.bam",
            region="1",
        )


def test_load_remote_bam() -> None:
    BAM_URL = "https://s3.amazonaws.com/igv.org.demo/SKBR3/SKBR3_550bp_pcrFREE_S1_L001_AND_L002_R1_001.101bp.bwamem.ill.mapped.sort.bam"
    p = lv.SequenceAlignment.from_remote(BAM_URL, region="17:64040802-64045633")
    assert len(p.segments) == 1800


def test_download_bam() -> None:
    BAM_URL = "https://s3.amazonaws.com/igv.org.demo/SKBR3/SKBR3_550bp_pcrFREE_S1_L001_AND_L002_R1_001.101bp.bwamem.ill.mapped.sort.bam"
    with tempfile.TemporaryDirectory() as d:
        output_file = os.path.join(d, "test.bam")
        lv.remote.download_bam(
            BAM_URL, region="17:64040802-64045633", output_bam_path=output_file
        )
        assert os.path.getsize(output_file) > 100000


def test_group_segments() -> None:
    CHROMOSOME = "17"
    p = lv.SequenceAlignment.from_file("tests/data/SKBR3_Illumina_550bp_pcrFREE.bam", CHROMOSOME)
    kw = dict(filter_by=None, link_by=None, color_by=None, sort_by=None)
    # Group by strand
    segments, links, groups, colors = p._parse_segment_parameters(
        group_by="strand", **kw
    )
    c = collections.Counter(groups)
    assert c["forward"] == 905
    assert c["reverse"] == 895


    # Group by proper pair
    segments, links, groups, colors = p._parse_segment_parameters(
        group_by="proper_pair", **kw
    )
    c = collections.Counter(groups)
    assert c[True] == 1692
    assert c[False] == 108

    # Group by custom function
    segments, links, groups, colors = p._parse_segment_parameters(
        group_by=lambda seg: seg.query_name[-1], **kw
    )
    c = collections.Counter(groups)
    assert groups[:10] == ["2", "7", "3", "0", "0", "1", "0", "0", "8", "1"]


def test_filter_segments() -> None:
    CHROMOSOME = "17"
    p = lv.SequenceAlignment.from_file("tests/data/SKBR3_Illumina_550bp_pcrFREE.bam", CHROMOSOME)
    kw = dict(group_by=None, link_by=None, color_by=None, sort_by=None)

    # Group by proper pair
    segments, links, groups, colors = p._parse_segment_parameters(
        filter_by=lambda seg: seg.is_proper_pair, **kw
    )
    assert len(segments) == 1692


def test_link_segments() -> None:
    CHROMOSOME = "17"
    p = lv.SequenceAlignment.from_file("tests/data/SKBR3_PacBio.bam", CHROMOSOME)
    kw = dict(group_by=None, filter_by=None, color_by=None, sort_by=None)

    # Link by query name
    segments, links, groups, colors = p._parse_segment_parameters(link_by="name", **kw)
    c = collections.Counter(links)
    assert len(c) == 184


def test_sort_segments() -> None:
    CHROMOSOME = "17"
    p = lv.SequenceAlignment.from_file("tests/data/SKBR3_PacBio.bam", CHROMOSOME)
    kw = dict(group_by=None, filter_by=None, color_by=None, link_by=None)

    # Link by segment length
    segments, links, groups, colors = p._parse_segment_parameters(
        sort_by="length", **kw
    )
    assert [seg.query_alignment_length for seg in segments[:5]] == [
        31541,
        29597,
        26793,
        23539,
        22501,
    ]
    assert (
        segments[0].query_alignment_length
        >= segments[10].query_alignment_length
        >= segments[50].query_alignment_length
        >= segments[-1].query_alignment_length
    )


def test_get_aligned_query_sequence() -> None:
    p = lv.SequenceAlignment.from_file(
        "tests/data/SKBR3_PacBio.bam", region="17:64,042,916-64,043,519"
    )
    counter = collections.Counter(
        seg.get_aligned_query_sequence(64_043_248) for seg in p.segments
    )
    assert counter["T"] == 91 and counter["G"] == 52
