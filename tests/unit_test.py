#!/usr/bin/env python
# coding: utf-8

import collections
import pytest
import lakeview as lv


# TODO: load user supplied reference sequence.


def test_load_bam():
    p = lv.SequenceAlignment.from_file("tests/data/SKBR3_Illumina_550bp_pcrFREE.bam")
    assert len(p.segments) == 1800

    p = lv.SequenceAlignment.from_file(
        "tests/data/SKBR3_Illumina_550bp_pcrFREE.bam",
        reference_name="17",
        start=64041800,
        end=64043800,
    )
    assert len(p.segments) == 912

    with pytest.warns(UserWarning):
        p = lv.SequenceAlignment.from_file(
            "tests/data/SKBR3_Illumina_550bp_pcrFREE.bam",
            reference_name="1",
        )
    assert len(p.segments) == 0


def test_group_segments():
    p = lv.SequenceAlignment.from_file("tests/data/SKBR3_Illumina_550bp_pcrFREE.bam")
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


def test_filter_segments():
    p = lv.SequenceAlignment.from_file("tests/data/SKBR3_Illumina_550bp_pcrFREE.bam")
    kw = dict(group_by=None, link_by=None, color_by=None, sort_by=None)

    # Group by proper pair
    segments, links, groups, colors = p._parse_segment_parameters(
        filter_by=lambda seg: seg.is_proper_pair, **kw
    )
    assert len(segments) == 1692


def test_link_segments():
    p = lv.SequenceAlignment.from_file("tests/data/SKBR3_PacBio.bam")
    kw = dict(group_by=None, filter_by=None, color_by=None, sort_by=None)

    # Link by query name
    segments, links, groups, colors = p._parse_segment_parameters(link_by="name", **kw)
    c = collections.Counter(links)
    assert len(c) == 184


def test_sort_segments():
    p = lv.SequenceAlignment.from_file("tests/data/SKBR3_PacBio.bam")
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
