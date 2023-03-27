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
    p = lv.SequenceAlignment.from_file(
        "tests/data/SKBR3_Illumina_550bp_pcrFREE.bam", CHROMOSOME
    )
    # Group by strand
    group_identifiers = p._get_group_identifiers(group_by="strand")
    c = collections.Counter(group_identifiers)
    assert c["forward"] == 905
    assert c["reverse"] == 895

    # Group by proper pair
    group_identifiers = p._get_group_identifiers(group_by="proper_pair")
    c = collections.Counter(group_identifiers)
    assert c[True] == 1692
    assert c[False] == 108

    # Group by custom function
    group_identifiers = p._get_group_identifiers(
        group_by=lambda seg: seg.query_name[-1]
    )
    c = collections.Counter(group_identifiers)
    assert group_identifiers[:10] == ["2", "7", "3", "0", "0", "1", "0", "0", "8", "1"]


def test_filter_segments() -> None:
    CHROMOSOME = "17"
    p = lv.SequenceAlignment.from_file(
        "tests/data/SKBR3_Illumina_550bp_pcrFREE.bam", CHROMOSOME
    )
    # Filter by proper pair
    filter_keys = p._get_filter_keys(filter_by=lambda seg: seg.is_proper_pair)
    assert sum(filter_keys) == 1692


def test_link_segments() -> None:
    CHROMOSOME = "17"
    p = lv.SequenceAlignment.from_file("tests/data/SKBR3_PacBio.bam", CHROMOSOME)

    # Link by query name
    link_identifiers = p._get_link_identifiers(link_by="name")
    c = collections.Counter(link_identifiers)
    assert len(c) == 184
    # Link by read pair
    link_identifiers = p._get_link_identifiers(link_by="pair")
    c = collections.Counter(link_identifiers)
    assert len(c) == 184


def test_sort_segments() -> None:
    CHROMOSOME = "17"
    p = lv.SequenceAlignment.from_file("tests/data/SKBR3_PacBio.bam", CHROMOSOME)

    # Sort by segment length
    sort_keys = p._get_sort_keys(sort_by="length")
    sorted_segments = lv._layout.key_sort(p.segments, sort_keys)
    assert [seg.query_alignment_length for seg in sorted_segments[:5]] == [
        31541,
        29597,
        26793,
        23539,
        22501,
    ]
    assert (
        sorted_segments[0].query_alignment_length
        >= sorted_segments[10].query_alignment_length
        >= sorted_segments[50].query_alignment_length
        >= sorted_segments[-1].query_alignment_length
    )


def test_get_aligned_query_sequence() -> None:
    p = lv.SequenceAlignment.from_file(
        "tests/data/SKBR3_PacBio.bam", region="17:64,042,916-64,043,519"
    )
    counter = collections.Counter(
        seg.get_aligned_query_sequence(64_043_248) for seg in p.segments
    )
    assert counter["T"] == 91 and counter["G"] == 52


# TODO: fix this test
# def test_genome_viewer_widget_interactivity() -> None:
#     p = lv.SequenceAlignment.from_file(
#         "tests/data/SKBR3_Illumina_550bp_pcrFREE.bam",
#         region=("17", (64041800, 64043800)),
#     )
#     gv = lv.GenomeViewer(2)
#     p.draw_pileup(gv.axes[0])
#     p.draw_alignment(gv.axes[1], show_mismatches=False)
#     gv.set_xlim(64041800, 64043400)
#     assert gv.get_xlim() == (64041800.0, 64043400.0)
#     widget = gv.widget
#     widget._shift(0.5)
#     assert gv.get_xlim() == (64042600.0, 64044200.0)
#     widget._shift(-0.5)
#     assert gv.get_xlim() == (64041800, 64043400)
#     widget._zoom(0.5)
#     assert gv.get_xlim() == (64042400.0, 64042800.0)
#     widget._zoom(3)
#     assert gv.get_xlim() == (64040800.0, 64044400.0)
#     widget._goto(64041800, 64043400)
#     assert gv.get_xlim() == (64041800, 64043400)
#     widget._goto_region("64041800-64043500")
#     assert gv.get_xlim() == (64041800, 64043500)
#     widget._goto_region(widget._region_text.value)
#     assert gv.get_xlim() == (64041800, 64043500)



def test_load_bigwig() -> None:
    BIGWIG_PATH = "tests/data/test.bigwig"
    painter = lv.Wiggle.from_bigwig(BIGWIG_PATH, "1")
    assert painter.intervals == [(0, 1), (1, 2), (2, 3), (100, 150), (150, 151)]

    BIGWIG_URL = "https://github.com/deeptools/pyBigWig/raw/96b951c9281bfbe5358677d532fba2342bd2323f/pyBigWigTest/test.bw"
    painter = lv.Wiggle.from_bigwig(BIGWIG_URL, "1")
    assert painter.intervals == [(0, 1), (1, 2), (2, 3), (100, 150), (150, 151)]



