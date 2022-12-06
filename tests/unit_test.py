import pytest
import lakeview as lv


def test_load_bam():
    p = lv.SequenceAlignment.from_file("data/SKBR3_Illumina_550bp_pcrFREE.bam")
    assert len(p.segments) == 1800

    p = lv.SequenceAlignment.from_file(
        "data/SKBR3_Illumina_550bp_pcrFREE.bam",
        reference_name="17",
        start=64041800,
        end=64043800,
    )
    assert len(p.segments) == 912

    with pytest.warns(UserWarning):
        p = lv.SequenceAlignment.from_file(
            "data/SKBR3_Illumina_550bp_pcrFREE.bam",
            reference_name="1",
        )
    assert len(p.segments) == 0


def test_filter_segments():
    pass

