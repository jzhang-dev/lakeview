import os
import pytest
import lakeview as lv

def test_SKBR3():
    ILLUMINA_BAM_PATH = "data/SKBR3_Illumina_550bp_pcrFREE.bam"
    PACBIO_BAM_PATH = "data/SKBR3_PacBio.bam"
    OUTPUT_PNG_PATH = "output/SKBR3_Illumina_PacBio.png"
    illumina_painter = lv.SequenceAlignment.from_file(ILLUMINA_BAM_PATH, "rb")
    pacbio_painter = lv.SequenceAlignment.from_file(PACBIO_BAM_PATH, "rb")

    gv = lv.GenomeViewer(4, figsize=(12, 15), height_ratios=(1, 8, 1, 8))
    illumina_painter.draw_pileup(gv.axes[0])
    illumina_painter.draw_alignment(
        gv.axes[1],
        color_by=lambda segment: "lightgray" if segment.is_proper_pair else "firebrick",
        group_by="strand",
        max_group_height=30,
    )
    pacbio_painter.draw_pileup(gv.axes[2])
    pacbio_painter.draw_alignment(gv.axes[3])
    
    gv.set_xlim(64040802, 64045633)
    gv.axes[1].set_ylabel("Illumina")
    gv.axes[3].set_ylabel("PacBio")
    gv.set_xlabel("chr17") # TODO: add x formatter
    gv.set_title("SKBR3")
    gv.savefig(OUTPUT_PNG_PATH, dpi=300)

    assert os.path.getsize(OUTPUT_PNG_PATH) > 200e3
