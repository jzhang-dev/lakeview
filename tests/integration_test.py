import os, gzip
import pytest
import matplotlib.pyplot as plt
import lakeview as lv


def test_SKBR3():
    """
    These two BAM files have MD tag but no =/X CIGAR operations.
    """
    ILLUMINA_BAM_PATH = "tests/data/SKBR3_Illumina_550bp_pcrFREE.bam"
    PACBIO_BAM_PATH = "tests/data/SKBR3_PacBio.bam"
    OUTPUT_PNG_PATH = "tests/output/SKBR3_Illumina_PacBio.png"
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
    gv.set_xlabel("chr17")
    gv.axes[-1].xaxis.set_major_formatter(lv.util.base_formatter(unit="mb", fmt="{:.3f}"))
    gv.set_title("SKBR3")
    gv.savefig(OUTPUT_PNG_PATH, dpi=300)
    assert os.path.getsize(OUTPUT_PNG_PATH) > 200e3


def test_GNAS_WES():
    EXON_BAM_PATH = "tests/data/HG002_GNAS_Illumina_WES.bam"
    OUTPUT_SVG_PATH = "tests/output/GNAS_WES.svg"
    painter = lv.SequenceAlignment.from_file(EXON_BAM_PATH, "rb")
    painter.segments = [seg for i, seg in enumerate(painter.segments) if i % 20 == 0]
    assert len(painter.segments) == 754

    gv = lv.GenomeViewer(2, height_ratios=(1, 8), figsize=(12, 8))
    painter.draw_pileup(gv.axes[0], show_mismatches=False)
    painter.draw_alignment(
        gv.axes[1], show_mismatches=False, show_arrowheads=False, max_group_height=30
    )

    gv.set_xlabel("chr20")
    gv.set_title("HG002 GNAS Illumina WES")
    gv.savefig(OUTPUT_SVG_PATH)
    assert os.path.getsize(OUTPUT_SVG_PATH) > 1e6


def test_IGH():
    CHROMOSOME = "chr14"
    START = 104586347
    END = 107043718
    GENCODE_GTF_PATH = "tests/data/gencode.v40.annotation.gtf.gz"
    PACBIO_BAM_PATH = "tests/data/HG002_IGH_PacBio_CCS.bam"
    OUTPUT_PNG_PATH = "tests/output/IGH_PacBio_Gencode.png"

    with gzip.open(GENCODE_GTF_PATH, "rt") as f:
        gencode_painter = lv.GeneAnnotation.from_file(
            file_object=f,
            format="gtf",
            sequence_name=CHROMOSOME,
            start=START,
            end=END,
        )
    gencode_painter.transcripts.sort(key=len, reverse=True)
    pacbio_painter = lv.SequenceAlignment.from_file(PACBIO_BAM_PATH, "rb")

    gv = lv.GenomeViewer(3, height_ratios=(1, 8, 2))
    pacbio_painter.draw_pileup(
        gv.axes[0],
        show_mismatches=False,
    )
    pacbio_painter.draw_alignment(
        gv.axes[1],
        show_mismatches=False,
        sort_by=lambda seg: -seg.query_alignment_length,
        max_group_height=50,
    )
    gencode_painter.draw_transcripts(
        gv.axes[2], max_group_height=4, labels=lambda t: t.attributes["gene_name"]
    )

    gv.set_xlim((105679000, 105776000))
    gv.set_xlabel(CHROMOSOME)
    gv.set_title("HG002 IGH PacBio + Gencode")
    gv.savefig(OUTPUT_PNG_PATH, dpi=300)
    assert os.path.getsize(OUTPUT_PNG_PATH) > 100e3


def test_SNURF_methylation():
    PACBIO_BAM_PATH = "tests/data/HG002_GRCh38_SNURF_haplotagged.bam"
    OUTPUT_PNG_PATH = "tests/output/SNURF_methylation.png"
    p = lv.SequenceAlignment.from_file(PACBIO_BAM_PATH, "rb")
    fig, ax = plt.subplots(figsize=(8, 5))
    p.draw_alignment(
        ax,
        group_by="haplotype",
        show_modified_bases=True,
        modified_bases_kw=dict(
            linewidth=1.5,
            colormaps={("C", "m", "+"): "cividis", ("C", "m", "-"): "cividis"},
        ),
    )
    ax.set_title("$\it{SNURF}$ differentially methylated region")
    ax.set_xlim(24.953e6, 24.958e6)
    ax.set_xlabel("chr15 (Mb)")
    ax.xaxis.set_major_formatter(lv.util.base_formatter(unit="mb", fmt="{:.3f}"))
    fig.savefig(OUTPUT_PNG_PATH, dpi=300)
    assert os.path.getsize(OUTPUT_PNG_PATH) > 100e3


def test_dot_plot():
    IGH_FASTA_PATH = "tests/data/IGH_reference_sequences.fasta.gz"
    OUTPUT_PNG_PATH = "tests/output/IGH_dot_plot.png"

    with gzip.open(IGH_FASTA_PATH, "rt") as f1, gzip.open(IGH_FASTA_PATH, "rt") as f2:
        painter = lv.DotPlot.from_files(
            f1,
            f2,
            x_sequence_name="NC_000014.9:c106879844-105586437",
            y_sequence_name="NT_187600.1:c1351393-54793",
            sample_fraction=0.2,
        )
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    painter.draw_dots(axes[0])
    painter.draw_heatmap(
        axes[1], bin_size=10e3, cmin=0
    )  # TODO: add a microsatellite region to demonstrate StainedGlass-style plot
    fig.savefig(OUTPUT_PNG_PATH, dpi=300)
    assert os.path.getsize(OUTPUT_PNG_PATH) > 100e3
