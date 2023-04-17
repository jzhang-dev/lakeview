#!/usr/bin/env python
# coding: utf-8

import os, gzip
import pytest
import matplotlib.pyplot as plt
import lakeview as lv
from lakeview.plot import BasePairFormatter


def test_readme_demo():
    # Import Matplotlib and Lakeview
    import matplotlib.pyplot as plt
    import lakeview as lv

    # Load aligned segments in a selected region from a BAM file
    painter = lv.SequenceAlignment.from_file(
        "tests/data/HG002_IGH_PacBio_CCS.bam", region="chr14:105,660,000-105,780,000"
    )
    # Create an empty GenomeViewer with one track
    gv = lv.GenomeViewer(tracks=2, figsize=(8, 5), height_ratios=(1, 4))
    # Plot alignment pileup
    painter.draw_pileup(
        gv.axes[0], # Plot on the first track of the GenomeViewer
        show_mismatches=False,  # Do not highlight mismatched bases
    )
    # Plot aligned segments
    painter.draw_alignment(
        gv.axes[1],  # Plot on the second track of the GenomeViewer
        show_mismatches=False,  # Do not highlight mismatched bases
        sort_by="length",  # Plot longer reads first
        link_by="name",  # Link primary and supplementary alignments of the same read
        max_rows=30,  # Only show the first 30 alignment rows
    )
    # Adjust x axis limits
    gv.set_xlim(105_670_000, 105_777_000)
    # Save the plot
    gv.savefig("tests/output/readme_demo.png")
    assert os.path.getsize("tests/output/readme_demo.png") > 100e3


def test_SKBR3():
    """
    These two BAM files have MD tag but no =/X CIGAR operations.
    """
    CHROMOSOME = "17"
    ILLUMINA_BAM_PATH = "tests/data/SKBR3_Illumina_550bp_pcrFREE.bam"
    PACBIO_BAM_PATH = "tests/data/SKBR3_PacBio.bam"

    illumina_painter = lv.SequenceAlignment.from_file(
        ILLUMINA_BAM_PATH, region=CHROMOSOME
    )
    pacbio_painter = lv.SequenceAlignment.from_file(PACBIO_BAM_PATH, region=CHROMOSOME)

    gv = lv.GenomeViewer(4, figsize=(12, 15), height_ratios=(1, 8, 1, 8))
    illumina_painter.draw_pileup(gv.axes[0])
    illumina_painter.draw_alignment(
        gv.axes[1],
        color_by="proper_pair",
        group_by="strand",
        max_rows=30,
    )
    pacbio_painter.draw_pileup(gv.axes[2])
    pacbio_painter.draw_alignment(gv.axes[3])

    gv.set_xlim(64040802, 64045633)
    gv.axes[1].set_ylabel("Illumina")
    gv.axes[3].set_ylabel("PacBio")
    gv.set_xlabel("chr17")
    gv.axes[-1].xaxis.set_major_formatter(BasePairFormatter("Mb"))
    gv.set_title("SKBR3")

    OUTPUT_PNG_PATH = "tests/output/SKBR3_Illumina_PacBio.png"
    gv.savefig(OUTPUT_PNG_PATH, dpi=300)
    assert os.path.getsize(OUTPUT_PNG_PATH) > 200e3

    OUTPUT_PNG_PATH = "tests/output/SKBR3_PacBio_SNP_tagging.png"
    gv = lv.GenomeViewer(2, figsize=(12, 12), height_ratios=(1, 8))
    pacbio_painter.draw_pileup(gv.axes[0])

    def snp_tag(segment, position):
        qry_seq = segment.get_aligned_query_sequence(position)
        if qry_seq in ("T", "G"):
            return qry_seq
        return ""

    pacbio_painter.draw_alignment(
        gv.axes[1],
        group_by=lambda seg: snp_tag(seg, 64_043_248),
        show_mismatches=True,
        mismatches_kw=dict(linewidth=3),
        max_rows=20,
    )
    gv.set_xlim(64042997, 64043297)
    gv.savefig(OUTPUT_PNG_PATH, dpi=300)
    assert os.path.getsize(OUTPUT_PNG_PATH) > 50e3


def test_GNAS_WES():
    CHROMOSOME = "20"
    EXON_BAM_PATH = "tests/data/HG002_GNAS_Illumina_WES.bam"
    OUTPUT_PNG_PATH = "tests/output/GNAS_WES.png"
    OUTPUT_SVG_PATH = "tests/output/GNAS_WES.svg"

    painter = lv.SequenceAlignment.from_file(EXON_BAM_PATH, CHROMOSOME)
    painter.segments = [seg for i, seg in enumerate(painter.segments) if i % 20 == 0]
    assert len(painter.segments) == 754

    gv = lv.GenomeViewer(2, height_ratios=(1, 8), figsize=(12, 8))
    painter.draw_pileup(gv.axes[0], show_mismatches=False, window_size=100)
    painter.draw_alignment(
        gv.axes[1], show_mismatches=False, show_arrowheads=False, max_rows=30
    )
    gv.set_xlabel("chr20")
    gv.set_title("HG002 GNAS Illumina WES")
    gv.savefig(OUTPUT_PNG_PATH, dpi=300)
    gv.savefig(OUTPUT_SVG_PATH, dpi=300)

    assert os.path.getsize(OUTPUT_PNG_PATH) > 10e3
    assert 10e3 < os.path.getsize(OUTPUT_SVG_PATH) < 10e6


def test_IGH():
    CHROMOSOME = "chr14"
    START = 104586347
    END = 107043718
    GENCODE_GTF_PATH = "tests/data/gencode.v40.annotation.gtf.gz"
    PACBIO_BAM_PATH = "tests/data/HG002_IGH_PacBio_CCS.bam"
    OUTPUT_PNG_PATH = "tests/output/IGH_PacBio_Gencode.png"

    with gzip.open(GENCODE_GTF_PATH, "rt") as f:
        gencode_painter = lv.GeneAnnotation.from_gencode(
            f, format_="gtf", region=(CHROMOSOME, (START, END))
        )
    pacbio_painter = lv.SequenceAlignment.from_file(
        PACBIO_BAM_PATH, region=(CHROMOSOME, (START, END))
    )

    gv = lv.GenomeViewer(3, height_ratios=(1, 8, 2))
    pacbio_painter.draw_pileup(
        gv.axes[0],
        show_mismatches=False,
    )
    pacbio_painter.draw_alignment(
        gv.axes[1],
        show_mismatches=False,
        sort_by="length",
        link_by="name",
        max_rows=50,
    )
    gencode_painter.draw_transcripts(
        gv.axes[2], max_rows=5, sort_by="length", arrows_kw=dict(style="single")
    )

    gv.set_xlim((105679000, 105776000))
    gv.set_xlabel(CHROMOSOME)
    gv.set_title("HG002 IGH PacBio + Gencode")
    gv.xaxis.set_major_formatter(BasePairFormatter("Mb"))
    gv.savefig(OUTPUT_PNG_PATH, dpi=300)
    assert os.path.getsize(OUTPUT_PNG_PATH) > 100e3

    gv = lv.GenomeViewer(1, figsize=(8, 2))
    gencode_painter.draw_genes(
        gv.axes[0], label_by=lambda gene: gene.attributes['gene_name'], labels_kw=dict(size=10, rotation=15, verticalalignment='top')
    )
    gv.set_xlim((105679000, 105776000))
    gv.set_title("IGH genes")
    gv.savefig("tests/output/IGH_genes.png", dpi=300)
    assert os.path.getsize("tests/output/IGH_genes.png") > 10e3



def test_GAPDH_RNAseq():
    CHROMOSOME, START, END = "chr12", int(6.643e6), int(6.648e6)
    RNA_BAM_PATH = "tests/data/GM12878_RNAseq_GAPDH.sample=0.002.bam"
    REFSEQ_GFF_PATH = "tests/data/Refseq_GRCh37_genomic_annotation.gff.gz"
    OUTPUT_PNG_PATH = "tests/output/GM12878_RNAseq_GAPDH.png"

    alignment_painter = lv.SequenceAlignment.from_file(RNA_BAM_PATH, CHROMOSOME)
    with gzip.open(REFSEQ_GFF_PATH, "rt") as f:
        annotation_painter = lv.GeneAnnotation.from_refseq(
            f, region=("NC_000012.11", (START, END))
        )

    gv = lv.GenomeViewer(3, figsize=(8, 8), height_ratios=(1, 7, 2))
    alignment_painter.draw_pileup(gv.axes[0], show_mismatches=False)
    alignment_painter.draw_alignment(
        gv.axes[1],
        show_arrowheads=False,
        show_soft_clippings=False,
        show_hard_clippings=False,
        show_mismatches=True,
        max_rows=50,
        show_group_labels=False,
        show_group_separators=False,
    )
    annotation_painter.draw_transcripts(gv.axes[2], arrows_kw=dict(style="fishbone"))
    gv.set_xlim(6.643e6, 6.648e6)
    gv.set_xlabel(f"{CHROMOSOME}")
    gv.axes[-1].xaxis.set_major_formatter(BasePairFormatter("Mb"))
    gv.set_title(r"GM12878 $\it{GAPDH}$ RNAseq")
    gv.savefig(OUTPUT_PNG_PATH, dpi=300)
    assert os.path.getsize(OUTPUT_PNG_PATH) > 80e3


def test_SNURF_methylation():
    CHROMOSOME = "chr15"
    PACBIO_BAM_PATH = "tests/data/HG002_GRCh38_SNURF_haplotagged.bam"
    OUTPUT_PNG_PATH = "tests/output/SNURF_methylation.png"
    p = lv.SequenceAlignment.from_file(PACBIO_BAM_PATH, CHROMOSOME, load_pileup=False)
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
    ax.set_title(r"$\it{SNURF}$ differentially methylated region")
    ax.set_xlim(24.953e6, 24.958e6)
    ax.set_xlabel("chr15 (Mb)")
    ax.xaxis.set_major_formatter(BasePairFormatter("Mb", show_suffix=False))
    fig.savefig(OUTPUT_PNG_PATH, dpi=300)
    assert os.path.getsize(OUTPUT_PNG_PATH) > 100e3


def test_dot_plot():
    IGH_FASTA_PATH = "tests/data/IGH_reference_sequences.fasta.gz"
    OUTPUT_PNG_PATH = "tests/output/IGH_dot_plot.png"

    x_sequence_name = "NC_000014.9:c106879844-105586437"
    y_sequence_name = "NT_187600.1:c1351393-54793"
    with gzip.open(IGH_FASTA_PATH, "rt") as f1:
        x_sequence = lv.sequence.load_sequence(f1, x_sequence_name)
    with gzip.open(IGH_FASTA_PATH, "rt") as f2:
        y_sequence = lv.sequence.load_sequence(f2, y_sequence_name)
    painter = lv.DotPlot.from_sequences(
        x_sequence, y_sequence, k=50, sample_fraction=0.2
    )
    fig, ax = plt.subplots(figsize=(6, 6))
    painter.draw_dots(ax)
    ax.set_xlabel(x_sequence_name)
    ax.set_ylabel(y_sequence_name)

    # TODO: add a microsatellite region to demonstrate StainedGlass-style plot
    fig.savefig(OUTPUT_PNG_PATH, dpi=300)
    assert os.path.getsize(OUTPUT_PNG_PATH) > 100e3


def test_bigwig():
    BIGWIG_PATH = "tests/data/test.bigwig"
    painter = lv.Wiggle.from_bigwig(BIGWIG_PATH, "1")
    gv = lv.GenomeViewer(figsize=(8, 2))
    painter.draw(gv.axes[0])
    gv.savefig("tests/output/bigwig.png")


