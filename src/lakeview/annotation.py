#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import (
    Optional,
    Union,
    Callable,
    Literal,
    TextIO,
)
from collections.abc import Iterable, Sequence, Mapping, Container
from dataclasses import dataclass, field, asdict
import warnings
import numpy as np
from matplotlib.collections import LineCollection
from ._region_string import (
    parse_region_string,
    normalize_region_string,
    get_region_string,
)
from .helpers import filter_by_keys, sort_by_keys, pack_intervals
from .plot import get_ax_size
from ._custom_types import GroupIdentifier, Color, Axes, NativeHashable


@dataclass
class AnnotationRecord:
    sequence_name: str
    source: str
    feature: str
    start: int
    end: int
    score: Optional[float] = None
    strand: Optional[str] = None
    frame: Optional[str] = None
    attributes: dict[str, str] = field(default_factory=dict)

    def __len__(self):
        return self.end - self.start + 1


class GeneRecord(AnnotationRecord):  # TODO
    pass


class TranscriptRecord(AnnotationRecord):
    transcript_id: str
    gene_name: Optional[str]

    def __init__(
        self, *args, transcript_id: str, gene_name: Optional[str] = None, **kw
    ):
        super().__init__(*args, **kw)
        self.transcript_id = transcript_id
        self.gene_name = gene_name


class ExonRecord(AnnotationRecord):
    transcript_id: str

    def __init__(self, *args, transcript_id: str, **kw):
        super().__init__(*args, **kw)
        self.transcript_id = transcript_id


class CdsRecord(AnnotationRecord):
    transcript_id: str

    def __init__(self, *args, transcript_id: str, **kw):
        super().__init__(*args, **kw)
        self.transcript_id = transcript_id


@dataclass(repr=False)
class GeneAnnotation:
    genes: list[AnnotationRecord]
    transcripts: list[TranscriptRecord]
    exons: list[ExonRecord]
    cdss: list[CdsRecord]

    @staticmethod
    def parse_attribute_string(
        attribute_string,
        *,
        field_separator,
        keyval_separator,
        multival_separator,
        quoted_values,
        trim_prefix=0,
        trim_suffix=0,
    ):
        # Ref: http://daler.github.io/gffutils/dialect.html
        attr_dict = {}
        if trim_prefix > 0:
            attribute_string = attribute_string[trim_prefix:]
        if trim_suffix > 0:
            attribute_string = attribute_string[:-trim_suffix]
        for field in attribute_string.split(field_separator):
            key, value = field.split(keyval_separator)
            if quoted_values:
                value = value.strip('"')
            items = value.split(multival_separator)
            if len(items) > 1:
                value = items
            attr_dict[key] = value
        return attr_dict

    @staticmethod
    def _parse_file(
        file_object: TextIO,
        region: str | tuple[str, tuple[int, int]] | tuple[str, None],
        *,
        format_: Literal["gtf", "gff", "gff3"],
        features: Container[str],
    ) -> list[AnnotationRecord]:
        # Parse region
        sequence_name: str
        interval: tuple[int, int] | None
        normalized_region: str
        if isinstance(region, str):
            sequence_name, interval = parse_region_string(region)
            normalized_region = normalize_region_string(region)
        elif isinstance(region, tuple):
            sequence_name, interval = region
            normalized_region = get_region_string(sequence_name, interval)
        else:
            raise TypeError(
                f"Invalid type for `region`: {region!r}. Expecting an instance of str | tuple[str, tuple[int, int]] | tuple[str, None]."
            )
        start: float
        end: float
        if interval is None:
            start, end = 0, float("inf")
        else:
            start, end = interval
        # Ref: http://daler.github.io/gffutils/dialect.html
        # Ref: https://mblab.wustl.edu/GTF22.html
        # Ref: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        if format_ == "gtf":
            attr_kw = dict(
                field_separator="; ",
                keyval_separator=" ",
                multival_separator=",",
                quoted_values=True,
                trim_prefix=0,
                trim_suffix=1,
            )
        elif format_ in ("gff", "gff3"):
            attr_kw = dict(
                field_separator=";",
                keyval_separator="=",
                multival_separator=",",
                quoted_values=True,
                trim_prefix=0,
                trim_suffix=0,
            )
        else:
            raise ValueError("Only GTF and GFF3 formats are supported.")

        records: list[AnnotationRecord] = []
        for line in file_object:
            if line.startswith("#"):
                continue
            data = line.strip("\n").split("\t")
            seqname = data[0]
            if sequence_name is not None and seqname != sequence_name:
                continue
            source = data[1]
            feature = data[2]
            if features is not None and feature not in features:
                continue
            feature_start = int(data[3])
            feature_end = int(data[4])
            if start is not None and feature_end < start:
                continue
            if end is not None and feature_start > end:
                continue
            score = None if data[5] == "." else float(data[5])
            strand = None if data[6] == "." else data[6]
            frame = None if data[7] == "." else data[7]
            attr_string = data[8]
            attr_dict = GeneAnnotation.parse_attribute_string(attr_string, **attr_kw)
            records.append(
                AnnotationRecord(
                    feature=feature,
                    source=source,
                    sequence_name=seqname,
                    start=feature_start,
                    end=feature_end,
                    strand=strand,
                    score=score,
                    frame=frame,
                    attributes=attr_dict,
                )
            )
        return records

    @classmethod
    def from_file(
        cls,
        file: str | TextIO,
        region: str | tuple[str, tuple[int, int]] | tuple[str, None],
        *,
        format_: Literal["gtf", "gff", "gff3"],
        sequence_name: Optional[str] = None,
        start: Optional[float] = None,
        end: Optional[float] = None,
        gene_features: Iterable[str] = ["gene"],
        transcript_features: Iterable[str] = [
            "transcript",
            "primary_transcript",
            "RNA",
        ],
        exon_features: Iterable[str] = ["exon"],
        cds_features: Iterable[str] = ["CDS"],
        transcript_key: str = "transcript_id",  # Fetch the transcript_id from a transcript record
        parent_transcript_key: str = "transcript_id",  # Fetch the parent transcript_id from a exon/CDS record
        gene_key: str = "gene_id",  # Fetch the gene_id from a gene record. Reserved for future use
        parent_gene_key: str = "gene_id",  # Fetch the parent gene_id from a transcript record
    ):
        features: set[str] = (
            set(gene_features)
            | set(transcript_features)
            | set(exon_features)
            | set(cds_features)
        )
        if isinstance(file, str):
            with open(file, "rt") as f:
                records = cls._parse_file(
                    f,
                    region=region,
                    format_=format_,
                    features=features,
                )
        else:
            records = cls._parse_file(
                file,
                region=region,
                format_=format_,
                features=features,
            )
        if not records:
            warnings.warn("No annotation records have been loaded.")

        genes, transcripts, exons, cdss = [], [], [], []
        for record in records:
            if record.feature in gene_features:
                genes.append(record)
            elif record.feature in transcript_features:
                transcript_id = record.attributes.get(transcript_key)
                if transcript_id is None:
                    raise ValueError(f"Invalid `transcript_key` {transcript_key!r}.")
                transcript_record = TranscriptRecord(
                    **asdict(record), transcript_id=transcript_id
                )
                transcripts.append(transcript_record)
            elif record.feature in exon_features:
                transcript_id = record.attributes.get(parent_transcript_key)
                if transcript_id is None:
                    raise ValueError(
                        f"Invalid `parent_transcript_key` {transcript_key!r}."
                    )
                exon_record = ExonRecord(**asdict(record), transcript_id=transcript_id)
                exons.append(exon_record)
            elif record.feature in cds_features:
                transcript_id = record.attributes.get(parent_transcript_key)
                if transcript_id is None:
                    raise ValueError(
                        f"Invalid `parent_transcript_key` {transcript_key!r}."
                    )
                cds_record = CdsRecord(**asdict(record), transcript_id=transcript_id)
                cdss.append(cds_record)
        return cls(genes, transcripts, exons, cdss)

    @classmethod
    def from_gencode_gtf(
        cls,
        file: str | TextIO,
        region: str | tuple[str, tuple[int, int]] | tuple[str, None],
        *,
        build: Literal["GRCh38"] = "GRCh38",  # TODO: check GRCh37
    ):
        if isinstance(region, str):
            sequence_name = region
        elif isinstance(region, tuple):
            sequence_name = region[0]
        else:
            raise TypeError(
                f"Invalid type for `region`: {type(region)!r}. Expecting str | tuple[str, tuple[int, int]] | tuple[str, None]."
            )
        gene_features = ["gene"]
        transcript_features = ["transcript"]
        exon_features = ["exon"]
        cds_features = ["CDS"]
        transcript_key = "transcript_id"
        parent_transcript_key = "transcript_id"
        gene_key = "gene_id"
        parent_gene_key = "gene_id"
        instance = cls.from_file(
            file=file,
            region=region,
            format_="gtf",
            gene_features=gene_features,
            transcript_features=transcript_features,
            exon_features=exon_features,
            cds_features=cds_features,
            transcript_key=transcript_key,
            parent_transcript_key=parent_transcript_key,
            gene_key=gene_key,
            parent_gene_key=parent_gene_key,
        )
        # Parse gene names
        for transcript in instance.transcripts:
            transcript.gene_name = transcript.attributes.get("gene_name")
        return instance

    @classmethod
    def from_refseq_gff(
        cls,
        file: str | TextIO,
        region: str | tuple[str, tuple[int, int]] | tuple[str, None],
        *,
        build: Literal["GRCh37", "GRCh38"] = "GRCh38",  # TODO: check GRCh38
        chromosome: Optional[str] = None,
        start: Optional[float] = None,
        end: Optional[float] = None,
    ):
        # Sequence name
        sequence_name: str
        interval: tuple[int, int] | None
        if isinstance(region, str):
            sequence_name, interval = parse_region_string(region)
        elif isinstance(region, tuple):
            sequence_name, interval = region
        else:
            raise TypeError(
                f"Invalid type for `region`: {type(region)!r}. Expecting str | tuple[str, tuple[int, int]] | tuple[str, None]."
            )
        chromosome_dict: dict[str, dict[str, str]] = dict(
            GRCh37=dict(
                chr1="NC_000001.10",
                chr2="NC_000002.11",
                chr3="NC_000003.11",
                chr4="NC_000004.11",
                chr5="NC_000005.9",
                chr6="NC_000006.11",
                chr7="NC_000007.13",
                chr8="NC_000008.10",
                chr9="NC_000009.11",
                chr10="NC_000010.10",
                chr11="NC_000011.9",
                chr12="NC_000012.11",
                chr13="NC_000013.10",
                chr14="NC_000014.8",
                chr15="NC_000015.9",
                chr16="NC_000016.9",
                chr17="NC_000017.10",
                chr18="NC_000018.9",
                chr19="NC_000019.9",
                chr20="NC_000020.10",
                chr21="NC_000021.8",
                chr22="NC_000022.10",
                chrX="NC_000023.10",
                chrY="NC_000024.9",
                chrM="NC_012920.1",
            ),
            GRCh38=dict(),
        )
        sequence_name = chromosome_dict[build][sequence_name]
        region_string = get_region_string(sequence_name, interval)
        # Features
        gene_features = [
            "gene",
            "C_gene_segment",
            "D_gene_segment",
            "J_gene_segment",
            "V_gene_segment",
            "pseudogene",
        ]
        transcript_features = ["transcript", "mRNA"]
        exon_features = ["exon"]
        cds_features = ["CDS"]
        transcript_key = "ID"
        parent_transcript_key = "Parent"
        gene_key = "ID"
        parent_gene_key = "Parent"
        instance = cls.from_file(
            file=file,
            region=region_string,
            format_="gff",
            gene_features=gene_features,
            transcript_features=transcript_features,
            exon_features=exon_features,
            cds_features=cds_features,
            transcript_key=transcript_key,
            parent_transcript_key=parent_transcript_key,
            gene_key=gene_key,
            parent_gene_key=parent_gene_key,
        )
        # Parse gene names
        for transcript in instance.transcripts:
            transcript.gene_name = transcript.attributes.get("gene")
        return instance

    def _parse_runtime_parameters(  # TODO: delete?
        self,
        *,
        color_by,
        label_by,
    ):
        genes = self.genes
        # Colors
        if color_by is None:
            colors = ["b"] * len(genes)
        elif isinstance(color_by, Iterable):
            colors = list(color_by)
        elif isinstance(color_by, Callable):
            colors = [color_by(g) for g in genes]
        # Labels
        if label_by == None:
            labels = [""] * len(genes)
        elif isinstance(label_by, Iterable):
            labels = list(label_by)
        elif isinstance(label_by, Callable):
            labels = [labels(g) for g in genes]
        return labels, colors

    def draw_genes(
        self,
        ax: Axes,
        *,
        allow_overlaps=False,
        group_by=None,  # TODO
        group_labels=None,  # TODO
        color_by: Union[
            Callable[[AnnotationRecord], Color],  # TODO: define a color type
            Iterable[Color],
            None,
        ] = None,
        sort_by=None,  # TODO
        label_by: Union[
            Callable[[AnnotationRecord], str],
            Iterable[str],
            None,
        ] = None,  # TODO: support label_by="name"
        gene_height=None,
        show_labels=True,
        labels_kw={},
    ):
        genes = self.genes
        intervals = [(g.start, g.end) for g in genes]
        if allow_overlaps:
            offsets = [0] * len(genes)
        else:
            offsets = pack_intervals(intervals)
        colors, labels = self._parse_runtime_parameters(
            color_by=color_by, label_by=label_by
        )
        self._draw_gene_blocks(ax, genes, offsets, height=5, colors=colors)
        if show_labels:
            self._draw_labels(ax, genes, labels, offsets, colors=colors, **labels_kw)
        ax.set_ylim(max(offsets) + 0.5, min(offsets) - 0.5)
        ax.set_ylabel("")

    def _draw_gene_blocks(self, ax, genes, offsets, height, *, colors, **kw):
        lines = [((g.start, y), (g.end, y)) for g, y in zip(genes, offsets)]
        ax.add_collection(
            LineCollection(
                lines,
                linewidths=height,
                colors=colors[0] if len(set(colors)) == 1 else colors,
                zorder=0,
                facecolors="none",
            )
        )

    def _draw_labels(self, ax, annotations, labels, offsets, *, colors, size=8, **kw):
        for g, l, y, c in zip(annotations, labels, offsets, colors):
            if l:
                x = (g.start + g.end) / 2
                ax.text(
                    x, y + 0.5, l, ha="center", va="bottom", size=size, clip_on=True
                )

    @staticmethod
    def _pack_transcripts(transcripts: Sequence[AnnotationRecord]) -> np.ndarray:
        intervals = [(t.start, t.end) for t in transcripts]
        return np.array(pack_intervals(intervals), dtype=np.float32)

    @staticmethod
    def _get_transcript_offsets(transcripts, groups, *, max_group_offset) -> list[int]:
        offsets = np.zeros(len(transcripts))
        y = 0
        for group in list(sorted(set(groups))):
            group_indices = [i for i, g in enumerate(groups) if g == group]
            group_transcripts = [transcripts[i] for i in group_indices]
            group_offsets = GeneAnnotation._pack_transcripts(group_transcripts)
            group_offsets[group_offsets > max_group_offset] = -np.inf
            offsets[group_indices] = group_offsets + y
            y = max(offsets) + 2
        return list(offsets)

    def _parse_transcript_parameters(
        self,
        *,
        sort_by,
        group_by,
        filter_by,
        color_by,
        label_by,
    ) -> tuple[list[TranscriptRecord], list[GroupIdentifier], list[Color], list[str]]:
        transcripts = self.transcripts
        n_transcripts = len(transcripts)
        # Groups
        groups: list[GroupIdentifier] = []
        if group_by is None:
            groups = [0] * n_transcripts
        elif callable(group_by):
            groups = [group_by(t) for t in transcripts]
        elif isinstance(group_by, Iterable):
            groups = list(group_by)

        # Colors
        colors: list[Color] = []
        if color_by is None:
            colors = ["b"] * n_transcripts
        elif callable(color_by):
            colors = [color_by(t) for t in transcripts]
        elif isinstance(color_by, Iterable):
            colors = list(color_by)

        # Labels
        labels: list[str] = []
        if label_by is None:
            labels = [""] * n_transcripts
        elif label_by == "gene_name":
            labels = [
                t.gene_name if t.gene_name is not None else "" for t in transcripts
            ]
        elif isinstance(label_by, str):
            raise TypeError()
        elif callable(label_by):
            labels = [label_by(t) for t in transcripts]
        elif isinstance(label_by, Iterable):
            labels = list(label_by)

        # Filter transcripts
        selection: list[bool] = []
        if isinstance(filter_by, Iterable):
            selection = list(filter_by)
            if len(selection) != n_transcripts:
                raise ValueError()
        elif callable(filter_by):
            selection = [filter_by(t) for t in transcripts]
        elif isinstance(filter_by, Iterable):
            selection = list(filter_by)
        if filter_by is not None:
            transcripts, groups, labels = filter_by_keys(
                transcripts, groups, labels, keys=selection
            )
            if not transcripts:
                warnings.warn("All segments removed after filtering.")

        # Sort transcripts
        keys: list[NativeHashable] = []
        if sort_by is None:
            pass
        elif sort_by == "length":
            keys = [-len(t) for t in transcripts]
        elif isinstance(sort_by, str):
            raise TypeError()
        elif callable(sort_by):
            keys = [sort_by(t) for t in transcripts]
        elif isinstance(sort_by, Iterable):
            keys = list(sort_by)
        if sort_by is not None:
            transcripts, groups, labels = sort_by_keys(
                transcripts, groups, labels, keys=keys
            )

        return transcripts, groups, colors, labels

    def draw_transcripts(
        self,
        ax: Axes,
        *,
        sort_by: Union[
            Callable[[AnnotationRecord], NativeHashable],
            Iterable[NativeHashable],
            Literal["length"],
            None,
        ] = None,
        group_by: Union[
            Callable[[AnnotationRecord], GroupIdentifier],
            Iterable[GroupIdentifier],
            Literal["gene_name"],
            None,
        ] = None,
        filter_by: Union[
            Callable[[AnnotationRecord], bool],
            Iterable[bool],
            None,
        ] = None,
        color_by: Union[
            Callable[[AnnotationRecord], Color],
            Iterable[Color],
            None,
        ] = None,
        label_by: Union[
            Callable[[AnnotationRecord], str],
            Iterable[str],
            Literal["gene_name"],
            None,
        ] = "gene_name",
        group_labels: Union[
            Callable[[GroupIdentifier], str], Mapping[GroupIdentifier, str], None
        ] = None,
        height: Optional[float] = None,
        max_group_height=100,
        transcripts_kw={},
        exons_kw={},
        cdss_kw={},
        labels_kw={},
    ):
        transcripts, groups, colors, labels = self._parse_transcript_parameters(
            sort_by=sort_by,
            group_by=group_by,
            color_by=color_by,
            filter_by=filter_by,
            label_by=label_by,
        )

        offsets = self._get_transcript_offsets(
            transcripts, groups, max_group_offset=max_group_height - 1
        )

        transcripts, colors, groups, labels, offsets = filter_by_keys(
            transcripts, colors, groups, labels, offsets, keys=[y >= 0 for y in offsets]
        )

        if height is None:
            _, ax_height = get_ax_size(ax)
            height = max(
                3,
                ax_height / (max(offsets) - min(offsets) + 2) * 0.5 * 72,
            )

        ax.set_xlim(
            min(t.start for t in transcripts),
            max(t.end for t in transcripts),
        )
        ax.set_ylim(max(offsets) + 0.5, min(offsets) - 0.5)
        ax.set_yticks([])

        self._draw_transcript_backbones(
            ax, transcripts, offsets, colors=colors, **transcripts_kw
        )

        offset_dict = {t.transcript_id: y for t, y in zip(transcripts, offsets)}
        exons = [x for x in self.exons if x.transcript_id in offset_dict]
        cdss = [x for x in self.cdss if x.transcript_id in offset_dict]
        exon_offsets = [offset_dict[x.transcript_id] for x in exons]
        cds_offsets = [offset_dict[x.transcript_id] for x in cdss]
        color_dict = {t.transcript_id: c for t, c in zip(transcripts, colors)}
        exon_colors = [color_dict[x.transcript_id] for x in exons]
        cds_colors = [color_dict[x.transcript_id] for x in cdss]

        self._draw_exons(
            ax, exons, exon_offsets, height=height / 2, colors=exon_colors, **exons_kw
        )
        self._draw_cdss(
            ax, cdss, cds_offsets, height=height, colors=cds_colors, **cdss_kw
        )
        if labels:
            if "size" not in labels_kw:
                labels_kw = labels_kw.copy()
                labels_kw["size"] = max(height / 2, 5)
            self._draw_labels(
                ax, transcripts, labels, offsets, colors=colors, **labels_kw
            )

    def _draw_transcript_backbones(self, ax, transcripts, offsets, height=1, *, colors):
        lines = [((t.start, y), (t.end, y)) for t, y in zip(transcripts, offsets)]
        ax.add_collection(
            LineCollection(
                lines,
                linewidths=height,
                colors=colors[0] if len(set(colors)) == 1 else colors,
                zorder=0,
                facecolors="none",
            )
        )

    def _draw_exons(self, ax, exons, offsets, height, *, colors, **kw):
        lines = [((x.start, y), (x.end, y)) for x, y in zip(exons, offsets)]
        ax.add_collection(
            LineCollection(
                lines,
                linewidths=height,
                colors=colors[0] if len(set(colors)) == 1 else colors,
                zorder=0,
                facecolors="none",
            )
        )

    def _draw_cdss(self, ax, cdss, offsets, height, *, colors, **kw):
        lines = [((x.start, y), (x.end, y)) for x, y in zip(cdss, offsets)]
        ax.add_collection(
            LineCollection(
                lines,
                linewidths=height,
                colors=colors[0] if len(set(colors)) == 1 else colors,
                zorder=0,
                facecolors="none",
            )
        )
