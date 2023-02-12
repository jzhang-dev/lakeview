#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Optional, Union, Callable, Literal, TextIO, TypeVar
from collections.abc import Iterable, Sequence, Mapping, Container
from dataclasses import dataclass, field, asdict
import warnings
from math import floor
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.path import Path
from .region_notation import (
    parse_region_notation,
    normalize_region_notation,
    get_region_notation,
)
from .plot import get_ax_size
from ._layout import key_filter, key_sort, pack_intervals
from ._type_alias import GroupIdentifier, Color, Axes, Identifier


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
    genes: Sequence[AnnotationRecord]
    transcripts: Sequence[TranscriptRecord]
    exons: Sequence[ExonRecord]
    cdss: Sequence[CdsRecord]

    @staticmethod
    def parse_attribute_string(
        attribute_string,
        *,
        field_separator,
        key_value_separator,
        multiple_value_separator,
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
            key, value = field.split(key_value_separator)
            if quoted_values:
                value = value.strip('"')
            items = value.split(multiple_value_separator)
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
    ) -> Sequence[AnnotationRecord]:
        # Parse region
        sequence_name: str
        interval: tuple[int, int] | None
        if isinstance(region, str):
            sequence_name, interval = parse_region_notation(region)
        elif isinstance(region, tuple):
            sequence_name, interval = region
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

        # Parse format
        # Ref: http://daler.github.io/gffutils/dialect.html
        # Ref: https://mblab.wustl.edu/GTF22.html
        # Ref: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        if format_ == "gtf":
            attr_kw = dict(
                field_separator="; ",
                key_value_separator=" ",
                multiple_value_separator=",",
                quoted_values=True,
                trim_prefix=0,
                trim_suffix=1,
            )
        elif format_ == "gff3":
            attr_kw = dict(
                field_separator=";",
                key_value_separator="=",
                multiple_value_separator=",",
                quoted_values=True,
                trim_prefix=0,
                trim_suffix=0,
            )
        else:
            raise ValueError(
                f"Invalid value for `format_`: {format!r}. Expecting one of ('gtf', 'gff3')."
            )


        records: list[AnnotationRecord] = []
        visited_sequence_names: set[str] = set()
        for line in file_object:
            if line.startswith("#"):
                continue
            data = line.strip("\n").split("\t")
            seqname = data[0]
            visited_sequence_names.add(seqname)
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
        if len(records) == 0 and sequence_name is not None and sequence_name not in visited_sequence_names:
            raise ValueError(f"Invalid sequence name: {sequence_name!r}. Found following sequence names in the file: {visited_sequence_names!r}")
        return records

    @classmethod
    def from_file(
        cls,
        file: str | TextIO,
        format_: Literal["gtf", "gff3"],
        region: str | tuple[str, tuple[int, int]] | tuple[str, None],
        *,
        gene_features: Iterable[str] = ["gene"],
        transcript_features: Iterable[str] = ["transcript"],
        exon_features: Iterable[str] = ["exon"],
        cds_features: Iterable[str] = ["CDS"],
    ):
        # Combine features
        features: set[str] = (
            set(gene_features)
            | set(transcript_features)
            | set(exon_features)
            | set(cds_features)
        )
        # Parse records
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
        # Parse gene_key and transcript_key
        transcript_key: str  # Fetch the transcript_id from a transcript record
        parent_transcript_key: str  # Fetch the parent transcript_id from a exon/CDS record
        gene_key: str  # Fetch the gene_id from a gene record. Reserved for future use
        parent_gene_key: str  # Fetch the parent gene_id from a transcript record
        if format_ == "gtf":
            transcript_key = "transcript_id"
            parent_transcript_key = "transcript_id"
            # gene_key = "gene_id"
            # parent_gene_key = "gene_id"
        elif format_ == "gff3":
            transcript_key = "ID"
            parent_transcript_key = "Parent"
            # gene_key = "ID"
            # parent_gene_key = "Parent"
        # Identify genes, transcripts, exons, cdss
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
    def from_gencode(
        cls,
        file: str | TextIO,
        format_: Literal["gtf", "gff3"],
        region: str | tuple[str, tuple[int, int]] | tuple[str, None],
    ):
        instance = cls.from_file(
            file=file,
            region=region,
            format_=format_,
            gene_features=["gene"],
            transcript_features=["transcript"],
            exon_features=["exon"],
            cds_features=["CDS"],
        )
        # Parse gene names
        for transcript in instance.transcripts:
            transcript.gene_name = transcript.attributes.get("gene_name")
        return instance

    @classmethod
    def from_refseq(
        cls,
        file: str | TextIO,
        format_: Literal["gtf", "gff3"],
        region: str | tuple[str, tuple[int, int]] | tuple[str, None],
    ):
        instance = cls.from_file(
            file=file,
            region=region,
            format_=format_,
            gene_features=[
                "gene",
                "C_gene_segment",
                "D_gene_segment",
                "J_gene_segment",
                "V_gene_segment",
                "pseudogene",
            ],
            transcript_features=["transcript", "mRNA"],
            exon_features=["exon"],
            cds_features=["CDS"],
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
        offsets: Sequence[int]
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
    def _get_transcript_offsets(
        transcripts, groups, *, max_group_offset
    ) -> Sequence[int]:
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
    ) -> tuple[
        Sequence[TranscriptRecord],
        Sequence[GroupIdentifier],
        Sequence[Color],
        Sequence[str],
    ]:
        transcripts: Sequence[TranscriptRecord] = self.transcripts
        n_transcripts = len(transcripts)
        # Groups
        groups: Sequence[GroupIdentifier] = []
        if group_by is None:
            groups = [0] * n_transcripts
        elif callable(group_by):
            groups = [group_by(t) for t in transcripts]
        elif isinstance(group_by, Iterable):
            groups = list(group_by)

        # Colors
        colors: Sequence[Color] = []
        if color_by is None:
            colors = [(0, 0, 178 / 255)] * n_transcripts
        elif callable(color_by):
            colors = [color_by(t) for t in transcripts]
        elif isinstance(color_by, Iterable):
            colors = list(color_by)

        # Labels
        labels: Sequence[str] = []
        if label_by is None:
            labels = [""] * n_transcripts
        if isinstance(label_by, str):
            if label_by == "gene_name":
                labels = [
                    t.gene_name if t.gene_name is not None else "" for t in transcripts
                ]
            else:
                raise TypeError()
        elif callable(label_by):
            labels = [label_by(t) for t in transcripts]
        elif isinstance(label_by, Iterable):
            labels = list(label_by)

        # Filter transcripts
        filter_keys: Sequence[bool] = []
        if isinstance(filter_by, str):
            pass
        if isinstance(filter_by, Iterable):
            filter_keys = list(filter_by)
            if len(filter_keys) != n_transcripts:
                raise ValueError()
        elif callable(filter_by):
            filter_keys = [filter_by(t) for t in transcripts]
        elif isinstance(filter_by, Iterable):
            filter_keys = list(filter_by)

        # Sort transcripts
        sort_keys: Sequence[Identifier] = []
        if sort_by is None:
            pass
        if isinstance(sort_by, str):
            if sort_by == "length":
                sort_keys = [-len(t) for t in transcripts]
            else:
                raise TypeError()
        elif callable(sort_by):
            sort_keys = [sort_by(t) for t in transcripts]
        elif isinstance(sort_by, Iterable):
            sort_keys = list(sort_by)

        if filter_by is not None:
            transcripts = key_filter(transcripts, filter_keys)
            groups = key_filter(groups, filter_keys)
            labels = key_filter(labels, filter_keys)
            if sort_by is not None:
                sort_keys = key_filter(sort_keys, filter_keys)
            if len(transcripts) == 0:
                warnings.warn("All segments removed after filtering.")

        if sort_by is not None:
            transcripts = key_sort(transcripts, sort_keys)
            groups = key_sort(groups, sort_keys)
            labels = key_sort(labels, sort_keys)

        return transcripts, groups, colors, labels

    def draw_transcripts(
        self,
        ax: Axes,
        *,
        sort_by: Union[
            Callable[[AnnotationRecord], Identifier],
            Iterable[Identifier],
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
        max_rows: int = 100,
        draw_arrows: bool = True,
        transcripts_kw={},
        arrows_kw={},
        exons_kw={},
        cdss_kw={},
        labels_kw={},
    ):
        """
        Draw sequence alignment patterns, in a style similar to `IGH feature track <https://software.broadinstitute.org/software/igv/feature_track_options>`_.

        :param max_rows: The maximum number of rows to layout transcripts. Excess transcripts will not be drawn. If multiple transcript groups exist, this parameter limits the maximum number of rows *per group*.
        """
        transcripts, groups, colors, labels = self._parse_transcript_parameters(
            sort_by=sort_by,
            group_by=group_by,
            color_by=color_by,
            filter_by=filter_by,
            label_by=label_by,
        )

        max_group_offset = max_rows - 1
        offsets = self._get_transcript_offsets(
            transcripts, groups, max_group_offset=max_group_offset
        )

        # Remove transcripts exceeding `max_rows`
        if not all(y >= 0 for y in offsets):
            filter_keys = [y >= 0 for y in offsets]
            transcripts = key_filter(transcripts, filter_keys)
            colors = key_filter(colors, filter_keys)
            groups = key_filter(groups, filter_keys)
            offsets = key_filter(offsets, filter_keys)
            labels = key_filter(labels, filter_keys)

        # Get transcript height
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
        if draw_arrows:
            self._draw_arrows(
                ax,
                transcripts,
                offsets,
                height=height,
                colors=colors,
                **arrows_kw,
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
            ax, exons, exon_offsets, height=height, colors=exon_colors, **exons_kw
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

    def _draw_transcript_backbones(
        self, ax, transcripts, offsets, colors, *, height=1, **kw
    ):
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

    def _draw_exons(self, ax, exons, offsets, height, colors, **kw):
        line_width = height * 0.5
        lines = [((x.start, y), (x.end, y)) for x, y in zip(exons, offsets)]
        ax.add_collection(
            LineCollection(
                lines,
                linewidths=line_width,
                colors=colors[0] if len(set(colors)) == 1 else colors,
                zorder=0.1,
                facecolors="none",
            )
        )

    def _draw_cdss(self, ax, cdss, offsets, height, colors, **kw):
        line_width = height
        lines = [((x.start, y), (x.end, y)) for x, y in zip(cdss, offsets)]
        ax.add_collection(
            LineCollection(
                lines,
                linewidths=line_width,
                colors=colors[0] if len(set(colors)) == 1 else colors,
                zorder=0.2,
                facecolors="none",
            )
        )

    def _draw_arrows(
        self,
        ax: Axes,
        transcripts: Sequence[TranscriptRecord],
        offsets: Sequence[float],
        height: float,
        colors: Sequence[Color],
        *,
        style: Literal["single", "fishbone"] = "single",
        **kw,
    ) -> None:
        if style == "single":
            self._draw_single_arrows(ax, transcripts, offsets, height, colors, **kw)
        elif style == "fishbone":
            self._draw_fishbone_arrows(ax, transcripts, offsets, height, colors, **kw)
        else:
            raise ValueError()

    def _draw_fishbone_arrows(
        self,
        ax: Axes,
        transcripts: Sequence[TranscriptRecord],
        offsets: Sequence[float],
        height: float,
        colors: Sequence[Color],
        *,
        min_spacing: float = 200,
        linewidth: float = 0.5,
        **kw,
    ) -> None:
        marker_size = height * 0.5
        markers = {
            "+": Path([(0, 0.5), (0.5, 0), (0, -0.5)], readonly=True),
            "-": Path([(0, 0.5), (-0.5, 0), (0, -0.5)], readonly=True),
        }
        xs: dict[str, list[float]] = {"+": [], "-": []}
        ys: dict[str, list[float]] = {"+": [], "-": []}
        cs: dict[str, list[Color]] = {"+": [], "-": []}

        for transcript, y, color in zip(transcripts, offsets, colors):
            if transcript.strand is None:
                continue
            n_marker = max(floor((transcript.end - transcript.start) / min_spacing), 1)
            marker_positions: list[float]
            if n_marker == 1:
                marker_positions = [(transcript.end + transcript.start) / 2]
            else:
                spacing: float = (transcript.end - transcript.start) / n_marker
                marker_positions = [
                    transcript.start + (i + 0.5) * spacing for i in range(n_marker)
                ]

            xs[transcript.strand] += marker_positions
            ys[transcript.strand] += [y] * n_marker
            cs[transcript.strand] += [color] * n_marker

        for strand in ("+", "-"):
            ax.scatter(
                xs[strand],
                ys[strand],
                s=marker_size**2,
                c=cs[strand],
                marker=markers[strand],
                linewidths=linewidth,
                zorder=0.5,
                **kw,
            )

    def _draw_single_arrows(
        self,
        ax: Axes,
        transcripts: Sequence[TranscriptRecord],
        offsets: Sequence[float],
        height: float,
        colors: Sequence[Color],
        *,
        linewidth: float = 0.5,
        **kw,
    ) -> None:
        arrow_size = 0.3
        marker_size = height * 2
        fwd_marker = Path(
            [
                (0, 0),
                (0, 1),
                (1, 1),
                (1 - arrow_size, 1 - arrow_size),
                (1 - arrow_size, 1 + arrow_size),
                (1, 1),
            ],
            [
                Path.MOVETO,
                Path.LINETO,
                Path.LINETO,
                Path.LINETO,
                Path.MOVETO,
                Path.LINETO,
            ],
        )
        rev_marker = Path(
            [
                (0, 0),
                (0, -1),
                (-1, -1),
                (-1 + arrow_size, -1 - arrow_size),
                (-1 + arrow_size, -1 + arrow_size),
                (-1, -1),
            ],
            [
                Path.MOVETO,
                Path.LINETO,
                Path.LINETO,
                Path.LINETO,
                Path.MOVETO,
                Path.LINETO,
            ],
        )
        markers = {
            "+": fwd_marker,
            "-": rev_marker,
        }
        xs: dict[str, list[float]] = {"+": [], "-": []}
        ys: dict[str, list[float]] = {"+": [], "-": []}
        cs: dict[str, list[Color]] = {"+": [], "-": []}

        for transcript, y, color in zip(transcripts, offsets, colors):
            if transcript.strand is None:
                continue
            elif transcript.strand == "+":
                marker_position = transcript.start
            elif transcript.strand == "-":
                marker_position = transcript.end
            else:
                raise ValueError()
            xs[transcript.strand].append(marker_position)
            ys[transcript.strand].append(y)
            cs[transcript.strand].append(color)

        for strand in ("+", "-"):
            ax.scatter(
                xs[strand],
                ys[strand],
                s=marker_size**2,
                edgecolors=cs[strand],
                facecolors="none",
                marker=markers[strand],
                linewidths=linewidth,
                zorder=0.5,
                **kw,
            )
