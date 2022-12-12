#!/usr/bin/env python
# coding: utf-8

from typing import (
    Optional,
    Dict,
    List,
    Set,
    Union,
    Callable,
    Collection,
    Iterable,
    Sequence,
    Literal,
    Tuple,
    Mapping,
    TextIO,
)
from dataclasses import dataclass, field
import warnings
import numpy as np
from matplotlib.collections import LineCollection

from . import helpers
from .custom_types import GroupIdentifier, Color, Axes, NativeHashable


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
    attributes: Dict[str, str] = field(default_factory=dict)
    id: Optional[str] = None
    parent: Optional[str] = None

    def __len__(self):
        return self.end - self.start + 1


@dataclass(repr=False)
class GeneAnnotation:
    genes: List[AnnotationRecord]
    transcripts: List[AnnotationRecord]
    exons: List[AnnotationRecord]
    cdss: List[AnnotationRecord]

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
        *,
        format_: Literal["gtf", "gff", "gff3"],
        features: Set[str],
        sequence_name: Optional[str] = None,
        start: Optional[float] = None,
        end: Optional[float] = None,
        assume_sorted: bool = False,
    ) -> List[AnnotationRecord]:
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

        records: List[AnnotationRecord] = []
        for line in file_object:
            if line.startswith("#"):
                continue
            data = line.strip("\n").split("\t")
            seqname = data[0]
            if sequence_name is not None and seqname != sequence_name:
                if assume_sorted and records:
                    break
                else:
                    continue
            source = data[1]
            feature = data[2]
            if features is not None and feature not in features:
                continue
            feature_start = int(data[3])
            if start is not None and feature_start < start:
                if assume_sorted and records:
                    break
                else:
                    continue
            feature_end = int(data[4])
            if end is not None and feature_end > end:
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
        file_path=None,
        *,
        file_object=None,
        format_: Literal["gtf", "gff", "gff3"],
        sequence_name: Optional[str] = None,
        start: Optional[float] = None,
        end: Optional[float] = None,
        gene_features: Collection[str] = ["gene"],
        transcript_features: Collection[str] = [
            "transcript",
            "primary_transcript",
            "RNA",
        ],
        exon_features: Collection[str] = ["exon"],
        cds_features: Collection[str] = ["CDS"],
        transcript_key: str = "transcript_id",
        gene_key: str = "",  # Reserve for future use
        assume_sorted: bool = False,
    ):
        features: Set[str] = (
            set(gene_features)
            | set(transcript_features)
            | set(exon_features)
            | set(cds_features)
        )
        if file_path is not None and file_object is None:
            with open(file_path, "r") as f:
                records = cls._parse_file(
                    f,
                    format_=format_,
                    features=features,
                    sequence_name=sequence_name,
                    start=start,
                    end=end,
                    assume_sorted=assume_sorted,
                )
        elif file_object is not None and file_path is None:
            records = cls._parse_file(
                file_object,
                format_=format_,
                features=features,
                sequence_name=sequence_name,
                start=start,
                end=end,
                assume_sorted=assume_sorted,
            )
        else:
            raise ValueError("Either `file_path` or `file_object` must be provided.")
        if not records:
            warnings.warn("No annotation records have been loaded.")

        genes, transcripts, exons, cdss = [], [], [], []
        for r in records:
            if r.feature in gene_features:
                genes.append(r)
            elif r.feature in transcript_features:
                transcript_id = r.attributes.get(transcript_key)
                if transcript_id is None:
                    raise ValueError(f"Invalid `transcript_key` {transcript_key!r}.")
                r.id = transcript_id
                transcripts.append(r)
            elif r.feature in exon_features or r.feature in cds_features:
                attr_dict = r.attributes
                parent = attr_dict.get(transcript_key)
                if parent is None:
                    raise ValueError(f"Invalid `transcript_key` {transcript_key!r}.")
                r.parent = parent
                if r.feature in exon_features:
                    exons.append(r)
                else:
                    cdss.append(r)
        return cls(genes, transcripts, exons, cdss)

    @classmethod
    def from_gencode_gtf(
        cls,
        file_path=None,
        *,
        file_object=None,
        build: Literal["GRCh38"] = "GRCh38",
        chromosome: Optional[str] = None,
        start: Optional[float] = None,
        end: Optional[float] = None,
    ):
        sequence_name = chromosome  # TODO: check GRCh37
        transcript_features: Collection[str] = ["transcript", "primary_transcript"]
        exon_features: Collection[str] = ["exon"]
        cds_features: Collection[str] = ["CDS"]
        transcript_key: str = "transcript_id"
        return cls.from_file(
            file_path=file_path,
            file_object=file_object,
            format_="gtf",
            sequence_name=sequence_name,
            start=start,
            end=end,
            transcript_features=transcript_features,
            exon_features=exon_features,
            cds_features=cds_features,
            transcript_key=transcript_key,
        )
        # TODO: add auto gene name parsing 

    @staticmethod
    def _infer_transcript_key(attr_dict, *, preset):
        if preset == "gencode":
            transcript_key = "transcript_id"
        elif preset == "refseq":
            transcript_key = "Parent"
        elif preset == "auto":
            if "Parent" in attr_dict:
                transcript_key = "Parent"
            elif "transcript_id" in attr_dict:
                transcript_key = "transcript_id"
            else:
                raise ValueError("Failed to automatically detect `transcript_key`.")
        return transcript_key

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
            offsets = helpers.pack_intervals(intervals)
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
        return np.array(helpers.pack_intervals(intervals), dtype=np.float32)

    @staticmethod
    def _get_transcript_offsets(transcripts, groups, *, max_group_offset) -> List[int]:
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
        sort_by: Union[
            Callable[[AnnotationRecord], NativeHashable],
            Iterable[NativeHashable],
            None,
        ] = None,
        group_by: Union[
            Callable[[AnnotationRecord], GroupIdentifier],
            Iterable[GroupIdentifier],
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
            None,
        ] = None,
    ) -> Tuple[List[AnnotationRecord], List[GroupIdentifier], List[Color], List[str]]:
        transcripts = self.transcripts
        n_transcripts = len(transcripts)

        # Groups
        groups: List[GroupIdentifier] = []
        if group_by is None:
            groups = [0] * n_transcripts
        elif callable(group_by):
            groups = [group_by(t) for t in transcripts]
        elif isinstance(group_by, Iterable):
            groups = list(group_by)

        # Colors
        colors: List[Color] = []
        if color_by is None:
            colors = ["b"] * n_transcripts
        elif callable(color_by):
            colors = [color_by(t) for t in transcripts]
        elif isinstance(color_by, Iterable):
            colors = list(color_by)

        # Labels
        labels: List[str] = []
        if label_by is None:
            labels = [""] * n_transcripts
        if callable(label_by):
            labels = [label_by(t) for t in transcripts]
        elif isinstance(label_by, Iterable):
            labels = list(label_by)

        # Filter transcripts
        selection: List[bool] = []
        if isinstance(filter_by, Iterable):
            selection = list(filter_by)
            if len(selection) != n_transcripts:
                raise ValueError()
        elif callable(filter_by):
            selection = [filter_by(t) for t in transcripts]
        if filter_by is not None:
            transcripts, groups, labels = helpers.filter_by(
                transcripts, groups, labels, by=selection
            )
            if not transcripts:
                warnings.warn("All segments removed after filtering.")

        # Sort transcripts
        if callable(sort_by):
            keys = [sort_by(t) for t in transcripts]
        if sort_by is not None:
            transcripts, groups, labels = helpers.sort_by(
                transcripts, groups, labels, by=keys
            )

        return transcripts, groups, colors, labels

    def draw_transcripts(  # TODO: mask = None for reversible subsetting
        self,
        ax: Axes,
        *,
        sort_by: Union[
            Callable[[AnnotationRecord], NativeHashable],
            Iterable[NativeHashable],
            None,
        ] = None,
        group_by: Union[
            Callable[[AnnotationRecord], GroupIdentifier],
            Iterable[GroupIdentifier],
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
            None,
        ] = None,
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
        transcripts = self.transcripts
        n_transcripts = len(transcripts)

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
        transcripts, colors, groups, labels, offsets = helpers.filter_by(
            transcripts, colors, groups, labels, offsets, by=[y >= 0 for y in offsets]
        )

        if height is None:
            _, ax_height = helpers.get_ax_size(ax)
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

        offset_dict = {t.id: y for t, y in zip(transcripts, offsets)}
        exons = [x for x in self.exons if x.parent in offset_dict]
        cdss = [x for x in self.cdss if x.parent in offset_dict]
        exon_offsets = [offset_dict[x.parent] for x in exons]
        cds_offsets = [offset_dict[x.parent] for x in cdss]
        color_dict = {t.id: c for t, c in zip(transcripts, colors)}
        exon_colors = [color_dict[x.parent] for x in exons]
        cds_colors = [color_dict[x.parent] for x in cdss]

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
