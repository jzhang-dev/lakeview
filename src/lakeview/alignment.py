#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
import collections

from typing import (
    Callable,
    Iterable,
    Optional,
    List,
    Tuple,
    Dict,
    Sequence,
    Union,
    Mapping,
    Literal,
    Collection,
)
from dataclasses import dataclass
import warnings
import functools
import itertools
import numpy as np
import matplotlib as mpl
from matplotlib.path import Path
from matplotlib.collections import LineCollection
import pysam

from . import helpers
from .custom_types import (
    NativeHashable,
    GroupIdentifier,
    LinkIdentifier,
    Color,
    Position,
    Axes,
    Base,
    Point,
    Line,
)

# TODO: label metadata
# Ref skips
# Get query sequence by position
# Max depth marker


class TrackPainter:
    pass


class Chromosome(TrackPainter):
    pass


class CoverageDepth(TrackPainter):
    @classmethod
    def from_samtools_depth_output(self, file):
        pass


@dataclass
class CigarOperation:
    segment: AlignedSegment
    reference_offset: int  # Start position relative to segment.reference_start
    size: int

    @property
    def reference_position(self) -> float:
        reference_start = self.segment.reference_start
        if reference_start is None:
            raise ValueError("Segment is not aligned.")
        return reference_start + self.reference_offset


@dataclass
class AlignmentMatch(CigarOperation):
    pass


@dataclass
class Insertion(CigarOperation):
    pass


@dataclass
class Deletion(CigarOperation):
    pass


@dataclass
class ModifiedBase:
    reference_position: int
    canonical_base: str
    modification: str
    strand: str
    probability: Optional[float] = None


@dataclass(init=False)
class MismatchedBase(CigarOperation):
    segment: AlignedSegment
    reference_offset: int
    size: int
    query_offset: int
    depth: int = 1

    def __init__(
        self,
        segment: AlignedSegment,
        reference_offset: int,
        query_offset: int,
        depth: int = 1,
    ):
        super().__init__(segment=segment, reference_offset=reference_offset, size=1)
        self.query_offset = query_offset

    @property
    def query_position(self) -> int:
        return self.query_offset

    @property
    def query_base(self) -> Base:
        return self.segment.query_sequence[self.query_position]

    @property
    def reference_base(self) -> str:
        return self.segment.reference_sequence[self.reference_offset]


@dataclass
class ReferenceSkip(CigarOperation):
    pass


@dataclass
class ClippedBases(CigarOperation):
    pass


class SoftClippedBases(ClippedBases):
    pass


class HardClippedBases(ClippedBases):
    pass


@dataclass
class CIGAR:
    alignment_matches: List[AlignmentMatch]
    insertions: List[Insertion]
    deletions: List[Deletion]
    reference_skips: List[ReferenceSkip]
    soft_clipping: List[SoftClippedBases]
    hard_clipping: List[HardClippedBases]
    mismatched_bases: List[MismatchedBase]

    @classmethod
    def from_aligned_segment(cls, segment: AlignedSegment):
        if segment.cigartuples is None:
            raise ValueError("Segment is not aligned.")
        cigartuples: List[Tuple[int, int]] = segment.cigartuples
        alignment_matches = []
        insertions = []
        deletions = []
        reference_skips = []
        soft_clipping = []
        hard_clipping = []
        mismatched_bases = []
        ref_offset = 0
        qry_offset = 0

        for operation, length in cigartuples:
            if operation == 0:  # Alignment match; may not be a sequence match
                pass
            elif operation == 1:  # Insertion
                insertions.append(
                    Insertion(
                        segment=segment,
                        reference_offset=ref_offset,
                        size=length,
                    )
                )
            elif operation == 2:  # Deletion
                deletions.append(
                    Deletion(
                        segment=segment,
                        reference_offset=ref_offset,
                        size=length,
                    )
                )
            elif operation == 3:  # Reference skip
                reference_skips.append(
                    ReferenceSkip(
                        segment=segment, reference_offset=ref_offset, size=length
                    )
                )
            elif operation == 4:  # Soft clipping
                soft_clipping.append(
                    SoftClippedBases(
                        segment=segment, reference_offset=ref_offset, size=length
                    )
                )
            elif operation == 5:  # Hard clipping
                hard_clipping.append(
                    HardClippedBases(
                        segment=segment,
                        size=length,
                        reference_offset=ref_offset,
                    )
                )
            elif operation == 8:  # Mismatched bases
                mismatched_bases += [
                    MismatchedBase(
                        segment=segment,
                        reference_offset=ref_offset + i,
                        query_offset=qry_offset + i,
                    )
                    for i in range(length)
                ]
            if operation in (0, 7, 8):  # Alignment matches
                alignment_matches.append(
                    AlignmentMatch(
                        segment=segment,
                        reference_offset=ref_offset,
                        size=length,
                    )
                )
            if operation in (0, 2, 3, 7, 8):
                # Only these operations 'consume reference'
                ref_offset += length
            if operation in (0, 1, 4, 7, 8):
                # Only these operations 'consume query'
                qry_offset += length
            if operation == 9:
                # TODO: check Operation 9 (BAM_CBACK)
                raise ValueError(f"operation={operation}")

        return cls(
            alignment_matches=alignment_matches,
            insertions=insertions,
            deletions=deletions,
            reference_skips=reference_skips,
            mismatched_bases=mismatched_bases,
            soft_clipping=soft_clipping,
            hard_clipping=hard_clipping,
        )


@dataclass
class MdMismatchedBase:  # TODO
    reference_position: int
    query_position: int
    reference_base: int
    query_base: int


# A wrapper around pysam.AlignedSegment
@dataclass(init=False)
class AlignedSegment:
    wrapped: pysam.AlignedSegment
    reference_start: int
    reference_end: int

    def __init__(self, wrapped: pysam.AlignedSegment):
        self.wrapped = wrapped
        if wrapped.reference_start is None or wrapped.reference_end is None:
            raise ValueError()
        self.reference_start = wrapped.reference_start
        self.reference_end = wrapped.reference_end
        if wrapped.query_name is None:
            raise ValueError()
        self.query_name: str = wrapped.query_name
        self.is_forward: bool = wrapped.is_forward  # type: ignore
        self.is_reverse: bool = not self.is_forward
        self.is_proper_pair: bool = wrapped.is_proper_pair
        self.is_secondary: bool = wrapped.is_secondary
        self.query_alignment_length: int = wrapped.query_alignment_length

        if self.cigartuples is not None:
            self.cigar = CIGAR.from_aligned_segment(self)

    # def __getattr__(self, name):
    #     return getattr(self.wrapped, name)

    def has_tag(self, tag: str):
        return self.wrapped.has_tag(tag)

    def get_tag(self, tag: str):
        return self.wrapped.get_tag(tag)

    @property
    def cigartuples(self):
        return self.wrapped.cigartuples

    @functools.cached_property
    def query_sequence(self):
        return self.wrapped.query_sequence

    def get_aligned_pairs(self, *args, **kw):
        return self.wrapped.get_aligned_pairs(*args, **kw)

    @functools.cached_property
    def reference_sequence(self) -> str:
        return self.wrapped.get_reference_sequence()

    @functools.cached_property
    def reference_position_dict(self) -> Dict[int, Optional[int]]:
        return {
            qry_pos: ref_pos
            for qry_pos, ref_pos in self.get_aligned_pairs()
            if qry_pos is not None
        }

    def get_aligned_reference_position(self, query_position: int) -> Optional[int]:
        return self.reference_position_dict[query_position]

    def get_aligned_sequence(self, reference_start, reference_end=None):  # TODO
        pass

    @property
    def _alignment_match_intervals(self) -> List[Tuple[float, float]]:  # TODO
        reference_start = self.reference_start - 0.5
        intervals: List[Tuple[float, float]] = []
        for m in self.alignment_matches:
            interval_start = m.reference_position - 0.5
            interval_end = interval_start + m.size
            if intervals:
                # Check previous interval and merge if possible
                previous_start, previous_end = intervals.pop(-1)
                if interval_start - previous_end <= 0:
                    interval_start = previous_start
                else:
                    intervals.append((previous_start, previous_end))
            intervals.append((interval_start, interval_end))
        return intervals

    @property
    def alignment_matches(self) -> List[AlignmentMatch]:
        return self.cigar.alignment_matches

    @property
    def insertions(self):
        return self.cigar.insertions

    @property
    def deletions(self):
        return self.cigar.deletions

    @property
    def soft_clipping(self):
        return self.cigar.soft_clipping

    @property
    def hard_clipping(self):
        return self.cigar.hard_clipping

    def _get_md_mismatched_bases(self) -> List[MdMismatchedBase]:
        query_sequence = self.query_sequence
        mismatched_bases = []
        for qry_pos, ref_pos, ref_base in self.get_aligned_pairs(
            matches_only=True, with_seq=True
        ):
            qry_base = query_sequence[qry_pos]
            if qry_base != ref_base:
                mismatched_bases.append(
                    MdMismatchedBase(
                        reference_position=ref_pos,
                        query_position=qry_pos,
                        reference_base=ref_base,
                        query_base=qry_base,
                    )
                )
        return mismatched_bases

    @property
    def mismatched_bases(self) -> List[MismatchedBase]:
        if self.cigar.mismatched_bases:
            mismatched_bases = self.cigar.mismatched_bases
        elif self.has_tag("MD"):
            mismatched_bases = self._get_md_mismatched_bases()
        else:
            mismatched_bases = []
        return mismatched_bases

    @functools.cached_property
    def modified_bases(self) -> List[ModifiedBase]:
        modified_bases = []
        for (
            (
                canonical_base,
                strand,
                modification,
            ),
            data,
        ) in self.wrapped.modified_bases.items():  # type: ignore
            strand = {0: "+", 1: "-"}[strand]
            for pos, qual in data:
                if qual == -1:
                    probability = None
                else:
                    probability = qual / 256
                reference_position = self.get_aligned_reference_position(pos)
                if reference_position is not None:
                    modified_bases.append(
                        ModifiedBase(
                            reference_position=reference_position,
                            canonical_base=canonical_base,
                            modification=modification,
                            strand=strand,
                            probability=probability,
                        )
                    )
        return modified_bases

    @property
    def reference_skips(self) -> List[ReferenceSkip]:
        return self.cigar.reference_skips


@dataclass
class LinkedSegment:
    segments: List[AlignedSegment]

    @property
    def reference_start(self) -> int:
        return min(seg.reference_start for seg in self.segments)

    @property
    def reference_end(self) -> int:
        return max(seg.reference_end for seg in self.segments)


@dataclass(repr=False)
class SequenceAlignment(TrackPainter):
    segments: List[AlignedSegment]
    pileup_depths: Optional[Dict[int, int]] = None
    pileup_bases: Optional[Dict[int, collections.Counter[str]]] = None
    reference_name: Optional[str] = None
    reference_sequence: Optional[str] = None

    def __post_init__(self):
        self.segments = [seg for seg in self.segments if seg.wrapped.is_mapped]
        for segment in self.segments:
            segment.alignment = self

    @classmethod
    def from_file(
        cls,
        file_path: str,
        mode: Optional[
            Literal["r", "w", "wh", "rb", "wb", "wbu", "wb0", "rc", "wc"]
        ] = None,
        reference_name: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        *,
        region: Optional[str] = None,
        reference_sequence: Optional[str] = None,
        load_alignment: bool = True,
        load_pileup: bool = True,
        check_sq: bool = True,
        **alignment_file_kw,
    ) -> SequenceAlignment:
        if not load_alignment and not load_pileup:
            raise ValueError("`load_alignment` and `load_pileup` cannot both be False.")
        with pysam.AlignmentFile(
            file_path, mode, check_sq=check_sq, **alignment_file_kw
        ) as alignment_file:
            # Check reference_name
            reference_name_tuple = alignment_file.references
            if start is not None or end is not None:
                if reference_name is None:
                    if len(reference_name_tuple) == 1:
                        reference_name = reference_name_tuple[0]
                    else:
                        raise ValueError(
                            f"Reference name is not provided. Valid values: {reference_name_tuple}"
                        )
            # Load segments
            if load_alignment:
                segment_list = [
                    seg
                    for seg in alignment_file.fetch(
                        contig=reference_name, start=start, stop=end, region=region
                    )
                    if seg.reference_start is not None and seg.reference_end is not None
                ]
                if not segment_list:
                    warnings.warn("No aligned segments loaded.")
            else:
                segment_list = None
            # Load pileup
            if load_pileup and segment_list:
                pileup_depths = {}
                pileup_bases = {}
                for col in alignment_file.pileup(
                    contig=reference_name, start=start, stop=end, region=region
                ):
                    position: int = col.reference_pos
                    query_bases: List[str] = [
                        b.upper() for b in col.get_query_sequences() if b
                    ]
                    pileup_depths[position] = len(
                        query_bases
                    )  # col.nsegments includes reference skips; use len(query_bases) instead
                    base_counter = collections.Counter(query_bases)
                    if len(base_counter) > 1:
                        pileup_bases[position] = base_counter
                # If a position has no coverage, there will be no columns corresponding to that position.
                # Need to manually add zeros to the pileup_depths for correct plotting
                sorted_pileup_depths = {}
                for position in range(min(pileup_depths) - 1, max(pileup_depths) + 2):
                    sorted_pileup_depths[position] = pileup_depths.get(position, 0)
                pileup_depths = sorted_pileup_depths
            else:
                pileup_depths = None
                pileup_bases = None

            return cls(
                segments=[AlignedSegment(seg) for seg in segment_list],
                pileup_depths=pileup_depths,
                pileup_bases=pileup_bases,
                reference_name=reference_name,
                reference_sequence=reference_sequence,
            )

    @staticmethod
    def _link_segments(
        segments: Sequence[AlignedSegment], links: Sequence[LinkIdentifier]
    ) -> Dict[LinkIdentifier, LinkedSegment]:
        link_seg_dict = collections.defaultdict(list)
        for seg, link in zip(segments, links):
            link_seg_dict[link].append(seg)
        link_ls_dict = {}
        for link, seg_list in link_seg_dict.items():
            linked_segment = LinkedSegment(seg_list)
            link_ls_dict[link] = linked_segment
        return link_ls_dict

    @staticmethod
    def _pack_segments(
        segments: Sequence[AlignedSegment],
        links: Sequence[LinkIdentifier],
        *,
        padding: float = 0,
    ) -> np.ndarray:
        assert len(segments) == len(links)
        # Link segments
        link_ls_dict = SequenceAlignment._link_segments(segments, links)
        # Get offset for each LinkedSegment
        intervals = []
        for link, ls in link_ls_dict.items():
            intervals.append((ls.reference_start - padding, ls.reference_end + padding))
        link_offset_dict = {
            link: offset
            for link, offset in zip(link_ls_dict, helpers.pack_intervals(intervals))
        }
        # Retrive offset for each individual segment
        segment_offsets = np.array(
            [link_offset_dict[link] for link in links], dtype=np.float32
        )
        return segment_offsets

    @staticmethod
    def _get_segment_offsets(
        segments: Sequence[AlignedSegment],
        links: Sequence[LinkIdentifier],
        groups: Sequence[GroupIdentifier],
        *,
        max_group_offset: float,
        min_spacing: float,
    ) -> np.ndarray:
        if not len(segments) == len(links) == len(groups):
            raise ValueError
        # TODO: check no segments from differenct groups are linked together
        offsets = np.zeros(len(segments))
        y = 0
        for group in list(sorted(set(groups))):
            group_indices = [i for i, g in enumerate(groups) if g == group]
            group_segments = [segments[i] for i in group_indices]
            group_links = [links[i] for i in group_indices]
            group_offsets = SequenceAlignment._pack_segments(
                group_segments, group_links, padding=min_spacing / 2
            )
            group_offsets[group_offsets > max_group_offset] = -np.inf
            offsets[group_indices] = group_offsets + y
            y = max(offsets) + 2
        return offsets

    def _parse_segment_parameters(
        self,
        *,
        sort_by: Union[
            Callable[[AlignedSegment], NativeHashable],
            Iterable[NativeHashable],
            str,
            None,
        ],
        link_by: Union[
            Callable[[AlignedSegment], LinkIdentifier],
            Iterable[LinkIdentifier],
            str,
            None,
        ],
        group_by: Union[
            Callable[[AlignedSegment], GroupIdentifier],
            Iterable[GroupIdentifier],
            Literal["haplotype", "proper_pair", "strand"],
            None,
        ],
        filter_by: Union[
            Callable[[AlignedSegment], bool],
            Iterable[bool],
            str,
            None,
        ],
        color_by: Union[
            Callable[[AlignedSegment], Color],  # TODO: define a color type
            Iterable[Color],
            str,
            None,
        ],
    ) -> Tuple[
        List[AlignedSegment],
        List[LinkIdentifier],
        List[GroupIdentifier],
        List[Color],
    ]:
        segments = self.segments
        n_segments = len(segments)

        if segments is None:
            raise ValueError("Alignment has not been loaded.")

        # Colors
        colors: List[Color] = []
        if color_by is None:
            colors = ["lightgray"] * n_segments
        elif color_by == "random":
            colors = ["lightgray"] * n_segments  # TODO: random colors
        elif color_by == "strand":
            colors = list(
                map(
                    lambda segment: "lightgray" if segment.is_forward else "darkgray",
                    segments,
                )
            )  # TODO: better colors
        elif isinstance(color_by, str):
            raise ValueError()
        elif isinstance(color_by, Iterable):
            colors = list(color_by)
            if len(colors) != n_segments:
                raise ValueError()
        elif callable(color_by):
            colors = [color_by(seg) for seg in segments]

        # Groups
        groups: List[GroupIdentifier] = []
        if group_by is None:
            groups = [0] * n_segments  # Assign all segments into the same group
        elif group_by == "haplotype":
            groups = [
                seg.get_tag("HP") if seg.has_tag("HP") else float("inf")
                for seg in segments
            ]
        elif group_by == "proper_pair":
            groups = [seg.is_proper_pair for seg in segments]
        elif group_by == "strand":
            groups = ["forward" if seg.is_forward else "reverse" for seg in segments]
        elif isinstance(group_by, str):
            raise ValueError()
        elif isinstance(group_by, Iterable):
            groups = list(group_by)
            if len(groups) != n_segments:
                raise ValueError()
        if callable(group_by):
            groups = [group_by(seg) for seg in segments]

        # Links
        links: List[LinkIdentifier] = []
        if link_by is None:
            links = list(range(n_segments))  # Do not link segments together
        elif link_by == "pair":
            links = list(range(n_segments))  # TODO
        elif link_by == "name":
            links = [seg.query_name for seg in segments]
        elif isinstance(link_by, str):
            raise ValueError()
        elif isinstance(link_by, Iterable):
            links = list(link_by)
            if len(links) != n_segments:
                raise ValueError()
        elif callable(link_by):
            links = [link_by(seg) for seg in segments]

        # Filter segments
        selection: List[bool] = []
        if filter_by == "no_secondary":
            selection = [not seg.is_secondary for seg in segments]
        elif isinstance(filter_by, str):
            raise ValueError()
        elif isinstance(filter_by, Iterable):
            selection = list(filter_by)  # type: ignore
            if len(selection) != n_segments:
                raise ValueError()
        elif callable(filter_by):
            selection = [filter_by(seg) for seg in segments]
        if filter_by is not None:
            segments, colors, links, groups = helpers.filter_by(
                segments, colors, links, groups, by=selection
            )
            if not segments:
                warnings.warn("All segments removed after filtering.")

        # Sort segments
        keys: List[NativeHashable] = []
        if sort_by == "start":
            keys = [seg.reference_start for seg in segments]
        elif sort_by == "length":
            keys = [-seg.query_alignment_length for seg in segments]
        elif isinstance(sort_by, str):
            raise ValueError()
        elif isinstance(sort_by, Iterable):
            keys = list(sort_by)
            if len(keys) != n_segments:
                raise ValueError()
        elif callable(sort_by):
            keys = [sort_by(seg) for seg in segments]
        if sort_by is not None:
            segments, colors, links, groups = helpers.sort_by(
                segments, colors, links, groups, by=keys
            )

        return segments, links, groups, colors

    def _get_default_segment_height(
        self, ax: Axes, offsets, *, min_height=2, max_height=10
    ) -> float:
        _, ax_height = helpers.get_ax_size(ax)
        height = max(
            min_height,
            ax_height / (max(offsets) - min(offsets) + 2) * 0.9 * 72,
        )
        height = min(height, max_height)
        return height

    def _get_default_spacing(self, segments: Sequence[AlignedSegment]) -> float:
        segment_lengths = [seg.query_alignment_length for seg in segments]
        spacing = float(np.median(segment_lengths) * 0.1)
        return spacing

    def draw_alignment(
        self,
        ax,
        *,
        sort_by: Union[
            Callable[[AlignedSegment], NativeHashable],
            Iterable[NativeHashable],
            str,
            None,
        ] = None,
        link_by: Union[
            Callable[[AlignedSegment], LinkIdentifier],
            Iterable[LinkIdentifier],
            str,
            None,
        ] = None,
        group_by: Union[
            Callable[[AlignedSegment], GroupIdentifier],
            Iterable[GroupIdentifier],
            Literal["haplotype", "proper_pair", "strand"],
            None,
        ] = None,
        filter_by: Union[
            Callable[[AlignedSegment], bool],
            Iterable[bool],
            str,
            None,
        ] = None,
        color_by: Union[
            Callable[[AlignedSegment], Color],  # TODO: define a color type
            Iterable[Color],
            str,
            None,
        ] = None,
        group_labels: Union[
            Callable[[GroupIdentifier], str], Mapping[GroupIdentifier, str], None
        ] = None,
        height: Optional[float] = None,
        min_spacing: Optional[float] = None,
        show_backbones=True,
        show_arrowheads=True,
        show_links=True,
        show_insertions=True,
        min_insertion_size=10,
        show_deletions=True,
        min_deletion_size=10,
        show_mismatches=True,  # TODO: show_mismatches=None -> draw if available
        show_modified_bases=False,
        show_soft_clipping=True,
        min_soft_clipping_size=10,
        show_hard_clipping=True,
        min_hard_clipping_size=10,
        show_reference_skips: bool = True,
        show_letters=False,  # TODO
        show_group_labels: Optional[bool] = None,
        show_group_separators=None,
        max_group_height=1000,
        backbones_kw={},
        arrowheads_kw={},
        links_kw={},
        insertions_kw={},
        deletions_kw={},
        mismatches_kw={},
        modified_bases_kw={},
        soft_clipping_kw={},
        hard_clipping_kw={},
        reference_skips_kw={},
        letters_kw={},
        group_labels_kw={},  # TODO
        group_separators_kw={},
    ):
        """
        Groups are ordered by group id.
        Linked reads are ordered by first read in the link.
        """
        segments = self.segments

        if segments is None:
            raise ValueError("Alignment has not been loaded.")

        # Get default spacing
        if min_spacing is None:
            min_spacing = self._get_default_spacing(segments)

        # Parse segment parameters
        segments, links, groups, colors = self._parse_segment_parameters(
            filter_by=filter_by,
            group_by=group_by,
            link_by=link_by,
            color_by=color_by,
            sort_by=sort_by,
        )

        # Get segment offsets
        offset_array = self._get_segment_offsets(
            segments,
            links,
            groups,
            max_group_offset=max_group_height,
            min_spacing=min_spacing,
        )

        # Remove segments exceeding `max_group_offset`
        segments, links, groups, offsets, colors = helpers.filter_by(
            segments, links, groups, offset_array, colors, by=offset_array >= 0
        )

        # Get segment height
        if height is None:
            height = self._get_default_segment_height(ax, offsets)

        # Parse group labels
        unique_groups = set(groups)
        group_label_dict: Dict[GroupIdentifier, str] = {}
        if group_labels is None:
            if group_by == "strand":
                group_label_dict = {
                    "forward": "Forward strand",
                    "reverse": "Reverse strand",
                }
            elif group_by == "haplotype":
                group_label_dict = {hp: f"Haplotype {hp}" for hp in groups}
                group_label_dict[float("inf")] = "Haplotype unknown"
            else:
                group_label_dict = {g: str(g) for g in unique_groups}
        elif callable(group_labels):
            group_label_dict = {g: group_labels(g) for g in unique_groups}
        elif isinstance(group_labels, Mapping):
            group_label_dict = {g: group_labels[g] for g in unique_groups}

        # Draw components
        if show_backbones:
            self._draw_backbones(
                ax,
                segments,
                offsets,
                height=height,
                colors=colors,
                **backbones_kw,
            )
        if show_arrowheads:
            self._draw_arrowheads(
                ax,
                segments,
                offsets,
                height=height,
                colors=colors,
                **arrowheads_kw,
            )
        if show_links:
            self._draw_links(ax, segments, offsets, links, **links_kw)
        if show_mismatches:
            self._draw_alignment_mismatches(
                ax, segments, offsets, height=height, **mismatches_kw
            )
        if show_insertions:
            self._draw_insertions(
                ax,
                segments,
                offsets,
                height=height,
                min_insertion_size=min_insertion_size,
                **insertions_kw,
            )
        if show_deletions:
            self._draw_deletions(
                ax,
                segments,
                offsets,
                height=height,
                min_deletion_size=min_deletion_size,
                **deletions_kw,
            )
        if show_reference_skips:
            self._draw_reference_skips(ax, segments, offsets, **reference_skips_kw)
        if show_soft_clipping:
            self._draw_soft_clipping(
                ax,
                segments,
                offsets,
                height=height,
                min_soft_clipping_size=min_soft_clipping_size,
                show_arrowheads=show_arrowheads,
                **soft_clipping_kw,
            )
        if show_hard_clipping:
            self._draw_hard_clipping(
                ax,
                segments,
                offsets,
                height=height,
                min_hard_clipping_size=min_hard_clipping_size,
                show_arrowheads=show_arrowheads,
                **hard_clipping_kw,
            )
        if show_modified_bases:
            self._draw_modified_bases(
                ax, segments, offsets, height=height, **modified_bases_kw
            )
        if show_group_labels is True or (
            show_group_labels is None and group_by is not None
        ):
            self._draw_group_labels(
                ax, groups, offsets, group_labels=group_label_dict, **group_labels_kw
            )
        if show_group_separators is True or (
            show_group_separators is None and group_by is not None
        ):
            self._draw_group_separators(ax, groups, offsets, **group_separators_kw)

        # Set axis limits
        ax.set_xlim(
            min(segment.reference_start for segment in segments),
            max(segment.reference_end for segment in segments),
        )
        ax.set_ylim(max(offsets) + 1, min(offsets) - 1)
        ax.set_yticks([])

    def _draw_backbones(self, ax, segments, offsets, height, *, colors, **kw):
        lines: List[Tuple[Point, Point]] = []
        for seg, y in zip(segments, offsets):
            for interval_start, interval_end in seg._alignment_match_intervals:
                start_point = (interval_start, y)
                end_point = (interval_end, y)
                lines.append((start_point, end_point))
        ax.add_collection(
            LineCollection(
                lines,
                linewidths=height,
                colors=colors[0] if len(set(colors)) == 1 else colors,
                zorder=0,
                facecolors="none",
            )
        )

    def _draw_arrowheads(self, ax, segments, offsets, height, *, colors, **kw):

        forward_xs = [seg.reference_end - 0.5 for seg in segments if seg.is_forward]
        forward_ys = [y for seg, y in zip(segments, offsets) if seg.is_forward]
        reverse_xs = [seg.reference_start - 0.5 for seg in segments if seg.is_reverse]
        reverse_ys = [y for seg, y in zip(segments, offsets) if seg.is_reverse]
        forward_marker = Path([(0, 0.5), (0.5, 0), (0, -0.5), (0, 0.5)], readonly=True)
        reverse_marker = Path([(0, 0.5), (-0.5, 0), (0, -0.5), (0, 0.5)], readonly=True)
        forward_colors = [c for seg, c in zip(segments, colors) if seg.is_forward]
        reverse_colors = [c for seg, c in zip(segments, colors) if seg.is_reverse]

        for xs, ys, marker, marker_colors in zip(
            (forward_xs, reverse_xs),
            (forward_ys, reverse_ys),
            (forward_marker, reverse_marker),
            (forward_colors, reverse_colors),
        ):
            if not xs:
                continue
            if len(set(forward_colors)) == 1:
                ax.plot(
                    xs,
                    ys,
                    marker=marker,
                    markersize=height,
                    color=marker_colors[0],
                    markeredgecolor="none",
                    ls="",
                    zorder=1,
                )
            else:
                ax.scatter(
                    xs,
                    ys,
                    c=marker_colors,
                    marker=marker,
                    s=height**2,
                    ec="none",
                    zorder=1,
                )

    def _draw_alignment_mismatches(
        self,
        ax,
        segments,
        offsets,
        height,
        *,
        palette={
            "A": "tab:green",
            "T": "tab:red",
            "C": "tab:blue",
            "G": "tab:brown",
            None: "k",
        },
        linewidth=1.5,
        **kw,
    ):
        xs_dict = collections.defaultdict(list)
        ys_dict = collections.defaultdict(list)
        marker = Path(
            [(0, 0.5), (0, -0.5)], readonly=True
        )  # Optimally, when zoomed-in at single-base level, the mismatche marker should extend to 1-base width. For simplicity, this is not currently supported.

        for seg, y in zip(segments, offsets):
            if seg.mismatched_bases is None:
                raise RuntimeError(
                    "Failed to obtain mismatched bases using CIGAR string or the MD tag. Please provide the reference sequence or use `show_mismatched_bases=False`."
                )
            for mb in seg.mismatched_bases:
                xs_dict[mb.query_base].append(mb.reference_position)
                ys_dict[mb.query_base].append(y)

        for query_base, color in palette.items():
            xs = xs_dict[query_base]
            ys = ys_dict[query_base]
            if xs:
                ax.plot(
                    xs,
                    ys,
                    marker=marker,
                    markersize=height,
                    markeredgecolor=color,
                    markerfacecolor="none",
                    markeredgewidth=linewidth,
                    ls="",
                    **kw,
                )

    def _draw_insertions(
        self,
        ax,
        segments,
        offsets,
        height,
        *,
        min_insertion_size,
        color="tab:purple",
        linewidth=1.5,
        **kw,
    ):
        marker = Path(
            [(-0.2, 0.5), (0.2, 0.5), (0, 0.5), (0, -0.5), (-0.2, -0.5), (0.2, -0.5)],
            [
                Path.MOVETO,
                Path.LINETO,
                Path.MOVETO,
                Path.LINETO,
                Path.MOVETO,
                Path.LINETO,
            ],
            readonly=True,
        )

        xs = []
        ys = []
        for seg, y in zip(segments, offsets):
            insertions = [i for i in seg.insertions if i.size >= min_insertion_size]
            xs += [i.reference_position + 0.5 for i in insertions]
            ys += [y] * len(insertions)

        ax.plot(
            xs,
            ys,
            marker=marker,
            markersize=height - linewidth,
            markeredgewidth=linewidth,
            markerfacecolor="none",
            markeredgecolor=color,
            ls="",
            zorder=2,
            **kw,
        )

    @functools.cached_property
    def _deletion_marker(self):
        # This marker prevents the deletion from being invisible when zoomed out.
        # The 'size' of the marker appears to be automatically normalized to the farest point from the origin.
        # Two placeholder points (-1, 0) and (1, 0) are added to make the marker 'look larger' before scaling
        marker = Path(
            [(-1, 0), (-0.1, 0), (0.1, 0), (1, 0)],
            [
                Path.MOVETO,
                Path.MOVETO,
                Path.LINETO,
                Path.MOVETO,
            ],
            readonly=True,
        )
        return marker

    def _draw_deletions(
        self,
        ax,
        segments,
        offsets,
        height,
        *,
        min_deletion_size,
        color="k",
        linewidth=1.5,
        **kw,
    ):
        lines = []
        xs = []
        ys = []
        for seg, y in zip(segments, offsets):
            deletions = [d for d in seg.deletions if d.size >= min_deletion_size]
            lines += [
                (
                    (d.reference_position - 0.5, y),
                    (d.reference_position - 0.5 + d.size, y),
                )  # The reference_position for a deletion is the first deleted base
                for d in deletions
            ]
            xs += [d.reference_position + d.size / 2 for d in deletions]
            ys += [y] * len(deletions)

        # Marker
        ax.plot(
            xs,
            ys,
            marker=self._deletion_marker,
            markersize=height,
            markeredgecolor=color,
            markerfacecolor="none",
            ls="",
            linewidth=linewidth,
            zorder=1.2,
            **kw,
        )
        # Black line
        ax.add_collection(
            LineCollection(
                lines,
                linewidths=1,
                colors="k",
                zorder=1.1,
                facecolors="none",
            )
        )

    def _draw_reference_skips(
        self, ax, segments, offsets, *, color="lightgray", linewidth=1.5, **kw
    ):
        skip_lines: List[Line] = []
        for seg, y in zip(segments, offsets):
            skip_lines += [
                (
                    (skip.reference_position - 0.5, y),
                    (skip.reference_position - 0.5 + skip.size, y),
                )
                for skip in seg.reference_skips
            ]
        ax.add_collection(
            LineCollection(
                skip_lines,
                linewidths=linewidth,
                colors=color,
                zorder=-0.2,
                facecolors="none",
                **kw,
            )
        )

    def _draw_modified_bases(
        self,
        ax,
        segments,
        offsets,
        height,
        *,
        colormaps: Dict[
            Tuple[Base, str, Literal["+", "-"]], Union[str, mpl.colors.Colormap]
        ] = {("C", "m", "+"): "Reds", ("C", "m", "-"): "Reds"},
        linewidth=1,
        zorder=3,
        **kw,
    ):
        marker = Path([(0, 0.5), (0, -0.5)], readonly=True)
        xs_dict = collections.defaultdict(list)
        ys_dict = collections.defaultdict(list)
        cs_dict = collections.defaultdict(list)

        for seg, y in zip(segments, offsets):
            for mb in seg.modified_bases:
                key = (mb.canonical_base, mb.modification, mb.strand)
                xs_dict[key].append(mb.reference_position)
                ys_dict[key].append(y)
                cs_dict[key].append(mb.probability)

        for key, cmap in colormaps.items():
            xs = xs_dict[key]
            ys = ys_dict[key]
            cs = cs_dict[key]
            if xs:
                ax.scatter(
                    xs,
                    ys,
                    s=height**2,
                    c=cs,
                    cmap=cmap,
                    vmin=0,
                    vmax=1,
                    marker=marker,
                    linewidths=linewidth,
                    zorder=zorder,
                    **kw,
                )

    def _draw_group_labels(
        self,
        ax: Axes,
        groups,
        offsets,
        group_labels: Dict[GroupIdentifier, str],
        *,
        dx=0.01,  # Axes unit
        dy=0.5,  # Data unit
        size=8,
        horizontalalignment="left",
        verticalalignment="bottom",
        **kw,
    ):
        max_group_offset_dict: Dict[GroupIdentifier, int] = collections.defaultdict(
            lambda: 0
        )
        for g, y in zip(groups, offsets):
            max_group_offset_dict[g] = max(max_group_offset_dict[g], y)
        for g, y in max_group_offset_dict.items():
            ax.text(
                x=dx,
                y=y + dy,
                s=group_labels[g],
                size=size,
                horizontalalignment=horizontalalignment,
                verticalalignment=verticalalignment,
                transform=ax.get_yaxis_transform(),
                **kw,
            )

    def _draw_group_separators(
        self, ax, groups, offsets, *, linewidth=1, color="gray", linestyle="-"
    ):
        max_group_offset_dict = collections.defaultdict(lambda: 0)
        for g, y in zip(groups, offsets):
            max_group_offset_dict[g] = max(max_group_offset_dict[g], y)
        for y in max_group_offset_dict.values():
            ax.axhline(y + 1, linewidth=linewidth, color=color, linestyle=linestyle)

    def _draw_links(
        self,
        ax,
        segments,
        offsets,
        links,
        *,
        linewidth=0.5,
        color="lightgray",
        linestyle="-",
        **kw,
    ):
        link_ls_dict: Dict[LinkIdentifier, LinkedSegment] = self._link_segments(
            segments, links
        )
        link_offset_dict: Dict[LinkIdentifier, int] = {}
        for link, y in zip(links, offsets):
            if link in link_offset_dict:
                if link_offset_dict[link] != y:
                    raise ValueError()
            else:
                link_offset_dict[link] = y

        link_lines: List[Line] = [
            (
                (ls.reference_start - 0.5, link_offset_dict[link]),
                (ls.reference_end - 0.5, link_offset_dict[link]),
            )
            for link, ls in link_ls_dict.items()
        ]
        ax.add_collection(
            LineCollection(
                link_lines,
                linewidths=linewidth,
                colors=color,
                zorder=-0.3,
                facecolors="none",
            )
        )

    def _draw_soft_clipping(
        self,
        ax,
        segments,
        offsets,
        height,
        *,
        min_soft_clipping_size,
        show_arrowheads=True,
        linewidth=1.5,
        color="deeppink",
    ):
        xs = []
        ys = []
        tail_marker = Path([(0, 0.5), (0, -0.5)], readonly=True)
        fwd_head_marker = Path([(0, 0.5), (0.5, 0), (0, -0.5)], readonly=True)
        rev_head_marker = Path([(0, 0.5), (-0.5, 0), (0, -0.5)], readonly=True)

        # Group clipping by type
        xs_dict: Dict[str, List[float]] = collections.defaultdict(list)
        ys_dict: Dict[str, List[float]] = collections.defaultdict(list)
        for seg, y in zip(segments, offsets):
            for clip in seg.soft_clipping:
                if clip.size >= min_soft_clipping_size:
                    start_offset = abs(seg.reference_start - clip.reference_position)
                    end_offset = abs(seg.reference_end - clip.reference_position)
                    if start_offset < end_offset:
                        if seg.is_forward:
                            clip_type = "fwd_tail"
                        else:
                            clip_type = "rev_head"
                    else:
                        if seg.is_forward:
                            clip_type = "fwd_head"
                        else:
                            clip_type = "rev_tail"
                    xs_dict[clip_type].append(clip.reference_position - 0.5)
                    ys_dict[clip_type].append(y)
        for clip_type in set(xs_dict):
            xs = xs_dict[clip_type]
            ys = ys_dict[clip_type]
            marker = dict(
                fwd_head=fwd_head_marker,
                fwd_tail=tail_marker,
                rev_head=rev_head_marker,
                rev_tail=tail_marker,
            )[clip_type]
            ax.plot(
                xs,
                ys,
                marker=marker if show_arrowheads else tail_marker,
                markersize=height,
                markeredgecolor=color,
                markerfacecolor="none",
                markeredgewidth=linewidth,
                ls="",
            )

    def _draw_hard_clipping(
        self,
        ax,
        segments,
        offsets,
        height,
        *,
        min_hard_clipping_size,
        show_arrowheads=True,
        linewidth=1.5,
        color="deeppink",
    ):
        xs = []
        ys = []
        tail_marker = Path([(0, 0.5), (0, -0.5)], readonly=True)
        fwd_head_marker = Path([(0, 0.5), (0.5, 0), (0, -0.5)], readonly=True)
        rev_head_marker = Path([(0, 0.5), (-0.5, 0), (0, -0.5)], readonly=True)

        # Group clipping by type
        xs_dict: Dict[str, List[float]] = collections.defaultdict(list)
        ys_dict: Dict[str, List[float]] = collections.defaultdict(list)
        for seg, y in zip(segments, offsets):
            for clip in seg.hard_clipping:
                if clip.size >= min_hard_clipping_size:
                    start_offset = abs(seg.reference_start - clip.reference_position)
                    end_offset = abs(seg.reference_end - clip.reference_position)
                    if start_offset < end_offset:
                        if seg.is_forward:
                            clip_type = "fwd_tail"
                        else:
                            clip_type = "rev_head"
                    else:
                        if seg.is_forward:
                            clip_type = "fwd_head"
                        else:
                            clip_type = "rev_tail"

                    xs_dict[clip_type].append(clip.reference_position - 0.5)
                    ys_dict[clip_type].append(y)
        for clip_type in set(xs_dict):
            xs = xs_dict[clip_type]
            ys = ys_dict[clip_type]
            marker = dict(
                fwd_head=fwd_head_marker,
                fwd_tail=tail_marker,
                rev_head=rev_head_marker,
                rev_tail=tail_marker,
            )[clip_type]
            ax.plot(
                xs,
                ys,
                marker=marker if show_arrowheads else tail_marker,
                markersize=height,
                markeredgecolor=color,
                markerfacecolor="none",
                markeredgewidth=linewidth,
                ls="",
            )

    @functools.cached_property
    def _reference_bases(self) -> Dict[int, Optional[str]]:
        reference_base_dict: collections.defaultdict = collections.defaultdict(
            lambda: None
        )
        for seg in self.segments:
            for mb in seg.mismatched_bases:
                if mb.reference_position not in reference_base_dict:
                    reference_base_dict[
                        mb.reference_position
                    ] = mb.reference_base.upper()
        return reference_base_dict

    def draw_pileup(
        self,
        ax: Axes,
        *,
        color: Color = "lightgray",
        show_mismatches: bool = True,  # TODO: show_mismatches=None -> draw if available
        min_alt_frequency: float = 0.2,
        min_alt_depth: float = 2,
        mismatch_kw={},
        **kw,
    ):
        if self.pileup_depths is None:
            raise ValueError("Pileup has not been loaded.")
        x = list(self.pileup_depths)
        y = list(self.pileup_depths.values())
        ax.fill_between(x, y1=y, y2=0, step="mid", facecolor=color, edgecolor="none")
        ax.set_ylim(bottom=0)
        if show_mismatches:
            self._draw_pileup_mismatches(
                ax,
                min_alt_frequency=min_alt_frequency,
                min_alt_depth=min_alt_depth,
                **mismatch_kw,
            )

    def _draw_pileup_mismatches(
        self,
        ax: Axes,
        *,
        min_alt_frequency: float,
        min_alt_depth: float,
        linewidth: float = 1.5,
        palette: Dict[Base, Color] = {
            "A": "tab:green",
            "T": "tab:red",
            "C": "tab:blue",
            "G": "tab:brown",
        },
    ):
        mismatch_positions = set()
        if self.pileup_bases is None or self.pileup_depths is None:
            raise ValueError("Pileup has not been loaded.")
        for position, base_counter in self.pileup_bases.items():
            reference_base = self._reference_bases[position]
            total_depth = self.pileup_depths[position]
            if reference_base is None:
                raise ValueError(
                    "Reference sequence is required for drawing pileup mismatches."
                )
            for base, depth in base_counter.items():
                if (
                    base != reference_base
                    and depth >= min_alt_depth
                    and depth / total_depth >= min_alt_frequency
                ):
                    mismatch_positions.add(position)
                    break
        mismatch_position_array = np.array(sorted(mismatch_positions), dtype=int)

        bottom = np.zeros(mismatch_position_array.shape)
        for base, color in palette.items():
            counts = np.array(
                [self.pileup_bases[p][base] for p in mismatch_position_array], dtype=int
            )
            nonzero = counts > 0
            xs = mismatch_position_array[nonzero]
            ys = counts[nonzero]
            bs = bottom[nonzero]

            ax.bar(
                xs,
                ys,
                width=1,
                linewidth=linewidth,
                facecolor=color,
                edgecolor=color,
                bottom=bs,
            )
            bottom += counts


class KmerDotPlot:
    pass


class OpticalMapAlignment:
    pass
