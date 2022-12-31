#!/usr/bin/env python
# coding: utf-8

"""Test docstring for annotation.py"""

from __future__ import annotations
from typing import (
    Callable,
    Optional,
    Union,
    Literal,
)
from collections.abc import (
    Iterable,
    Sequence,
    Mapping,
)
import collections
import os
import tempfile
from dataclasses import dataclass
import warnings
import functools
import itertools
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.collections import LineCollection
import pysam

from .helpers import filter_by_keys, sort_by_keys, pack_intervals
from .plot import get_ax_size
from ._region_notation import (
    parse_region_notation,
    normalize_region_notation,
    get_region_notation,
)
from ._custom_types import (
    Identifier,
    GroupIdentifier,
    LinkIdentifier,
    Color,
    Position,
    Axes,
    Base,
    Point,
    Line,
)


# TODO: Get query sequence by position
# TODO: Max depth marker


class TrackPainter:
    pass


class Chromosome(TrackPainter):
    pass


class CoverageDepth(TrackPainter):
    @classmethod
    def from_samtools_depth_output(self, file):
        pass


@dataclass(frozen=True)
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


class AlignmentMatch(CigarOperation):
    pass


class Insertion(CigarOperation):
    pass


class Deletion(CigarOperation):
    pass


@dataclass
class ModifiedBase:
    reference_position: int
    canonical_base: str
    modification: str
    strand: str
    probability: Optional[float] = None


class _MismatchedBase:
    @property
    def reference_position(self) -> int:
        raise NotImplementedError

    @property
    def query_position(self) -> int:
        raise NotImplementedError

    @property
    def reference_base(self) -> str:
        raise NotImplementedError

    @property
    def query_base(self) -> str:
        raise NotImplementedError


@dataclass(init=False)
class MdMismatchedBase(_MismatchedBase):
    """
    A mismatched base inferred from MD tag.
    """

    _reference_position: int
    _query_position: int
    _reference_base: Base
    _query_base: Base

    def __init__(self, reference_position, query_position, reference_base, query_base):
        self._reference_position = reference_position
        self._query_position = query_position
        self._reference_base = reference_base
        self._query_base = query_base

    @property
    def reference_position(self) -> int:
        return self._reference_position

    @property
    def query_position(self) -> int:
        return self._query_position

    @property
    def reference_base(self) -> str:
        return self._reference_base

    @property
    def query_base(self) -> str:
        return self._query_base


@dataclass(frozen=True)
class CigarMismatchedBase(_MismatchedBase):
    """
    A mismatched base inferred from CIGAR.
    """

    segment: AlignedSegment
    reference_offset: int
    query_offset: int

    @property
    def reference_position(self) -> int:
        return self.segment.reference_start + self.reference_offset

    @property
    def query_position(self) -> int:
        return self.query_offset

    @property
    def reference_base(self) -> str:
        return self.segment.reference_sequence[self.reference_offset]

    @property
    def query_base(self) -> str:
        return self.segment.query_sequence[self.query_position]


class ReferenceSkip(CigarOperation):
    pass


class ClippedBases(CigarOperation):
    pass


class SoftClippedBases(ClippedBases):
    pass


class HardClippedBases(ClippedBases):
    pass


@dataclass
class CIGAR:
    alignment_matches: list[AlignmentMatch]
    insertions: list[Insertion]
    deletions: list[Deletion]
    reference_skips: list[ReferenceSkip]
    soft_clippings: list[SoftClippedBases]
    hard_clippings: list[HardClippedBases]
    mismatched_bases: list[CigarMismatchedBase]

    @classmethod
    def from_aligned_segment(cls, segment: AlignedSegment):
        if segment.cigartuples is None:
            raise ValueError("Segment is not aligned.")
        cigartuples: list[tuple[int, int]] = segment.cigartuples
        alignment_matches = []
        insertions = []
        deletions = []
        reference_skips = []
        soft_clippings = []
        hard_clippings = []
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
                soft_clippings.append(
                    SoftClippedBases(
                        segment=segment, reference_offset=ref_offset, size=length
                    )
                )
            elif operation == 5:  # Hard clipping
                hard_clippings.append(
                    HardClippedBases(
                        segment=segment,
                        size=length,
                        reference_offset=ref_offset,
                    )
                )
            elif operation == 8:  # Mismatched bases
                mismatched_bases += [
                    CigarMismatchedBase(
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
            soft_clippings=soft_clippings,
            hard_clippings=hard_clippings,
        )


@dataclass(init=False)
class AlignedSegment:
    """
    A wrapper around pysam.AlignedSegment
    """
    wrapped: pysam.AlignedSegment
    reference_start: int
    reference_end: int

    def __init__(self, wrapped: pysam.AlignedSegment):
        """
        :param wrapped: an instance of pysam.AlignedSegment.
        """
        self.wrapped = wrapped
        "The wrapped pysam.AlignedSegment object"
        if wrapped.reference_start is None or wrapped.reference_end is None:
            raise ValueError()
        self.reference_start = wrapped.reference_start
        "0-based leftmost coordinate of the aligned reference position of the segment; alias for `wrapped.reference_start`."
        self.reference_end = wrapped.reference_end
        "0-based coordinate of one past the last aligned residue of the segment; alias for `wrapped.reference_end`."
        if wrapped.query_name is None:
            raise ValueError()
        self.query_name: str = wrapped.query_name
        self.is_forward: bool = wrapped.is_forward  # type: ignore
        self.is_reverse: bool = not self.is_forward
        self.is_proper_pair: bool = wrapped.is_proper_pair
        self.is_secondary: bool = wrapped.is_secondary
        self.is_supplementary: bool = wrapped.is_supplementary
        self.is_mapped: bool = wrapped.is_mapped  # type: ignore
        self.query_alignment_length: int = wrapped.query_alignment_length

        if self.cigartuples is not None:
            self.cigar = CIGAR.from_aligned_segment(self)

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
    def reference_position_dict(self) -> dict[int, Optional[int]]:
        return {
            qry_pos: ref_pos
            for qry_pos, ref_pos in self.get_aligned_pairs()
            if qry_pos is not None
        }

    def get_aligned_reference_position(self, query_position: int) -> Optional[int]:
        return self.reference_position_dict[query_position]

    def get_aligned_query_sequence(self, reference_position: int) -> Optional[str]:
        # Note: reference_position is 0-indexed.
        if (
            reference_position < self.reference_start
            or reference_position >= self.reference_end
        ):
            return None  # Outside the alignment
        query_positions: list[int] = []
        for qry_pos, ref_pos in self.get_aligned_pairs():
            if ref_pos == reference_position:
                if qry_pos is None:
                    return "-"  # Deletion
                query_positions.append(qry_pos)
            elif query_positions:
                break
        return "".join(self.query_sequence[i] for i in query_positions)

    def get_blocks(self) -> list[tuple[int, int]]:
        return self.wrapped.get_blocks()

    def get_normalized_blocks(self) -> list[tuple[int, int]]:
        blocks: list[tuple[int, int]] = []
        for block_start, block_end in self.get_blocks():
            if blocks:
                # Check previous block and merge if possible
                previous_start, previous_end = blocks.pop(-1)
                if block_start - previous_end <= 0:
                    block_start = previous_start
                else:
                    blocks.append((previous_start, previous_end))
            blocks.append((block_start, block_end))
        return blocks

    @property
    def alignment_matches(self) -> list[AlignmentMatch]:
        return self.cigar.alignment_matches

    @property
    def insertions(self):
        return self.cigar.insertions

    @property
    def deletions(self):
        return self.cigar.deletions

    @property
    def soft_clipping(self):
        return self.cigar.soft_clippings

    @property
    def hard_clipping(self):
        return self.cigar.hard_clippings

    def _get_md_mismatched_bases(self) -> list[MdMismatchedBase]:
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

    @functools.cached_property
    def mismatched_bases(self) -> list[_MismatchedBase]:
        if self.cigar.mismatched_bases:
            mismatched_bases = self.cigar.mismatched_bases
        elif self.has_tag("MD"):
            mismatched_bases = self._get_md_mismatched_bases()
        else:
            mismatched_bases = []
        return mismatched_bases

    @functools.cached_property
    def modified_bases(self) -> list[ModifiedBase]:
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
    def reference_skips(self) -> list[ReferenceSkip]:
        return self.cigar.reference_skips


@dataclass
class _LinkedSegment:
    segments: list[AlignedSegment]

    @property
    def reference_start(self) -> int:
        return min(seg.reference_start for seg in self.segments)

    @property
    def reference_end(self) -> int:
        return max(seg.reference_end for seg in self.segments)


class SequenceAlignment(TrackPainter):
    """
    Plot sequence alignment from BAM files.
    """
    reference_name: str
    segments: list[AlignedSegment]
    pileup_depths: dict[int, int]
    pileup_bases: dict[int, collections.Counter[str]]

    def __init__(
        self,
        reference_name: str,
        segments: list[AlignedSegment],
        pileup_depths: dict[int, int],
        pileup_bases: dict[int, collections.Counter[str]],
    ):
        self.segments: list[AlignedSegment] = [seg for seg in segments if seg.is_mapped]
        self.pileup_depths = pileup_depths
        self.pileup_bases = pileup_bases
        self.reference_name = reference_name

    @classmethod
    def _from_pysam(
        cls,
        alignment_file: pysam.AlignmentFile,
        region: str | tuple[str, tuple[int, int]] | tuple[str, None],
    ):
        # Parse region
        reference_name: str
        normalized_region: str
        if isinstance(region, str):
            reference_name, interval = parse_region_notation(region)
            normalized_region = normalize_region_notation(region)
        elif isinstance(region, tuple):
            reference_name, interval = region
            normalized_region = get_region_notation(reference_name, interval)
        else:
            raise TypeError(
                f"Invalid type for `region`: {region!r}. Expecting an instance of str | tuple[str, tuple[int, int]] | tuple[str, None]."
            )
        # Check `reference_name`
        reference_names: set[str] = set(alignment_file.references)
        if reference_name not in reference_names:
            raise ValueError(
                f"Reference name {reference_name!r} is not found. Expecting one of {reference_names!r}"
            )
        # Load segments
        segment_list: list[AlignedSegment] = [
            AlignedSegment(seg)
            for seg in alignment_file.fetch(region=normalized_region)
            if seg.is_mapped
        ]
        if not segment_list:
            raise ValueError(f"No aligned segments found in {normalized_region!r}.")
        # Load pileup
        pileup_depths: dict[int, int] = {}
        pileup_bases: dict[int, collections.Counter[str]] = {}
        pileup_column: pysam.PileupColumn
        for pileup_column in alignment_file.pileup(region=normalized_region):
            position: int = pileup_column.reference_pos
            query_bases: list[str] = [
                b.upper() for b in pileup_column.get_query_sequences() if b
            ]
            pileup_depths[position] = len(
                query_bases
            )  # col.nsegments includes reference skips; use len(query_bases) instead
            base_counter = collections.Counter(query_bases)
            if len(base_counter) > 1:
                pileup_bases[position] = base_counter
        # If a position has no coverage, there will be no columns corresponding to that position.
        # Need to manually add zeros to the pileup_depths for correct plotting
        # TODO: Performance: move this to draw_pileup()
        sorted_pileup_depths: dict[int, int] = {}
        for position in range(min(pileup_depths) - 1, max(pileup_depths) + 2):
            sorted_pileup_depths[position] = pileup_depths.get(position, 0)
        pileup_depths = sorted_pileup_depths

        return cls(
            reference_name=reference_name,
            segments=segment_list,
            pileup_depths=pileup_depths,
            pileup_bases=pileup_bases,
        )

    @classmethod
    def from_file(
        cls,
        file_path: str,  # Pysam does not support reading and writing from true python file objects. See https://pysam.readthedocs.io/en/latest/usage.html#using-streams
        region: str | tuple[str, tuple[int, int]] | tuple[str, None],
        **kw,
    ):
        """
        region: reference_name compulsory even if only one reference sequence exists for explicity
        start and end coordinates are optional
        samtools-compatible region string, or a two-element tuple containing reference name and a coordinate interval (start, end). If the coordinate interval is `None`, it will be assumed to be from the first base to the last base of the reference sequence
        """
        with pysam.AlignmentFile(
            file_path, mode="rb", require_index=True, **kw
        ) as alignment_file:
            return cls._from_pysam(alignment_file, region=region)

    @classmethod
    def from_remote(
        cls,
        url: str,
        region: str | tuple[str, tuple[int, int]] | tuple[str, None],
        *,
        index_url: str | None = None,
        **kw,
    ):
        workdir = os.getcwd()
        with tempfile.TemporaryDirectory() as d:
            try:
                os.chdir(d)  # Change workdir before index file is downloaded
                with pysam.AlignmentFile(
                    url, mode="rb", require_index=True, index_filename=index_url, **kw
                ) as alignment_file:
                    return cls._from_pysam(alignment_file, region=region)
            finally:
                os.chdir(workdir)

    @staticmethod
    def _link_segments(
        segments: Sequence[AlignedSegment], links: Sequence[LinkIdentifier]
    ) -> dict[LinkIdentifier, _LinkedSegment]:
        link_seg_dict = collections.defaultdict(list)
        for seg, link in zip(segments, links):
            link_seg_dict[link].append(seg)
        link_ls_dict = {}
        for link, seg_list in link_seg_dict.items():
            linked_segment = _LinkedSegment(seg_list)
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
            for link, offset in zip(link_ls_dict, pack_intervals(intervals))
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
    ) -> list[int]:
        # Check each LinkIdentifier only mpas to one GroupIdentifier.
        link_group_dict: dict[LinkIdentifier, GroupIdentifier] = {}
        for link, group in zip(links, groups):
            recorded_group = link_group_dict.get(link, None)
            if recorded_group is not None and recorded_group != group:
                raise ValueError(
                    f"LinkIdentifier {link!r} maps to multiple GroupIdentifier values. Expecting each LinkIdentifier values maps to only one GroupIdentifier value."
                )

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
        return list(offsets)

    def _parse_segment_parameters(
        self,
        *,
        sort_by,
        link_by,
        group_by,
        filter_by,
        color_by,
    ) -> tuple[
        list[AlignedSegment],
        list[LinkIdentifier],
        list[GroupIdentifier],
        list[Color],
    ]:
        segments = self.segments
        n_segments = len(segments)

        if segments is None:
            raise ValueError("Alignment has not been loaded.")

        # Colors
        colors: list[Color] = []
        if color_by is None:
            colors = ["lightgray"] * n_segments
        elif color_by == "proper_pair":
            colors = [
                "lightgray" if segment.is_proper_pair else "firebrick"
                for segment in segments
            ]
        elif color_by == "strand":
            colors = list(
                map(
                    lambda segment: "lightgray" if segment.is_forward else "darkgray",
                    segments,
                )
            )
        elif isinstance(color_by, str):
            raise ValueError()
        elif isinstance(color_by, Iterable):
            colors = list(color_by)
            if len(colors) != n_segments:
                raise ValueError()
        elif callable(color_by):
            colors = [color_by(seg) for seg in segments]

        # Groups
        groups: list[GroupIdentifier] = []
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
        links: list[LinkIdentifier] = []
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
        selection: list[bool] = []
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
            segments, colors, links, groups = filter_by_keys(
                segments, colors, links, groups, keys=selection
            )
            if not segments:
                warnings.warn("All segments removed after filtering.")

        # Sort segments
        keys: list[Identifier] = []
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
            segments, colors, links, groups = sort_by_keys(
                segments, colors, links, groups, keys=keys
            )

        return segments, links, groups, colors

    def _get_default_segment_height(
        self, ax: Axes, offsets, *, min_height=2, max_height=10
    ) -> float:
        _, ax_height = get_ax_size(ax)
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
            Callable[[AlignedSegment], Identifier],
            Iterable[Identifier],
            Literal["length", "start"],
            None,
        ] = None,
        link_by: Union[
            Callable[[AlignedSegment], LinkIdentifier],
            Iterable[LinkIdentifier],
            Literal["name", "pair"],  # TODO: link by pair
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
            Callable[[AlignedSegment], Color],
            Iterable[Color],
            Literal["proper_pair", "strand"],
            None,
        ] = None,
        group_labels: Union[
            Callable[[GroupIdentifier], str], Mapping[GroupIdentifier, str], None
        ] = None,
        height: Optional[float] = None,
        min_spacing: Optional[float] = None,
        show_backbones: bool = True,
        show_arrowheads: bool = True,
        show_links: bool = True,
        show_insertions: bool = True,
        min_insertion_size: float = 10,
        show_deletions: bool = True,
        min_deletion_size: float = 10,
        show_mismatches: bool = True,
        show_modified_bases: bool = False,
        show_soft_clippings: bool = True,
        min_soft_clipping_size: float = 10,
        show_hard_clippings: bool = True,
        min_hard_clipping_size: float = 10,
        show_reference_skips: bool = True,
        show_group_labels: Optional[bool] = None,
        show_group_separators: bool = True,
        max_group_height: float = 1000,
        show_overheight_markers: bool = False,
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
        group_labels_kw={},
        group_separators_kw={},
        overheight_markers_kw={},
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
        offsets = self._get_segment_offsets(
            segments,
            links,
            groups,
            max_group_offset=max_group_height,
            min_spacing=min_spacing,
        )

        # Draw overheight markers
        if show_overheight_markers:
            self._draw_overheight_markers(
                ax, segments, groups, offsets, **overheight_markers_kw
            )

        # Remove segments exceeding `max_group_offset` (overheight segments)
        segments, links, groups, offsets, colors = filter_by_keys(
            segments, links, groups, offsets, colors, keys=[y >= 0 for y in offsets]
        )

        # Get segment height
        if height is None:
            height = self._get_default_segment_height(ax, offsets)

        # Parse group labels
        unique_groups = set(groups)
        group_label_dict: dict[GroupIdentifier, str] = {}
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
        if show_soft_clippings:
            self._draw_clippings(
                ax,
                segments,
                offsets,
                operation="soft",
                height=height,
                min_clipping_size=min_soft_clipping_size,
                show_arrowheads=show_arrowheads,
                **soft_clipping_kw,
            )
        if show_hard_clippings:
            self._draw_clippings(
                ax,
                segments,
                offsets,
                operation="hard",
                height=height,
                min_clipping_size=min_hard_clipping_size,
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
        if show_group_separators is True:
            self._draw_group_separators(ax, groups, offsets, **group_separators_kw)

        # Set axis limits
        ax.set_xlim(
            min(segment.reference_start for segment in segments),
            max(segment.reference_end for segment in segments),
        )
        ax.set_ylim(max(offsets) + 1, min(offsets) - 1)
        ax.set_yticks([])

    def _draw_backbones(
        self,
        ax: Axes,
        segments: Sequence[AlignedSegment],
        offsets: Sequence[float],
        height: float,
        *,
        colors: Sequence[Color],
        **kw,
    ) -> None:
        lines: list[tuple[Point, Point]] = []
        for seg, y in zip(segments, offsets):
            for block_start, block_end in seg.get_normalized_blocks():
                start_point = (block_start - 0.5, y)
                end_point = (block_end - 0.5, y)
                lines.append((start_point, end_point))
        # Match line
        ax.add_collection(
            LineCollection(
                lines,
                linewidths=height,
                colors=colors[0] if len(set(colors)) == 1 else colors,
                zorder=0,
                facecolors="none",
            )
        )

    def _draw_arrowheads(
        self,
        ax: Axes,
        segments: Sequence[AlignedSegment],
        offsets: Sequence[float],
        height: float,
        *,
        colors: Sequence[Color],
        **kw,
    ) -> None:
        forward_xs = [seg.reference_end - 0.5 for seg in segments if seg.is_forward]
        forward_ys = [y for seg, y in zip(segments, offsets) if seg.is_forward]
        reverse_xs = [seg.reference_start - 0.5 for seg in segments if seg.is_reverse]
        reverse_ys = [y for seg, y in zip(segments, offsets) if seg.is_reverse]
        forward_marker = Path([(0, 0.5), (0.5, 0), (0, -0.5), (0, 0.5)], readonly=True)
        reverse_marker = Path([(0, 0.5), (-0.5, 0), (0, -0.5), (0, 0.5)], readonly=True)
        forward_colors = [c for seg, c in zip(segments, colors) if seg.is_forward]
        reverse_colors = [c for seg, c in zip(segments, colors) if seg.is_reverse]

        # Arrorhead marker
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
                    zorder=0.1,
                )
            else:
                ax.scatter(
                    xs,
                    ys,
                    c=marker_colors,
                    marker=marker,
                    s=height**2,
                    ec="none",
                    zorder=0.1,
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
        # Insertion marker
        ax.plot(
            xs,
            ys,
            marker=marker,
            markersize=height - linewidth,
            markeredgewidth=linewidth,
            markerfacecolor="none",
            markeredgecolor=color,
            ls="",
            zorder=0.4,
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

        # Deletion marker
        ax.plot(
            xs,
            ys,
            marker=self._deletion_marker,
            markersize=height,
            markeredgecolor=color,
            markerfacecolor="none",
            ls="",
            linewidth=linewidth,
            zorder=0.3,
            **kw,
        )
        # Deletion line
        ax.add_collection(
            LineCollection(
                lines,
                linewidths=1,
                colors="k",
                zorder=-0.1,
                facecolors="none",
            )
        )

    def _draw_reference_skips(
        self, ax, segments, offsets, *, color: Color = "#96b8c8", linewidth=1, **kw
    ):
        skip_lines: list[Line] = []
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
        colormaps: dict[
            tuple[Base, str, Literal["+", "-"]], Union[str, mpl.colors.Colormap]
        ] = {("C", "m", "+"): "Reds", ("C", "m", "-"): "Reds"},
        linewidth=1,
        zorder=0.6,
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
        group_labels: dict[GroupIdentifier, str],
        *,
        dx=0.01,  # Axes unit
        dy=0.5,  # Data unit
        size=8,
        horizontalalignment="left",
        verticalalignment="bottom",
        **kw,
    ):
        max_group_offset_dict: dict[GroupIdentifier, int] = collections.defaultdict(
            lambda: 0
        )
        for g, y in zip(groups, offsets):
            max_group_offset_dict[g] = max(max_group_offset_dict[g], y)
        for g, y in max_group_offset_dict.items():
            # TODO: add path effect for better visibility
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
        self, ax, groups, offsets, *, linewidth=1, color="gray", linestyle="-", **kw
    ):
        if len(set(groups)) <= 1:  # Only one group
            return
        max_group_offset_dict = collections.defaultdict(lambda: 0)
        for g, y in zip(groups, offsets):
            max_group_offset_dict[g] = max(max_group_offset_dict[g], y)
        for y in max_group_offset_dict.values():
            ax.axhline(
                y + 1, linewidth=linewidth, color=color, linestyle=linestyle, **kw
            )

    def _draw_links(
        self,
        ax,
        segments,
        offsets,
        links,
        *,
        linewidth=1.5,
        color="lightgray",
        linestyle="-",
        **kw,
    ):
        link_ls_dict: dict[LinkIdentifier, _LinkedSegment] = self._link_segments(
            segments, links
        )
        link_offset_dict: dict[LinkIdentifier, int] = {}
        for link, y in zip(links, offsets):
            if link in link_offset_dict:
                if link_offset_dict[link] != y:
                    raise ValueError()
            else:
                link_offset_dict[link] = y

        link_lines: list[Line] = [
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

    def _draw_clippings(
        self,
        ax,
        segments,
        offsets,
        height,
        *,
        operation: Literal["hard", "soft"],
        min_clipping_size: float,
        show_arrowheads: bool,
        linewidth=1.5,
        color: Color = "deeppink",
        **kw,
    ):
        xs = []
        ys = []
        tail_marker = Path([(0, 0.5), (0, -0.5)], readonly=True)
        fwd_head_marker = Path([(0, 0.5), (0.5, 0), (0, -0.5)], readonly=True)
        rev_head_marker = Path([(0, 0.5), (-0.5, 0), (0, -0.5)], readonly=True)

        # Group clipping by type
        xs_dict: dict[str, list[float]] = collections.defaultdict(list)
        ys_dict: dict[str, list[float]] = collections.defaultdict(list)
        for seg, y in zip(segments, offsets):
            if operation == "soft":
                clippings = seg.soft_clipping
            elif operation == "hard":
                clippings = seg.soft_clipping
            else:
                raise TypeError(f"Got {operation=}, expecting Literal['hard', 'soft']")
            for clip in clippings:
                if clip.size >= min_clipping_size:
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
        # Draw clippings
        for clip_type in set(xs_dict):
            xs = xs_dict[clip_type]
            ys = ys_dict[clip_type]
            if show_arrowheads:
                marker = dict(
                    fwd_head=fwd_head_marker,
                    fwd_tail=tail_marker,
                    rev_head=rev_head_marker,
                    rev_tail=tail_marker,
                )[clip_type]
            else:
                marker = tail_marker
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

    def _draw_overheight_markers(
        self,
        ax: Axes,
        segments: Sequence[AlignedSegment],
        groups: Sequence[GroupIdentifier],
        offsets: Sequence[int],
        linewidth: float = 3,
        color: Color = "black",
        **kw,
    ) -> None:
        # TODO: remove this?
        # Get upper limit (minimum offset) of each group
        group_min_offset_dict: dict[GroupIdentifier, int] = {}
        for group, y in zip(groups, offsets):
            if y < 0:
                continue
            if group not in group_min_offset_dict:
                group_min_offset_dict[group] = y
            else:
                group_min_offset_dict[group] = min(group_min_offset_dict[group], y)
        # Get overheight segments
        overheight_segments, overheight_groups = filter_by_keys(
            segments, groups, keys=[y < 0 for y in offsets]
        )
        # Draw overheight markers
        overheight_lines: list[Line] = []
        for segment, group in zip(overheight_segments, overheight_groups):
            y = group_min_offset_dict[group]
            x0 = segment.reference_start
            x1 = segment.reference_end
            overheight_lines.append(((x0, y), (x1, y)))
        print(overheight_lines)
        ax.add_collection(
            LineCollection(
                overheight_lines,
                linewidths=linewidth,
                colors=color,
                zorder=10,
                facecolors="none",
                **kw,
            )
        )

    @functools.cached_property
    def _reference_bases(self) -> dict[int, Optional[str]]:
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
        show_mismatches: bool = True,
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
        palette: dict[Base, Color] = {
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
