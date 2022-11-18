#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
import collections
from multiprocessing.sharedctypes import Value

from typing import Any, Callable, Iterable, Optional, List, Tuple, Dict, Sequence, Union
from numbers import Real
from dataclasses import dataclass, field
import warnings
import functools
import itertools
import numpy as np
from matplotlib.path import Path
from matplotlib.collections import LineCollection
import pysam

from . import helpers, util

# TODO: label metadata
# Integer mode; ref skips
# Binder
# Get query sequence by position
# Max depth marker
# TODO: aln end postions are exclusive in pysam


class TrackPainter:
    pass


class Chromosome(TrackPainter):
    pass


class CoverageDepth(TrackPainter):
    @classmethod
    def from_samtools_depth_output(self, file):
        pass


@dataclass
class Insertion:
    reference_offset: int
    size: int
    reference_position: Optional[int] = None


@dataclass
class Deletion:
    reference_offset: int
    size: int
    reference_position: Optional[int] = None


@dataclass
class AlignedBase:
    reference_offset: Optional[int] = None
    query_offset: Optional[int] = None
    reference_position: Optional[int] = None
    query_position: Optional[int] = None
    reference_base: Optional[str] = None
    query_base: Optional[str] = None
    depth: int = 1


@dataclass
class ModifiedBase:
    reference_position: int
    canonical_base: str
    modification: str
    strand: str
    probability: Optional[float] = None


@dataclass
class MismatchedBase:
    reference_offset: Optional[int] = None
    query_offset: Optional[int] = None
    reference_position: Optional[int] = None
    query_position: Optional[int] = None
    reference_base: Optional[str] = None
    query_base: Optional[str] = None
    depth: int = 1


# A wrapper around pysam.AlignedSegment
@dataclass(repr=False)
class AlignedSegment:
    wrapped: pysam.AlignedSegment
    alignment: Optional[SequenceAlignment] = None

    def __getattr__(self, name):
        return getattr(self.wrapped, name)

    @functools.cached_property
    def reference_position_dict(self):
        return {
            qry_pos: ref_pos
            for qry_pos, ref_pos in self.get_aligned_pairs()
            if qry_pos is not None
        }

    def get_aligned_reference_position(self, query_position: int) -> Optional[int]:
        return self.reference_position_dict[query_position]

    def get_query_base(self, reference_position: int) -> int:
        pass  # TODO

    @staticmethod
    def parse_cigartuples(cigartuples):
        insertions = []
        deletions = []
        mismatched_bases = []
        ref_offset = 0
        qry_offset = 0
        found_op_match = False
        for operation, length in cigartuples:
            if operation == 0:  # Alignment match; may not be a sequence match
                found_op_match = True
            elif operation == 1:  # Insertion
                insertions.append(Insertion(reference_offset=ref_offset, size=length))
            elif operation == 2:  # Deletion
                deletions.append(Deletion(reference_offset=ref_offset, size=length))
            # elif operation == 7:  # Sequence match
            #     found_op_equal = True
            elif operation == 8:  # Mismatched bases
                mismatched_bases += [
                    AlignedBase(
                        reference_offset=ref_offset + i, query_offset=qry_offset + i
                    )
                    for i in range(length)
                ]
            if operation in (0, 2, 3, 7, 8):
                # Only these operations 'consume reference'
                ref_offset += length
            if operation in (0, 1, 4, 7, 8):
                # Only these operations 'consume query'
                qry_offset += length
            if operation == 9:
                # TODO: check Operation 9 (BAM_CBACK)
                raise ValueError(f"operation={operation}")
            if found_op_match:
                # In the presence of the ambiguous "M" operation, disgard the "X" operations to be safe
                mismatched_bases = None
        return insertions, deletions, mismatched_bases

    @functools.cached_property
    def _cigar_differences(self):
        insertions, deletions, mismatched_bases = self.parse_cigartuples(
            self.cigartuples
        )
        reference_start = self.reference_start
        for i in insertions:
            i.reference_position = reference_start + i.reference_offset
        for d in deletions:
            d.reference_position = reference_start + d.reference_offset
        return insertions, deletions, mismatched_bases

    @property
    def insertions(self):
        return self._cigar_differences[0]

    @property
    def deletions(self):
        return self._cigar_differences[1]

    @functools.cached_property
    def mismatched_bases(self):
        query_sequence = self.query_sequence
        mismatched_bases = self._cigar_differences[2]
        if mismatched_bases is not None:
            reference_start = self.reference_start
            for m in mismatched_bases:
                m.reference_position = reference_start + m.reference_offset
                m.query_position = m.query_offset
                m.query_base = query_sequence[m.query_position]
        elif self.has_tag("MD"):
            # M operations present in CIGAR string instead of =/X operations.
            # Fall back to parsing the MD tag.
            mismatched_bases = []
            for qry_pos, ref_pos, ref_base in self.get_aligned_pairs(
                matches_only=True, with_seq=True
            ):
                qry_base = query_sequence[qry_pos]
                if qry_base != ref_base:
                    mismatched_bases.append(
                        AlignedBase(
                            reference_position=ref_pos,
                            query_position=qry_pos,
                            reference_base=ref_base,
                            query_base=qry_base,
                        )
                    )
        elif self.alignment and self.alignment.reference_sequence:
            # MD tag not present.
            # Fall back to checking user-supplied reference sequence.
            reference_sequence = self.alignment.reference_sequence
            for qry_pos, ref_pos in self.get_aligned_pairs(
                matches_only=True, with_seq=False
            ):
                ref_base = reference_sequence[ref_pos]
                qry_base = query_sequence[qry_pos]
                if qry_base != ref_base:
                    mismatched_bases.append(
                        AlignedBase(
                            reference_position=ref_pos,
                            query_position=qry_pos,
                            reference_base=ref_base,
                            query_base=qry_base,
                        )
                    )
        else:
            # All methods failed. Notify the user.
            raise RuntimeError(
                "Failed to obtain mismatched bases using CIGAR string or the MD tag. Please provide the reference sequence."
            )
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
        ) in self.wrapped.modified_bases.items():
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
                            reference_position=self.get_aligned_reference_position(pos),
                            canonical_base=canonical_base,
                            modification=modification,
                            strand=strand,
                            probability=probability,
                        )
                    )
        return modified_bases


@dataclass
class ClippedBases:
    reference_offset: int
    size: int
    reference_position: Optional[int] = None


class SoftClippedBases(ClippedBases):
    pass


class HardClippedBases(ClippedBases):
    pass


@dataclass
class CIGAR:
    insertions: List[Insertion]
    deletions: List[Deletion]
    soft_clipping: List[SoftClippedBases]
    hard_clipping: List[HardClippedBases]
    mismatched_bases: Optional[List[MismatchedBase]] = None

    @classmethod
    def from_cigartuples(
        cls, cigartuples: List[Tuple[str, int]], *, reference_start=None
    ):
        insertions = []
        deletions = []
        soft_clipping = []
        hard_clipping = []
        mismatched_bases = []
        ref_offset = 0
        qry_offset = 0
        found_op_match = False
        for operation, length in cigartuples:
            if operation == 0:  # Alignment match; may not be a sequence match
                found_op_match = True
            elif operation == 1:  # Insertion
                insertions.append(
                    Insertion(reference_offset=ref_offset + 0.5, size=length)
                )
            elif operation == 2:  # Deletion
                deletions.append(Deletion(reference_offset=ref_offset, size=length))
            elif operation == 4:  # Soft clipping
                soft_clipping.append(
                    SoftClippedBases(size=length, reference_offset=ref_offset + 0.5)
                )
            elif operation == 5:  # Hard clipping
                hard_clipping.append(
                    HardClippedBases(size=length, reference_offset=ref_offset + 0.5)
                )
            elif operation == 8:  # Mismatched bases
                mismatched_bases += [
                    MismatchedBase(
                        reference_offset=ref_offset + i, query_offset=qry_offset + i
                    )
                    for i in range(length)
                ]
            if operation in (0, 2, 3, 7, 8):
                # Only these operations 'consume reference'
                ref_offset += length
            if operation in (0, 1, 4, 7, 8):
                # Only these operations 'consume query'
                qry_offset += length
            if operation == 9:
                # TODO: check Operation 9 (BAM_CBACK)
                raise ValueError(f"operation={operation}")
            if found_op_match:
                # In the presence of the ambiguous "M" operation, disgard the "X" operations to be safe
                mismatched_bases = None

        if reference_start is not None:
            for x in itertools.chain(
                insertions,
                deletions,
                soft_clipping,
                hard_clipping,
                mismatched_bases or [],
            ):
                x.reference_position = reference_start + x.reference_offset

        return cls(
            insertions=insertions,
            deletions=deletions,
            mismatched_bases=mismatched_bases,
            soft_clipping=soft_clipping,
            hard_clipping=hard_clipping,
        )


# A wrapper around pysam.AlignedSegment
@dataclass
class AlignedSegment:
    wrapped: pysam.AlignedSegment
    alignment: Optional[SequenceAlignment] = None

    def __getattr__(self, name):
        return getattr(self.wrapped, name)

    @functools.cached_property
    def reference_position_dict(self):
        return {
            qry_pos: ref_pos
            for qry_pos, ref_pos in self.get_aligned_pairs()
            if qry_pos is not None
        }

    def get_aligned_reference_position(self, query_position: int) -> Optional[int]:
        return self.reference_position_dict[query_position]

    def get_aligned_sequence(self, reference_start, reference_end=None):  # TODO
        pass

    @functools.cached_property
    def cigar(self):
        return CIGAR.from_cigartuples(
            self.cigartuples, reference_start=self.reference_start
        )

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

    @functools.cached_property
    def mismatched_bases(self):
        query_sequence = self.query_sequence
        mismatched_bases = self.cigar.mismatched_bases
        if mismatched_bases is not None:
            # Use CIGAR mismatches if available
            for m in mismatched_bases:
                m.query_position = m.query_offset
                m.query_base = query_sequence[m.query_position]
        elif self.has_tag("MD"):
            # M operations present in CIGAR string instead of =/X operations.
            # Fall back to parsing the MD tag.
            mismatched_bases = []
            for qry_pos, ref_pos, ref_base in self.get_aligned_pairs(
                matches_only=True, with_seq=True
            ):
                qry_base = query_sequence[qry_pos]
                if qry_base != ref_base:
                    mismatched_bases.append(
                        MismatchedBase(
                            reference_position=ref_pos,
                            query_position=qry_pos,
                            reference_base=ref_base,
                            query_base=qry_base,
                        )
                    )
        elif self.alignment and self.alignment.reference_sequence:
            # MD tag not present.
            # Fall back to checking user-supplied reference sequence.
            reference_sequence = self.alignment.reference_sequence
            for qry_pos, ref_pos in self.get_aligned_pairs(
                matches_only=True, with_seq=False
            ):
                ref_base = reference_sequence[ref_pos]
                qry_base = query_sequence[qry_pos]
                if qry_base != ref_base:
                    mismatched_bases.append(
                        MismatchedBase(
                            reference_position=ref_pos,
                            query_position=qry_pos,
                            reference_base=ref_base,
                            query_base=qry_base,
                        )
                    )
        else:
            # All methods failed. Notify the user.
            raise RuntimeError(
                "Failed to obtain mismatched bases using CIGAR string or the MD tag. Please provide the reference sequence."
            )
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
        ) in self.wrapped.modified_bases.items():
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
                            reference_position=self.get_aligned_reference_position(pos),
                            canonical_base=canonical_base,
                            modification=modification,
                            strand=strand,
                            probability=probability,
                        )
                    )
        return modified_bases


@dataclass
class LinkedSegment:
    segments: List[AlignedSegment]

    @property
    def reference_start(self):
        return min(seg.reference_start for seg in self.segments)

    @property
    def reference_end(self):
        return max(seg.reference_end for seg in self.segments)

    # def __getattr__(self, name):
    #     if name in (
    #         "insertions",
    #         "deletions",
    #         "soft_clipping",
    #         "hard_clipping",
    #         "mismatched_bases",
    #         "modified_bases",
    #     ):
    #         seg_attr_list = [get_attr(seg, name) for seg in self.segments]
    #         return list(itertools.chain(*seg_attr_list))
    #     else:
    #         return super().__getattr__(self, name)


@dataclass(repr=False)
class SequenceAlignment(TrackPainter):
    segments: List[AlignedSegment]
    pileup_depths: Optional[Dict[int, int]] = None
    pileup_bases: Optional[Dict[int, Dict[str:int]]] = None
    reference_name: Optional[str] = None
    reference_sequence: Optional[str] = None

    def __post_init__(self):
        self.segments = [seg for seg in self.segments if seg.is_mapped]
        for segment in self.segments:
            segment.alignment = self

    @classmethod
    def from_file(
        cls,
        file_path,
        mode=None,
        reference_name=None,
        start=None,
        end=None,
        *,
        region=None,
        reference_sequence=None,
        load_alignment=True,
        load_pileup=True,
        check_sq=True,
        **alignment_file_kw,
    ):
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
            if any((reference_name, start, end, region)):
                region_kw = dict(
                    contig=reference_name, start=start, stop=end, region=region
                )
            else:
                region_kw = {}
            # Load segments
            if load_alignment:
                segment_list = list(alignment_file.fetch(**region_kw))
                if not segment_list:
                    warnings.warn(f"No aligned segments loaded.")
            else:
                segment_list = None
            # Load pileup
            if load_pileup:
                pileup_depths = {}
                pileup_bases = {}
                for col in alignment_file.pileup(**region_kw):
                    position = col.reference_pos
                    query_bases = [b.upper() for b in col.get_query_sequences() if b]
                    base_counter = collections.Counter(query_bases)
                    if len(base_counter) > 1:
                        pileup_bases[position] = base_counter
                    pileup_depths[position] = col.nsegments
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
        segments: Sequence[AlignedSegment], links: Sequence
    ) -> Dict[Any, LinkedSegment]:
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
        segments: Sequence[AlignedSegment], links: Sequence, *, padding=0
    ) -> np.array:
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
        segments, links, groups, *, max_group_offset, min_spacing
    ) -> List[int]:
        assert len(segments) == len(links) == len(groups)
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

    def _parse_alignment_parameters(
        self,
        *,
        filter_by: Optional[Union[Callable, Iterable, str]] = None,  # TODO
        group_by: Optional[Union[Callable, Iterable, str]] = None,
        link_by: Optional[Union[Callable, Iterable, str]] = None,  # TODO
        color_by: Optional[Union[Callable, Iterable, str]] = None,
        sort_by: Optional[Union[Callable, Iterable, str]] = None,
        max_group_offset: Optional[Real] = float("inf"),
        min_spacing: Optional[Real] = None,
    ):
        segments = self.segments
        n_segments = len(segments)

        if segments is None:
            raise ValueError("Alignment has not been loaded.")

        # Colors
        if color_by is None:
            colors = ["lightgray"] * n_segments
        elif color_by == "random":
            color_by = lambda segment: "lightgray"  # TODO: random colors
        elif color_by == "strand":
            color_by = (
                lambda segment: "lightgray" if segment.is_forward else "darkgray"
            )  # TODO: better colors
        elif isinstance(color_by, str):
            raise ValueError()
        elif isinstance(color_by, Iterable):
            colors = list(color_by)
        if isinstance(color_by, Callable):
            colors = [color_by(seg) for seg in segments]
        if len(colors) != n_segments:
            raise ValueError()

        # Groups
        if group_by is None:
            groups = [0] * n_segments  # Assign all segments into the same group
        elif group_by == "strand":
            group_by = lambda segment: 0 if segment.is_forward else 1
        elif isinstance(group_by, str):
            raise ValueError()
        elif isinstance(group_by, Iterable):
            groups = list(group_by)
        if isinstance(group_by, Callable):
            groups = [group_by(seg) for seg in segments]
        if len(groups) != n_segments:
            raise ValueError()

        # Links
        if link_by is None:
            links = list(range(n_segments))  # Do not link segments together
        elif link_by == "pair":
            pass  # TODO
        elif link_by == "name":
            link_by = lambda seg: seg.query_name
        elif isinstance(link_by, str):
            raise ValueError()
        elif isinstance(link_by, Iterable):
            links = list(link_by)
        if isinstance(link_by, Callable):
            links = [link_by(seg) for seg in segments]

        # Filter segments
        if filter_by is None:
            selection = [True] * n_segments  # No filtering
        elif filter_by == "no_secondary":
            filter_by = lambda seg: not seg.is_secondary
        elif isinstance(filter_by, str):
            raise ValueError()
        elif isinstance(filter_by, Iterable):
            selection = list(filter_by)
        if isinstance(filter_by, Callable):
            selection = [filter_by(seg) for seg in segments]
        # TODO: implement filtering

        # Sort segments
        if sort_by == "start":
            sort_by = lambda seg: seg.reference_start
        elif sort_by == "length":
            sort_by = lambda seg: seg.query_alignment_length
        elif isinstance(sort_by, str):
            raise ValueError()
        elif isinstance(sort_by, Iterable):
            keys = list(sort_by)
        if isinstance(sort_by, Callable):
            keys = [sort_by(seg) for seg in segments]
        if sort_by is not None:
            if keys is not None and len(keys) != n_segments:
                raise ValueError()
            segments, colors, links, groups = helpers.sort_by(
                segments, colors, links, groups, by=keys
            )

        # Get default spacing
        if min_spacing is None:
            min_spacing = self._get_default_spacing(segments)

        # Get segment offsets
        offsets = self._get_segment_offsets(
            segments,
            links,
            groups,
            max_group_offset=max_group_offset,
            min_spacing=min_spacing,
        )

        # Remove segments exceeding `max_group_offset`
        segments, colors, links, groups, offsets = helpers.filter_by(
            segments, colors, links, groups, offsets, by=offsets >= 0
        )

        return segments, links, groups, offsets, colors

    def _get_default_segment_height(self, ax, offsets, *, min_height=2, max_height=10):
        _, ax_height = helpers.get_ax_size(ax)
        height = max(
            min_height,
            ax_height / (max(offsets) - min(offsets) + 2) * 0.9 * 72,
        )
        height = min(height, max_height)
        return height

    def _get_default_spacing(self, segments) -> Real:
        segment_lengths = [seg.query_alignment_length for seg in segments]
        spacing = np.median(segment_lengths) * 0.1
        return spacing

    def draw_alignment(
        self,
        ax,
        *,
        filter_by: Optional[Union[Callable, Iterable, str]] = None,  # TODO
        group_by: Optional[Union[Callable, Iterable, str]] = None,
        group_labels: Optional[Union[Callable, Iterable]] = None,  # TODO
        link_by: Optional[Union[Callable, Iterable, str]] = None,  # TODO
        color_by: Optional[Union[Callable, Iterable, str]] = None,
        sort_by: Optional[Union[Callable, Iterable, str]] = None,
        height=None,
        min_spacing: Optional[Real] = None,
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
        show_letters=False,
        show_group_separators=True,
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
        letters_kw={},
        group_separators_kw={},
    ):
        """
        Groups are ordered by group id.
        Linked reads are ordered by first read in the link.
        """
        segments = self.segments
        n_segments = len(segments)

        if segments is None:
            raise ValueError("Alignment has not been loaded.")

        # Get runtime parameters
        segments, links, groups, offsets, colors = self._parse_alignment_parameters(
            filter_by=filter_by,
            group_by=group_by,
            link_by=link_by,
            color_by=color_by,
            sort_by=sort_by,
            max_group_offset=max_group_height - 1,
            min_spacing=min_spacing,
        )

        # Get segment height
        if height is None:
            height = self._get_default_segment_height(ax, offsets)

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
        if show_group_separators:
            self._draw_group_separators(ax, groups, offsets, **group_separators_kw)

        # Set axis limits
        ax.set_xlim(
            min(segment.reference_start for segment in segments),
            max(segment.reference_end for segment in segments),
        )
        ax.set_ylim(max(offsets) + 1, min(offsets) - 1)
        ax.set_yticks([])

    def _draw_backbones(self, ax, segments, offsets, height, *, colors, **kw):
        lines = [
            ((seg.reference_start, y), (seg.reference_end, y))
            for seg, y in zip(segments, offsets)
        ]
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

        forward_xs = [seg.reference_end for seg in segments if seg.is_forward]
        forward_ys = [y for seg, y in zip(segments, offsets) if seg.is_forward]
        reverse_xs = [seg.reference_start for seg in segments if seg.is_reverse]
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
        marker = Path([(0, 0.5), (0, -0.5)], readonly=True)

        for seg, y in zip(segments, offsets):
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
            xs += [i.reference_position for i in insertions]
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

    def _draw_deletions(
        self,
        ax,
        segments,
        offsets,
        height,
        *,
        min_deletion_size,
        color="w",
        linewidth=1.5,
        **kw,
    ):
        marker = Path([(0, 0.5), (0, -0.5)], readonly=True)
        lines = []
        xs = []
        ys = []

        for seg, y in zip(segments, offsets):
            deletions = [d for d in seg.deletions if d.size >= min_deletion_size]
            lines += [
                ((d.reference_position, y), (d.reference_position + d.size, y))
                for d in deletions
            ]
            xs += [d.reference_position + d.size / 2 for d in deletions]
            ys += [y] * len(deletions)

        ax.add_collection(
            LineCollection(
                lines,
                linewidths=height,
                colors=color,
                zorder=1,
                facecolors="none",
            )
        )
        ax.plot(
            xs,
            ys,
            marker=marker,
            markersize=height,
            markeredgecolor=color,
            markerfacecolor="none",
            ls="",
            linewidth=linewidth,
            zorder=1.1,
            **kw,
        )
        ax.add_collection(
            LineCollection(
                lines,
                linewidths=1,
                colors="k",
                zorder=1.2,
                facecolors="none",
            )
        )

    def _draw_modified_bases(
        self,
        ax,
        segments,
        offsets,
        height,
        *,
        colormaps={("C", "m", "+"): "Reds", ("C", "m", "-"): "Reds"},
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
        # Link segments
        link_ls_dict = self._link_segments(segments, links)
        link_offset_dict = {}
        for link, y in zip(links, offsets):
            if link in link_offset_dict:
                if link_offset_dict[link] != y:
                    raise ValueError()
            else:
                link_offset_dict[link] = y

        lines = [
            (
                (ls.reference_start, link_offset_dict[link]),
                (ls.reference_end, link_offset_dict[link]),
            )
            for link, ls in link_ls_dict.items()
        ]
        ax.add_collection(
            LineCollection(
                lines,
                linewidths=linewidth,
                colors=color,
                zorder=-1,
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

        xs_dict = collections.defaultdict(list)
        ys_dict = collections.defaultdict(list)

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

                    xs_dict[clip_type].append(clip.reference_position)
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

        xs_dict = collections.defaultdict(list)
        ys_dict = collections.defaultdict(list)

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

                    xs_dict[clip_type].append(clip.reference_position)
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

    def _draw_pileup_mismatches(
        self,
        ax,
        *,
        min_alt_frequency,
        min_alt_depth,
        linewidth=1.5,
        palette={
            "A": "tab:green",
            "T": "tab:red",
            "C": "tab:blue",
            "G": "tab:brown",
        },
    ):

        mismatch_positions = set()
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
        mismatch_positions = np.array(sorted(mismatch_positions), dtype=int)

        bottom = np.zeros(mismatch_positions.shape)
        for base, color in palette.items():
            counts = np.array(
                [self.pileup_bases[p][base] for p in mismatch_positions], dtype=int
            )
            nonzero = counts > 0
            xs = mismatch_positions[nonzero]
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

    @functools.cached_property
    def _reference_bases(self) -> Dict[int, str]:
        reference_base_dict = collections.defaultdict(lambda: None)
        for seg in self.segments:
            for mb in seg.mismatched_bases:
                if (
                    mb.reference_base is not None
                    and mb.reference_position not in reference_base_dict
                ):
                    reference_base_dict[
                        mb.reference_position
                    ] = mb.reference_base.upper()
        return reference_base_dict

    def draw_pileup(
        self,
        ax,
        *,
        color="lightgray",
        show_mismatches=True,  # TODO: show_mismatches=None -> draw if available
        min_alt_frequency=0.2,
        min_alt_depth=2,
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
        ax,
        *,
        min_alt_frequency,
        min_alt_depth,
        linewidth=1.5,
        palette={
            "A": "tab:green",
            "T": "tab:red",
            "C": "tab:blue",
            "G": "tab:brown",
        },
    ):

        mismatch_positions = set()
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
        mismatch_positions = np.array(sorted(mismatch_positions), dtype=int)

        bottom = np.zeros(mismatch_positions.shape)
        for base, color in palette.items():
            counts = np.array(
                [self.pileup_bases[p][base] for p in mismatch_positions], dtype=int
            )
            nonzero = counts > 0
            xs = mismatch_positions[nonzero]
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
