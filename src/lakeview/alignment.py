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

from .plot import get_ax_size
from ._layout import key_filter, key_sort, pack_intervals
from .region_notation import (
    parse_region_notation,
    normalize_region_notation,
    get_region_notation,
)
from ._type_alias import (
    Identifier,
    GroupIdentifier,
    LinkIdentifier,
    Color,
    Axes,
    Base,
    Point,
    Line,
)


class CoverageDepth:
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
class _CIGAR:
    insertions: list[Insertion]
    deletions: list[Deletion]
    reference_skips: list[ReferenceSkip]
    soft_clippings: list[SoftClippedBases]
    hard_clippings: list[HardClippedBases]
    mismatched_bases: list[CigarMismatchedBase]

    @classmethod
    def from_aligned_segment(cls, segment: AlignedSegment):
        if segment._cigartuples is None:
            raise ValueError("Segment is not aligned.")
        cigartuples: list[tuple[int, int]] = segment._cigartuples
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
        :param wrapped: an instance of :external:py:class:`pysam.AlignedSegment`
        """
        self.wrapped = wrapped
        "the wrapped :external:py:class:`pysam.AlignedSegment` instance."
        if wrapped.reference_start is None or wrapped.reference_end is None:
            raise ValueError()
        self.reference_start = wrapped.reference_start
        "0-based coordinate of the first aligned reference position of the segment; alias of :external:py:attr:`pysam.AlignedSegment.reference_start`"
        self.reference_end = wrapped.reference_end
        "0-based coordinate of one past the last aligned reference position of the segment; alias of :external:py:attr:`pysam.AlignedSegment.reference_end`"
        if wrapped.query_name is None:
            raise ValueError()
        self.query_name: str = wrapped.query_name
        "alias of :external:py:attr:`pysam.AlignedSegment.query_name`"
        self.is_forward: bool = wrapped.is_forward  # type: ignore
        "alias of :external:py:attr:`pysam.AlignedSegment.is_forward`"
        self.is_reverse: bool = not self.is_forward
        "alias of :external:py:attr:`pysam.AlignedSegment.is_reverse`"
        self.is_proper_pair: bool = wrapped.is_proper_pair
        "alias of :external:py:attr:`pysam.AlignedSegment.is_proper_pair`"
        self.is_secondary: bool = wrapped.is_secondary
        "alias of :external:py:attr:`pysam.AlignedSegment.is_secondary`"
        self.is_supplementary: bool = wrapped.is_supplementary
        "alias of :external:py:attr:`pysam.AlignedSegment.is_secondary`"
        self.is_mapped: bool = wrapped.is_mapped  # type: ignore
        "alias of :external:py:attr:`pysam.AlignedSegment.is_mapped`"
        self.mapping_quality: int = wrapped.mapping_quality
        "alias of :external:py:attr:`pysam.AlignedSegment.mapping_quality`"
        self.query_alignment_length: int = wrapped.query_alignment_length
        "alias of :external:py:attr:`pysam.AlignedSegment.query_alignment_length`"

        if self._cigartuples is not None:
            self._cigar: _CIGAR = _CIGAR.from_aligned_segment(self)
            "alignment patterns recorded in the CIGAR string"
        else:
            raise ValueError("CIGAR not found.")

    def has_tag(self, tag: str) -> bool:
        "alias of :external:py:meth:`pysam.AlignedSegment.has_tag`"
        return self.wrapped.has_tag(tag)

    def get_tag(self, tag: str):
        "alias of :external:py:meth:`pysam.AlignedSegment.get_tag`"
        return self.wrapped.get_tag(tag)

    @property
    def _cigartuples(self):
        return self.wrapped.cigartuples

    @functools.cached_property
    def query_sequence(self):
        "alias of :external:py:attr:`pysam.AlignedSegment.query_sequence`"
        return self.wrapped.query_sequence

    def get_aligned_pairs(self, *args, **kw):
        "alias of :external:py:meth:`pysam.AlignedSegment.get_aligned_pairs`"
        return self.wrapped.get_aligned_pairs(*args, **kw)

    @functools.cached_property
    def reference_sequence(self) -> str:
        "reference sequence in the wrapped :external:py:class:`pysam.AlignedSegment`; requires the MD tag; see :external:py:meth:`pysam.AlignedSegment.get_reference_sequence`"
        return self.wrapped.get_reference_sequence()

    @functools.cached_property
    def _reference_position_dict(self) -> dict[int, int | None]:
        "a :py:class:`dict` mapping each query position to the corresponding reference position"
        return {
            qry_pos: ref_pos
            for qry_pos, ref_pos in self.get_aligned_pairs()
            if qry_pos is not None
        }

    def get_aligned_reference_position(self, query_position: int) -> int | None:
        """
        Given the `query_position`, returns the corresponding reference position , or `None` if the `query_position` is part of an insertion.

        :param query_position: 0-based query position
        """
        return self._reference_position_dict[query_position]

    def get_aligned_query_sequence(self, reference_position: int) -> str | None:
        """
        Given the `reference_position`, returns the corresponding query sequence if ``self.reference_start <= reference_position < self.reference_end``, or ``None`` otherwise.

        If the `reference_position` is part of a deletion or reference skip, ``"-"`` is returned.

        This function is particularly useful for phasing sequencing reads into halpotype groups.

        :param reference_position: 0-based reference position
        """
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

    def _get_blocks(self) -> list[tuple[int, int]]:
        return self.wrapped.get_blocks()

    def _get_normalized_blocks(self) -> list[tuple[int, int]]:
        blocks: list[tuple[int, int]] = []
        for block_start, block_end in self._get_blocks():
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
    def insertions(self) -> Sequence[Insertion]:
        "insertions identified from CIGAR"
        return self._cigar.insertions

    @property
    def deletions(self) -> Sequence[Deletion]:
        "deletions identified from CIGAR"
        return self._cigar.deletions

    @property
    def soft_clipping(self) -> Sequence[SoftClippedBases]:
        "soft clipping identified from CIGAR"
        return self._cigar.soft_clippings

    @property
    def hard_clipping(self) -> Sequence[HardClippedBases]:
        "hard clipping identified from CIGAR"
        return self._cigar.hard_clippings

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
    def mismatched_bases(self) -> Sequence[_MismatchedBase]:
        """
        mismatched bases identified from CIGAR X (BAM_CDIFF) operation or the MD tag; empty if neither is available; see :external:attr:`pysam.AlignedSegment.cigartuples` and :external:meth:`pysam.AlignedSegment.get_aligned_pairs`

        .. note::
           For `Minimap2 <https://www.python.org/>`_, use ``--eqx`` to include mismatch information in the CIGAR string, ``--MD`` to include mismatch information in the MD tag.
        """
        mismatched_bases: Sequence[_MismatchedBase]
        if self._cigar.mismatched_bases:
            mismatched_bases = self._cigar.mismatched_bases
        elif self.has_tag("MD"):
            mismatched_bases = self._get_md_mismatched_bases()
        else:
            mismatched_bases = []
        return mismatched_bases

    @functools.cached_property
    def modified_bases(self) -> list[ModifiedBase]:
        "modified bases identified from the Ml and Mm tags; see :external:attr:`pysam.AlignedSegment.modified_bases`"
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
        """
        reference skips identified from CIGAR

        .. note::
           Reference skips often appear in RNAseq read alignment to represent introns.
        """
        return self._cigar.reference_skips

    @functools.cached_property
    def _pair_identifier(self) -> tuple[int, int] | None:
        if self.wrapped.mate_is_unmapped:  # Mate is not mapped
            return None
        if (
            self.wrapped.reference_id != self.wrapped.next_reference_id
        ):  # Mate is mapped to another reference sequence
            return None
        self_reference_start: int = self.reference_start
        mate_reference_start: int = self.wrapped.next_reference_start
        min_reference_start: int = min([self_reference_start, mate_reference_start])
        max_reference_start: int = max([self_reference_start, mate_reference_start])
        return (min_reference_start, max_reference_start)


@dataclass
class _LinkedSegment:
    segments: list[AlignedSegment]

    @property
    def reference_start(self) -> int:
        return min(seg.reference_start for seg in self.segments)

    @property
    def reference_end(self) -> int:
        return max(seg.reference_end for seg in self.segments)


class SequenceAlignment:
    """
    Plot sequence alignment from BAM files.
    """

    reference_name: str
    """Name of the reference sequence."""
    segments: list[AlignedSegment]
    """
    A list of aligned segments loaded, in the order that appears in the source file.
    
    .. note::
       For advanced users, it is possible to modify this attribute prior to plotting to apply complex filtering/sorting operations.
    """
    pileup_depths: dict[int, int]
    """A dict mapping reference positions to pileup depths."""
    pileup_bases: dict[int, collections.Counter[str]]
    """A dict mapping reference positions to a counter, which in turn maps each base to its count."""

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
        Load sequence alignment from a local BAM file.

        .. note::
           Python file objects are not supported due to `a known limitation of pysam <https://pysam.readthedocs.io/en/latest/usage.html#using-streams>`_.

        .. note::
           The SAM format is not supported as it does not allow random access.

        :param file_path: Path to the BAM file.
        :param region: region to load alignment records from. See :py:mod:`lakeview.region_notation` for additional details.
        :param kw: Keyword arguments passed to :external:class:`pysam.AlignmentFile`.
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
        """
        Load sequence alignment from a remote BAM file via HTTP, HTTPS, or FTP protocol.

        .. note::
           The SAM format is not supported as it does not allow random access.

        :param url: URL to the remote BAM file
        :param region: region to load alignment records from. See :py:mod:`lakeview.region_notation` for additional details.
        :param index_url: URL to the remote BAM index (.bai) file. The default is ``url + '.bai'``
        :param kw: keyword arguments passed to :external:class:`pysam.AlignmentFile`
        """
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
        segments: Sequence[AlignedSegment], link_identifiers: Sequence[LinkIdentifier]
    ) -> dict[LinkIdentifier, _LinkedSegment]:
        link_seg_dict = collections.defaultdict(list)
        for seg, link in zip(segments, link_identifiers):
            link_seg_dict[link].append(seg)
        link_ls_dict = {}
        for link, seg_list in link_seg_dict.items():
            linked_segment = _LinkedSegment(seg_list)
            link_ls_dict[link] = linked_segment
        return link_ls_dict

    @staticmethod
    def _pack_segments(
        segments: Sequence[AlignedSegment],
        link_identifiers: Sequence[LinkIdentifier],
        *,
        padding: float = 0,
        max_offset: float = float("inf"),
    ) -> np.ndarray:
        assert len(segments) == len(link_identifiers)
        # Link segments
        link_ls_dict = SequenceAlignment._link_segments(segments, link_identifiers)
        # Get offset for each LinkedSegment
        intervals = []
        for link, ls in link_ls_dict.items():
            intervals.append((ls.reference_start - padding, ls.reference_end + padding))
        link_offset_dict = {
            link: offset
            for link, offset in zip(
                link_ls_dict, pack_intervals(intervals, max_offset=max_offset)
            )
        }
        # Retrive offset for each individual segment
        segment_offsets = np.array(
            [link_offset_dict[link] for link in link_identifiers], dtype=np.float32
        )
        return segment_offsets

    @staticmethod
    def _get_segment_offsets(
        segments: Sequence[AlignedSegment],
        link_identifiers: Sequence[LinkIdentifier],
        group_identifiers: Sequence[GroupIdentifier],
        *,
        max_group_offset: float,
        min_spacing: float,
    ) -> Sequence[float]:
        # Check each LinkIdentifier only mpas to one GroupIdentifier.
        link_group_dict: dict[LinkIdentifier, GroupIdentifier] = {}
        for link, group in zip(link_identifiers, group_identifiers):
            recorded_group = link_group_dict.get(link, None)
            if recorded_group is not None and recorded_group != group:
                raise ValueError(
                    f"LinkIdentifier {link!r} maps to multiple GroupIdentifier values. Expecting each LinkIdentifier values maps to only one GroupIdentifier value."
                )

        offsets = np.zeros(len(segments))
        y = 0
        for group in list(sorted(set(group_identifiers))):
            group_indices = [i for i, g in enumerate(group_identifiers) if g == group]
            group_segments = [segments[i] for i in group_indices]
            group_links = [link_identifiers[i] for i in group_indices]
            group_offsets = SequenceAlignment._pack_segments(
                group_segments,
                group_links,
                padding=min_spacing / 2,
                max_offset=max_group_offset,
            )
            group_offsets[group_offsets < 0] = float("-inf")
            offsets[group_indices] = group_offsets + y
            y = max(offsets) + 2
        return list(offsets)

    def _get_filter_keys(
        self, filter_by: Callable[[AlignedSegment], bool] | Iterable[bool] | str | None
    ) -> Sequence[bool]:
        segments = self.segments
        filter_keys: Sequence[bool] = []
        if filter_by is None:
            pass
        elif isinstance(filter_by, str):
            if filter_by == "no_secondary":
                filter_keys = [not seg.is_secondary for seg in segments]
            else:
                raise ValueError()
        elif isinstance(filter_by, Iterable):
            filter_keys = list(filter_by)  # type: ignore
            if len(filter_keys) != len(segments):
                raise ValueError()
        elif callable(filter_by):
            filter_keys = [filter_by(seg) for seg in segments]
        else:
            raise TypeError(f"Invalid type for `filter_by`: f{type(filter_by)!r}.")
        return filter_keys

    def _get_sort_keys(
        self,
        sort_by: Callable[[AlignedSegment], Identifier]
        | Iterable[Identifier]
        | Literal["length", "start"]
        | None,
    ) -> Sequence[Identifier]:
        segments = self.segments
        sort_keys: Sequence[Identifier] = []
        if sort_by is None:
            sort_keys = [0] * len(segments)
        elif isinstance(sort_by, str):
            if sort_by == "start":
                sort_keys = [seg.reference_start for seg in segments]
            elif sort_by == "length":
                sort_keys = [-seg.query_alignment_length for seg in segments]
            else:
                raise ValueError()
        elif isinstance(sort_by, Iterable):
            sort_keys = list(sort_by)
            if len(sort_keys) != len(segments):
                raise ValueError()
        elif callable(sort_by):
            sort_keys = [sort_by(seg) for seg in segments]
        else:
            raise TypeError(f"Invalid type for `sort_by`: f{type(sort_by)!r}.")
        return sort_keys

    def _get_segment_colors(
        self,
        color_by: Callable[[AlignedSegment], Color]
        | Iterable[Color]
        | Literal["proper_pair", "strand"]
        | None,
    ) -> Sequence[Color]:
        segments = self.segments
        colors: Sequence[Color] = []
        if color_by is None:
            colors = ["lightgray"] * len(segments)
        if isinstance(color_by, str):
            if color_by == "proper_pair":
                colors = [
                    "lightgray" if segment.is_proper_pair else "firebrick"
                    for segment in segments
                ]
            elif color_by == "strand":
                colors = list(
                    map(
                        lambda segment: "lightgray"
                        if segment.is_forward
                        else "darkgray",
                        segments,
                    )
                )
            else:
                raise ValueError()
        elif isinstance(color_by, Iterable):
            colors = list(color_by)
            if len(colors) != len(segments):
                raise ValueError()
        elif callable(color_by):
            colors = [color_by(seg) for seg in segments]
        return colors

    def _get_link_identifiers(
        self,
        link_by: Callable[[AlignedSegment], LinkIdentifier]
        | Iterable[LinkIdentifier]
        | Literal["name", "pair"]
        | None,
    ) -> Sequence[LinkIdentifier]:
        segments = self.segments
        link_identifiers: Sequence[LinkIdentifier] = []
        if link_by is None:
            link_identifiers = list(
                range(len(segments))
            )  # Do not link segments together
        elif isinstance(link_by, str):
            if link_by == "pair":
                link_identifiers = []
                link_id: tuple[int, int]
                for i, segment in enumerate(segments):
                    pair_id = segment._pair_identifier
                    if pair_id is None:
                        link_id = (
                            -i - 1,
                            -i - 1,
                        )  # Arbitrary id to make sure non-paired reads are not linked together
                    else:
                        link_id = pair_id
                    link_identifiers.append(link_id)
            elif link_by == "name":
                link_identifiers = [seg.query_name for seg in segments]
            else:
                raise ValueError()
        elif isinstance(link_by, Iterable):
            link_identifiers = list(link_by)
            if len(link_identifiers) != len(segments):
                raise ValueError()
        elif callable(link_by):
            link_identifiers = [link_by(seg) for seg in segments]
        return link_identifiers

    def _get_group_identifiers(
        self,
        group_by: Callable[[AlignedSegment], GroupIdentifier]
        | Iterable[GroupIdentifier]
        | Literal["haplotype", "name", "proper_pair", "strand"]
        | None,
    ) -> Sequence[GroupIdentifier]:
        segments = self.segments
        group_identifiers: Sequence[GroupIdentifier] = []
        if group_by is None:
            group_identifiers = [0] * len(
                segments
            )  # Assign all segments into the same group
        if isinstance(group_by, str):
            if group_by == "haplotype":
                group_identifiers = [
                    seg.get_tag("HP") if seg.has_tag("HP") else float("inf")
                    for seg in segments
                ]
            elif group_by == "name":
                group_identifiers = [seg.query_name for seg in segments]
            elif group_by == "proper_pair":
                group_identifiers = [seg.is_proper_pair for seg in segments]
            elif group_by == "strand":
                group_identifiers = [
                    "forward" if seg.is_forward else "reverse" for seg in segments
                ]
            else:
                raise ValueError()
        elif isinstance(group_by, Iterable):
            group_identifiers = list(group_by)
            if len(group_identifiers) != len(segments):
                raise ValueError()
        if callable(group_by):
            group_identifiers = [group_by(seg) for seg in segments]
        return group_identifiers

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
        ax: Axes,
        *,
        filter_by: Callable[[AlignedSegment], bool]
        | Iterable[bool]
        | str
        | None = None,
        sort_by: Callable[[AlignedSegment], Identifier]
        | Iterable[Identifier]
        | Literal["length", "start"]
        | None = None,
        link_by: Callable[[AlignedSegment], LinkIdentifier]
        | Iterable[LinkIdentifier]
        | Literal["name", "pair"]
        | None = None,  # TODO: link by pair
        group_by: Callable[[AlignedSegment], GroupIdentifier]
        | Iterable[GroupIdentifier]
        | Literal["haplotype", "name", "proper_pair", "strand"]
        | None = None,
        color_by: Callable[[AlignedSegment], Color]
        | Iterable[Color]
        | Literal["proper_pair", "strand"]
        | None = None,
        group_labels: Callable[[GroupIdentifier], str]
        | Mapping[GroupIdentifier, str]
        | None = None,
        height: Optional[float] = None,
        max_rows: int = 200,
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
    ):
        """
        Draw sequence alignment patterns, in a style similar to `IGV alignment track <https://software.broadinstitute.org/software/igv/alignmentdata#alignmenttrack>`_.

        Segments are drawn from top to bottom on the given ``ax``.

        :param ax: Matplotlib :external:class:`matplotlib.axes.Axes` instance to plot on.
        :param filter_by: Specify which aligned segments plotted. The default is to plot all aligned segments.
        :param sort_by: Specify how aligned segments are sorted before plotting. The default is to keep the order of segments in the source file.
        :param link_by: Specify which segments are linked together before plotting. When two or more segments are linked, they will be plotted in the same row, even if they overlap with each other. The default is to link no segments together.
        :param group_by: Specify which segments are groupped together before plotting. Groupped segments are plotted next to each other. The default is not to link any segments. The default is to group no segments together.
        :param color_by: Specify segments colors. The default is to color all segments in light gray.
        :param group_labels: Text labels for each group. The default is to use group identifiers as labels. This parameter is only valid when a custom value for `group_by` is used.
        :param height: Segment height in points. The default is to infer automatically.
        :param min_spacing: The minimum horizontal spacing between two adjacent segments in the same row, in terms of number of bases. The default is to infer automatically.
        :param max_rows: The maximum number of rows to layout segments. Excess segments will not be drawn. If multiple segment groups exist, this parameter limits the maximum number of rows *per group*.
        :param show_backbones: Whether to show segment backbones (i.e. horizontal bars representing aligned positions of each segment).
        :param show_arrowheads: Whether to show triangular arrowheads to denote segment orientations.
        :param show_links: Whether to show horizontal lines connecting linked segments. Only has effect when ``link_by`` is specified.
        :param show_insertions: Whether to show the "I"-shaped markers to denote insertions.
        :param min_insertion_size: Insertions below ``min_insertion_size`` will not be shown.
        :param show_deletions: Whether to show the horizontal dash markers to denote deletions.
        :param min_deletion_size: Deletions below ``min_deletion_size`` will not be shown.
        :param show_mismatches: Whether to show vertical lines denoting mismatched bases. Requires the CIGAR X (BAM_CDIFF) operation or the MD tag. See :py:attr:`lakeview.alignment.AlignedSegment.mismatched_bases` for details.
        :param show_modified_bases: Whether to show vertical lines denoting DNA modification. Requires the Ml and Mm tags. See :py:attr:`lakeview.alignment.AlignedSegment.modified_bases` for details.
        :param show_soft_clippings: Whether to show markers at segment ends denoting soft-clipped bases.
        :param min_soft_clipping_size: Soft-clipping below ``min_soft_clipping_size`` will not be shown.
        :param show_hard_clippings: Whether to show markers at segment ends denoting hard-clipped bases.
        :param min_hard_clipping_size: Hard-clipping below ``min_hard_clipping_size`` will not be shown.
        :param show_reference_skips: Whether to show horizontal lines denoting reference skips. Reference skips are commonly found in intron regions of RNA-seq reads.
        :param show_group_labels: Whether to show a text label for each segment group. If ``None``, group labels will only be shown when ``group_by`` is specified.
        :param show_group_separators: Whether to show horizontal separator lines between adjacent segment groups.

        .. note::
           For a detailed explaination on the usage of ``filter_by``, ``sort_by``, ``link_by``, ``group_by``, and ``color_by``, see :ref:`Customising alignment layout`.

        """
        segments: Sequence[AlignedSegment] = self.segments

        colors = self._get_segment_colors(color_by)
        group_identifiers = self._get_group_identifiers(group_by)
        link_identifiers = self._get_link_identifiers(link_by)
        filter_keys = self._get_filter_keys(filter_by)
        sort_keys = self._get_sort_keys(sort_by)

        # Filter segments
        if filter_by is not None:
            segments = key_filter(segments, filter_keys)
            colors = key_filter(colors, filter_keys)
            link_identifiers = key_filter(link_identifiers, filter_keys)
            group_identifiers = key_filter(group_identifiers, filter_keys)
            if sort_by is not None:
                sort_keys = key_filter(sort_keys, filter_keys)
            if not segments:
                warnings.warn("No segment remains after filtering.")

        # Sort segments
        if sort_by is not None:
            segments = key_sort(segments, sort_keys)
            colors = key_sort(colors, sort_keys)
            link_identifiers = key_sort(link_identifiers, sort_keys)
            group_identifiers = key_sort(group_identifiers, sort_keys)

        # Get default spacing
        if min_spacing is None:
            min_spacing = self._get_default_spacing(segments)

        # Get segment offsets
        max_group_offset = max_rows - 1
        offsets = self._get_segment_offsets(
            segments,
            link_identifiers,
            group_identifiers,
            max_group_offset=max_group_offset,
            min_spacing=min_spacing,
        )

        # Remove segments exceeding `max_rows`
        if not all(y >= 0 for y in offsets):
            filter_keys = [y >= 0 for y in offsets]
            segments = key_filter(segments, filter_keys)
            link_identifiers = key_filter(link_identifiers, filter_keys)
            group_identifiers = key_filter(group_identifiers, filter_keys)
            offsets = key_filter(offsets, filter_keys)
            colors = key_filter(colors, filter_keys)

        # Parse group labels
        unique_group_identifiers = set(group_identifiers)
        group_label_dict: dict[GroupIdentifier, str] = {}
        if group_labels is None:
            if group_by == "strand":
                group_label_dict = {
                    "forward": "Forward strand",
                    "reverse": "Reverse strand",
                }
            elif group_by == "haplotype":
                group_label_dict = {hp: f"Haplotype {hp}" for hp in group_identifiers}
                group_label_dict[float("inf")] = "Haplotype unknown"
            else:
                group_label_dict = {g: str(g) for g in unique_group_identifiers}
        elif callable(group_labels):
            group_label_dict = {g: group_labels(g) for g in unique_group_identifiers}
        elif isinstance(group_labels, Mapping):
            group_label_dict = {g: group_labels[g] for g in unique_group_identifiers}

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
            self._draw_links(ax, segments, offsets, link_identifiers, **links_kw)
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
                ax,
                group_identifiers,
                offsets,
                group_labels=group_label_dict,
                **group_labels_kw,
            )
        if show_group_separators is True:
            self._draw_group_separators(
                ax, group_identifiers, offsets, **group_separators_kw
            )

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
        markeredgewidth=1,
        **kw,
    ) -> None:
        MARKER = Path([(0, 0.5), (0, -0.5)], readonly=True)
        lines: list[tuple[Point, Point]] = []
        line_colors: list[Color] = []
        marker_xs: list[float] = []
        marker_ys: list[float] = []
        marker_colors: list[Color] = []
        for seg, y, color in zip(segments, offsets, colors):
            for block_start, block_end in seg._get_normalized_blocks():
                start_point = (block_start - 0.5, y)
                end_point = (block_end - 0.5, y)
                lines.append((start_point, end_point))
                line_colors.append(color)
                marker_xs.append((block_start + block_end) / 2 - 0.5)
                marker_ys.append(y)
        ax.add_collection(
            LineCollection(
                lines,
                linewidths=height,
                colors=colors[0] if len(set(colors)) == 1 else line_colors,
                zorder=0,
                facecolors="none",
            )
        )
        MARKER_ZORDER = 0.05
        if len(set(marker_colors)) == 1:
            ax.plot(
                marker_xs,
                marker_ys,
                marker=MARKER,
                markersize=height,
                color=marker_colors[0],
                linestyle="",
                markeredgewidth=markeredgewidth,
                zorder=MARKER_ZORDER,
            )
        else:
            ax.scatter(
                marker_xs,
                marker_ys,
                c=marker_colors,
                marker=MARKER,
                s=height**2,
                linewidth=markeredgewidth,
                zorder=MARKER_ZORDER,
                **kw,
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
        group_identifiers,
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
        for g, y in zip(group_identifiers, offsets):
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
        self,
        ax,
        group_identifiers,
        offsets,
        *,
        linewidth=1,
        color="gray",
        linestyle="-",
        **kw,
    ):
        if len(set(group_identifiers)) <= 1:  # Only one group
            return
        max_group_offset_dict = collections.defaultdict(lambda: 0)
        for g, y in zip(group_identifiers, offsets):
            max_group_offset_dict[g] = max(max_group_offset_dict[g], y)
        for y in max_group_offset_dict.values():
            ax.axhline(
                y + 1, linewidth=linewidth, color=color, linestyle=linestyle, **kw
            )

    def _draw_links(
        self,
        ax: Axes,
        segments: Sequence[AlignedSegment],
        offsets: Sequence[float],
        link_identifiers: Sequence[LinkIdentifier],
        *,
        linewidth: float = 1.5,
        color: Color = "lightgray",
        linestyle: str = "-",
        **kw,
    ) -> None:
        link_ls_dict: dict[LinkIdentifier, _LinkedSegment] = self._link_segments(
            segments, link_identifiers
        )
        link_offset_dict: dict[LinkIdentifier, float] = {}
        for link, y in zip(link_identifiers, offsets):
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
        show_mismatches: bool = True,
        min_alt_frequency: float = 0.2,
        min_alt_depth: float = 2,
        window_size: int = 1,
        pileup_kw={},
        mismatch_kw={},
    ):
        """
        Draw pileup depths for alignment, in a style similar to `IGV coverage track <https://software.broadinstitute.org/software/igv/alignmentdata#coverage>`_.

        :param ax: Matplotlib :external:class:`matplotlib.axes.Axes` instance to plot on.
        :param show_mismatches: If ``True``, mismatched bases will be shown as colored bars.
        :param min_alt_frequency: Mismatched bases with frequency less than `min_alt_frequency` will not be shown.
        :param min_alt_depth: Mismatched bases with depth less than `min_alt_frequency` will not be shown.
        :param window_size: Show mean depth values in non-overlapping windows. Windows start from integer multiples of `window_size`.
        :param pileup_kw: Keyword arguments passed to :external:meth:`matplotlib.axes.Axes.fill_between`.
        :param mismatch_kw: Keyword arguments passed to :external:meth:`matplotlib.axes.Axes.bar`.
        """
        self._draw_pileup_fill(ax, window_size=window_size, **pileup_kw)
        if show_mismatches:
            self._draw_pileup_mismatches(
                ax,
                min_alt_frequency=min_alt_frequency,
                min_alt_depth=min_alt_depth,
                **mismatch_kw,
            )

    @staticmethod
    def _get_mean_depths_per_window(
        positions: Iterable[int], depths: Iterable[int], window_size: int = 1
    ) -> tuple[Sequence[float], Sequence[float]]:
        """
        Get per-window mean values of `depths` across given `positions`.
        Windows start from integer multiples of `window_size`.

        >>> depths    = [3, 4, 5, 0, 3, 0 ]
        >>> positions = [0, 1, 2, 8, 9, 10]
        >>> # windows :  -------  -  -----
        >>> SequenceAlignment._get_mean_depths_per_window(positions, depths, window_size=3)
        ([1.0, 7.0, 10.0], [4.0, 0.0, 1.0])
        """

        window_centers: list[float] = []
        mean_depths: list[float] = []

        if window_size <= 0:
            raise ValueError(
                f"Invalid value for `window_size`: {window_size:r}. Expecting a positive integer."
            )
        elif window_size == 1:
            window_centers = list(positions)
            mean_depths = list(depths)
        else:
            window_centers = []
            mean_depths = []
            window_start_position: int = 0
            window_depths: list[int] = []
            for position, depth in zip(positions, depths):
                if (
                    position - window_start_position < window_size
                ):  # Continue current window
                    window_depths.append(depth)
                else:  # Calculate mean depth and start a new window
                    window_centers.append(window_start_position + (window_size - 1) / 2)
                    mean_depth = sum(window_depths) / window_size
                    mean_depths.append(mean_depth)
                    window_depths = [depth]
                    window_start_position = position // window_size * window_size
            # Calculate mean depth for the last window
            window_centers.append(window_start_position + (window_size - 1) / 2)
            mean_depth = sum(window_depths) / window_size
            mean_depths.append(mean_depth)

        return window_centers, mean_depths

    def _draw_pileup_fill(
        self,
        ax: Axes,
        *,
        window_size: int = 1,
        facecolor: Color = "lightgray",
        edgecolor: Color = "none",
        linewidth: float = 1,
        **kw,
    ):
        window_centers, mean_depths = self._get_mean_depths_per_window(
            self.pileup_depths, self.pileup_depths.values(), window_size=window_size
        )
        x: list[float] = []
        y: list[float] = []
        for position, depth in zip(window_centers, mean_depths):
            if not x:  # First window
                x.append(position)
                y.append(depth)
            else:
                # Positions with zero coverage need to be inserted explicitly for correct plotting
                while x[-1] < position - window_size:
                    x.append(x[-1] + window_size)
                    y.append(0)
                x.append(position)
                y.append(depth)

        # Remove redundant values to reduce file size when saved as vector graphic
        excluded_indices: set[int] = set(
            i for i in range(1, len(y) - 1) if y[i - 1] == y[i] == y[i + 1]
        )
        x = [xi for i, xi in enumerate(x) if i not in excluded_indices]
        y = [yi for i, yi in enumerate(y) if i not in excluded_indices]

        ax.fill_between(
            x,
            y1=y,
            y2=0,
            step="mid",
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            **kw,
        )
        ax.set_ylim(bottom=0)

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
        **kw,
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
                **kw,
            )
            bottom += counts


class OpticalMapAlignment:
    pass
