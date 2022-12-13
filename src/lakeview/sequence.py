#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
import collections
import math
import itertools
from dataclasses import dataclass
from typing import Optional
import numpy as np
from Bio import SeqIO

from . import util


Dot = collections.namedtuple("Dot", ["x", "y"])


@dataclass(repr=False)
class DotPlot:
    dots: list[Dot]
    x_sequence_name: Optional[str] = None
    y_sequence_name: Optional[str] = None
    x_sequence_size: Optional[str] = None
    y_sequence_size: Optional[str] = None

    @classmethod
    def from_sequences(
        cls,
        x_sequence,
        y_sequence,
        k=50,
        *,
        sample_fraction=1,
        x_sequence_name=None,
        y_sequence_name=None,
    ):
        if sample_fraction > 1:
            raise ValueError("sample_fraction > 1")
        elif sample_fraction == 1:
            kmer_sample_fraction = 1
        else:
            kmer_sample_fraction = math.sqrt(sample_fraction)
        kmer_indices_x = cls.get_canonical_kmer_indices(
            x_sequence, k, sample_fraction=kmer_sample_fraction
        )
        kmer_indices_y = cls.get_canonical_kmer_indices(
            y_sequence, k, sample_fraction=kmer_sample_fraction
        )
        shared_kmers = set(kmer_indices_x) & set(kmer_indices_y)
        dots = []
        for kmer in shared_kmers:
            x_indices = kmer_indices_x[kmer]
            y_indices = kmer_indices_y[kmer]
            for x, y in itertools.product(x_indices, y_indices):
                dots.append(Dot(x, y))

        instance = cls(
            dots,
            x_sequence_name=x_sequence_name,
            y_sequence_name=y_sequence_name,
            x_sequence_size=len(x_sequence),
            y_sequence_size=len(y_sequence),
        )
        return instance

    @classmethod
    def from_files(
        cls,
        x_file_path=None,
        y_file_path=None,
        k=50,
        *,
        x_file_object=None,
        y_file_object=None,
        x_format="fasta",
        y_format="fasta",
        sample_fraction=1,
        x_sequence_name=None,
        y_sequence_name=None,
    ):
        x_name, x_sequence = cls.get_sequence(
            x_file_path, x_file_object, format=x_format, sequence_name=x_sequence_name
        )
        y_name, y_sequence = cls.get_sequence(
            y_file_path, y_file_object, format=y_format, sequence_name=y_sequence_name
        )
        return cls.from_sequences(
            x_sequence,
            y_sequence,
            k=k,
            sample_fraction=sample_fraction,
            x_sequence_name=x_name,
            y_sequence_name=y_name,
        )

    # TODO: Fix reading from the same handle twice
    @staticmethod
    def get_sequence(
        file_path=None, file_object=None, format="fasta", sequence_name=None
    ):
        sequence = None
        for i, record in enumerate(
            SeqIO.parse(file_object or file_path, format=format)
        ):
            if sequence_name is None:
                if i == 0:
                    name = record.id
                    sequence = str(record.seq)
                else:
                    raise ValueError(
                        "Found multiple sequences. Please specify the `sequence_name`."
                    )
            else:
                if record.id == sequence_name:
                    name = sequence_name
                    sequence = str(record.seq)
                    break
        if sequence is None:
            raise ValueError(f"Sequence {sequence_name!r} not found.")
        return name, sequence

    @staticmethod
    def get_canonical_kmer_indices(sequence, k, *, sample_fraction=1, seed=0):
        rng = np.random.default_rng(seed=seed)
        kmer_indices = collections.defaultdict(list)
        for index, kmer in util.enumerate_kmers(sequence, k, strand="single"):
            if sample_fraction == 1 or rng.random() <= sample_fraction:
                canonical_kmer = util.canonicalize(kmer)
                kmer_indices[canonical_kmer].append(index)
        return kmer_indices

    # TODO: max_dots
    def draw_dots(
        self, ax, *, x_offset=0, y_offset=0, s=0.5, max_dots=1e6, edgecolor="none"
    ):
        xs = np.array([dot.x for dot in self.dots]) + x_offset
        ys = np.array([dot.y for dot in self.dots]) + y_offset
        ax.scatter(xs, ys, s=s, edgecolor=edgecolor)
        ax.set_xlim(x_offset, x_offset + self.x_sequence_size - 1)
        ax.set_ylim(y_offset, y_offset + self.y_sequence_size - 1)
        if self.x_sequence_name is not None:
            ax.set_xlabel(self.x_sequence_name)
        if self.y_sequence_name is not None:
            ax.set_ylabel(self.y_sequence_name)

    def draw_heatmap(
        self,
        ax,
        *,
        x_offset=0,
        y_offset=0,
        bin_size=1000,
        density=True,
        cmap="cividis",
        cmin=1,
        **kw,
    ):
        xs = np.array([dot.x for dot in self.dots]) + x_offset
        ys = np.array([dot.y for dot in self.dots]) + y_offset
        bins = [
            np.arange(xs.min(), xs.max() + 1, bin_size),
            np.arange(ys.min(), ys.max() + 1, bin_size),
        ]
        ax.hist2d(xs, ys, bins=bins, density=density, cmap=cmap, cmin=cmin, **kw)
        ax.set_xlim(x_offset, x_offset + self.x_sequence_size - 1)
        ax.set_ylim(y_offset, y_offset + self.y_sequence_size - 1)
        if self.x_sequence_name is not None:
            ax.set_xlabel(self.x_sequence_name)
        if self.y_sequence_name is not None:
            ax.set_ylabel(self.y_sequence_name)

    def draw_self_dots(self, axis="x"):
        pass

    def draw_self_heatmap(self, axis="x"):
        pass
