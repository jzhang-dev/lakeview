#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
import collections
import math
import itertools
from dataclasses import dataclass, field
from typing import Literal, Sequence, Collection, Mapping
import numpy as np
from Bio import SeqIO


# Type aliases
Kmer = str


def init_reverse_complement():
    TRANSLATION_TABLE = str.maketrans("ACTGactg", "TGACtgac")

    def reverse_complement(sequence):
        """
        >>> reverse_complement("AATC")
        'GATT'
        >>> reverse_complement("CCANT")
        'ANTGG'
        """
        sequence = str(sequence)
        return sequence.translate(TRANSLATION_TABLE)[::-1]

    return reverse_complement


reverse_complement = init_reverse_complement()


def is_canonical(kmer):
    """
    >>> is_canonical("ATC")
    True
    >>> is_canonical("GGC")
    False
    """
    # This function uses the alphabetical/lexicographical ("ACGT") order to determine canonical kmers. 
    # JellyFish also uses the "ACGT" order, while Meryl uses the "ACTG" order by default. See https://github.com/marbl/meryl/blob/69c839cf43b774916b9a004f8143a00d4f9aac86/src/meryl2/merylOp-nextMer.C#L30
    return kmer <= reverse_complement(kmer)


def canonicalize(kmer):
    canonical_kmer = min(kmer, reverse_complement(kmer))
    return canonical_kmer


def enumerate_kmers(sequence, k, strand="single", drop_N=True):
    """
    >>> list(enumerate_kmers("AATGANGGG", 3, 'canonical'))
    [(0, 'AAT'), (1, 'ATG'), (2, 'TCA'), (6, 'CCC')]
    """
    if len(sequence) < k:
        raise ValueError
    if strand not in ("single", "canonical"):
        raise ValueError

    sequence = str(sequence).upper()
    l = len(sequence)
    for i in range(l):
        if i + k > l:
            return
        kmer = sequence[i : i + k]
        if drop_N and "N" in kmer:
            continue
        if strand == "single":
            yield i, kmer
        elif strand == "canonical":
            yield i, canonicalize(kmer)


def load_multiple_sequences(
    file_object,
    sequence_names: Collection[str],
    format_: Literal["fasta", "fastq"] = "fasta",
) -> Mapping[str, str]:
    sequence_dict: dict[str, str] = {}
    missing_names: list[str] = list(sequence_names)
    for record in SeqIO.parse(file_object, format=format_):
        if record.id in sequence_names:
            sequence_dict[record.id] = str(record.seq)
            missing_names = [name for name in missing_names if name not in sequence_dict]
            if not missing_names:
                break
    if missing_names:
        raise ValueError(f"Missing sequences: {missing_names!r}")
    return sequence_dict


def load_sequence(
    file_object, sequence_name: str, format_: Literal["fasta", "fastq"] = "fasta"
) -> str:
    sequence_dict = load_multiple_sequences(
        file_object, sequence_names=[sequence_name], format_=format_
    )
    return sequence_dict[sequence_name]


Dot = collections.namedtuple("Dot", ["x", "y"])


@dataclass
class DotPlot:
    dots: Sequence[Dot] = field(repr=False)
    x_sequence_size: int = field()
    y_sequence_size: int = field()

    @classmethod
    def from_sequences(
        cls,
        x_sequence,
        y_sequence,
        k=50,
        *,
        sample_fraction=1,
        x_start=0,  # TODO
        x_end=float("inf"),
        y_start=0,
        y_end=float("inf"),
    ):
        if not (0 < sample_fraction <= 1):
            raise ValueError(
                f"Invalid value for `sample_fraction`: {sample_fraction!r}. Expecting `0 < sample_fraction <= 1`."
            )
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
            dots, x_sequence_size=len(x_sequence), y_sequence_size=len(y_sequence)
        )
        return instance

    @staticmethod
    def get_canonical_kmer_indices(
        sequence: str,
        k: int,
        *,
        sample_fraction: float = 1,
        seed: int = 0,
        start: int = 0,  # TODO
        end: int = int(1e12),
    ) -> Mapping[Kmer, Sequence[int]]:
        rng = np.random.default_rng(seed=seed)
        kmer_indices = collections.defaultdict(list)
        for index, kmer in enumerate_kmers(sequence, k, strand="single"):
            if sample_fraction == 1 or rng.random() <= sample_fraction:
                canonical_kmer = canonicalize(kmer)
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
