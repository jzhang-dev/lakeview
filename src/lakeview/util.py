#!/usr/bin/env python
# coding: utf-8

import matplotlib as mpl


# Optional: remove this module

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
    return kmer <= reverse_complement(kmer)


def canonicalize(kmer):
    canonical_kmer = min(kmer, reverse_complement(kmer))
    return canonical_kmer


def iterate_kmers(sequence, k, strand="single", drop_N=True):
    """
    >>> list(iterate_kmers("AATGANGGG", 3, 'canonical'))
    ['AAT', 'ATG', 'TCA', 'CCC']
    """
    for __, kmer in enumerate_kmers(sequence, k, strand, drop_N):
        yield kmer


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