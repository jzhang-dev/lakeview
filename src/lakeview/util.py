#!/usr/bin/env python
# coding: utf-8

import matplotlib as mpl

def base_formatter(unit="mb", fmt="{:.2f}"):
    n = dict(bp=1, kb=int(1e3), mb=int(1e6), gb=int(1e9))[unit.lower()]

    @mpl.ticker.FuncFormatter
    def unit_formatter(x, pos):
        return fmt.format(x / n)

    return unit_formatter


def complement_base(base):
    if base == "A":
        return "T"
    if base == "T":
        return "A"
    if base == "G":
        return "C"
    if base == "C":
        return "G"
    if base == "N":
        return "N"
    else:
        raise ValueError