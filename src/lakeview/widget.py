#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Optional, Sequence, Literal, cast
from math import log10, ceil
from warnings import warn
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from IPython.display import display
import ipywidgets
from ._type_alias import Figure, Axes
from .plot import BasePairFormatter


class GenomeViewer:
    """
    An interactive widget for viewing genomic data in Jupyter Notebooks.
    """

    def __init__(
        self,
        tracks: int = 1,
        *,
        height_ratios: Optional[Sequence[float]] = None,
        figsize: tuple[float, float] = (10, 10),
    ) -> None:
        with plt.ioff():
            fig, axes = plt.subplots(
                nrows=tracks,
                ncols=1,
                figsize=figsize,
                sharex=True,
                sharey=False,
                squeeze=False,
                gridspec_kw=dict(height_ratios=height_ratios),
                constrained_layout=True,
            )
        self._figure = fig
        self._axes = axes[:, 0]
        self._app: ipywidgets.AppLayout | None = None
        self._initial_xlim: Optional[tuple[float, float]] = None

    @property
    def figure(self) -> Figure:
        return self._figure

    @property
    def axes(self) -> Sequence[Axes]:
        return self._axes

    def _init_app(self) -> None:
        center_widget = ipywidgets.Output()

        zoom_in_button = ipywidgets.Button(description="Zoom in")
        zoom_out_button = ipywidgets.Button(description="Zoom out")
        zoom_in_button.on_click(lambda button: self._zoom(1 / 1.2))
        zoom_out_button.on_click(lambda button: self._zoom(1.2))

        shift_left_button = ipywidgets.Button(description="<<")
        shift_right_button = ipywidgets.Button(description=">>")
        shift_left_button.on_click(lambda __: self._shift(-0.15))
        shift_right_button.on_click(lambda __: self._shift(0.15))

        reset_button = ipywidgets.Button(description="Reset")
        reset_button.on_click(lambda _: self._reset_xlim())

        region_text = ipywidgets.Text(
            value="",
            placeholder="",
            description="",
            disabled=False,
        )
        self._region_text = region_text
        go_button = ipywidgets.Button(description="Go")
        go_button.on_click(lambda __: self._goto_region(self._region_text.value))

        footer_widget = ipywidgets.VBox(
            [
                ipywidgets.HBox([region_text, go_button, reset_button]),
                ipywidgets.HBox(
                    [
                        shift_left_button,
                        shift_right_button,
                        zoom_in_button,
                        zoom_out_button,
                    ]
                ),
            ]
        )
        self._center_widget = center_widget
        self._footer_widget = footer_widget
        # self.update_display()
        self._app = ipywidgets.AppLayout(
            center=center_widget, footer=footer_widget, pane_heights=[0, 1, "100px"]
        )

    @property
    def app(self) -> ipywidgets.AppLayout:
        if self._app is None:
            self._init_app()
        return self._app

    @property
    def xaxis(self) -> mpl.axis.XAxis:
        return self.axes[-1].xaxis

    def set_xlim(self, *args, **kw) -> tuple[float, float]:
        return self.axes[0].set_xlim(*args, **kw)

    def get_xlim(self) -> tuple[float, float]:
        return self.axes[0].get_xlim()

    def set_xlabel(self, *args, **kw) -> mpl.text.Text:
        self.axes[-1].set_xlabel(*args, **kw)

    def set_title(self, *args, **kw) -> mpl.text.Text:
        self.axes[0].set_title(*args, **kw)

    def savefig(self, *args, **kw) -> None:
        self.figure.savefig(*args, **kw)

    def _ipython_display_(self) -> None:
        self._initial_xlim = self.get_xlim()
        self._update_display()
        display(self.app)

    def _zoom(self, scale, /):
        for ax in self.figure.axes:
            old_xlim = ax.get_xlim()
            center = sum(old_xlim) / 2
            radius = (old_xlim[1] - old_xlim[0]) / 2
            new_radius = radius * scale
            new_xlim = (center - new_radius, center + new_radius)
            ax.set_xlim(new_xlim)
        self._update_display()

    def _shift(self, distance, /):
        for ax in self.figure.axes:
            old_xlim = ax.get_xlim()
            radius = (old_xlim[1] - old_xlim[0]) / 2
            new_xlim = (
                old_xlim[0] + distance * radius,
                old_xlim[1] + distance * radius,
            )
            ax.set_xlim(new_xlim)
        self._update_display()

    def _goto(self, start: float, end: float) -> None:
        self.set_xlim(start, end)
        self._update_display()

    def _goto_region(self, region: str) -> None:
        start, end = region.split("-")
        if start.isnumeric() and end.isnumeric():
            self._goto(float(start), float(end))
        else:
            raise ValueError(f"Invalid region {region!r}")

    def _reset_xlim(self) -> None:
        initial_xlim = self._initial_xlim
        if initial_xlim is not None:
            self.set_xlim(initial_xlim)
        self._update_display()

    def _update_xaxis_ticklabels(self) -> None:
        start, end = self.get_xlim()
        span = end - start
        unit_divisor: int
        for unit in ("Tb", "Gb", "Mb", "kb", "bp"):
            unit_divisor = BasePairFormatter._get_unit_divisor(
                cast(Literal["bp", "kb", "Mb", "Gb", "Tb"], unit)
            )  # mypy does not infer a string as a literal unless assigned directly
            if unit_divisor <= start:
                break
        decimals: int = max(0, 2 - ceil(log10(span / unit_divisor)))
        formatter = BasePairFormatter(
            cast(Literal["bp", "kb", "Mb", "Gb", "Tb"], unit),
            decimals,
            show_suffix=True,
        )
        self.xaxis.set_major_formatter(formatter)

    def _update_region_text(self) -> None:
        start, end = self.get_xlim()
        self._region_text.value = f"{int(start):,} - {int(end):,}"

    def _update_display(self):
        self._update_region_text()
        self._update_xaxis_ticklabels()
        self._center_widget.clear_output(wait=True)
        with self._center_widget:
            display(self.figure)
