#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Optional, Sequence, Literal, cast
from dataclasses import dataclass
from math import log10, ceil
from warnings import warn
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from IPython.display import display
import ipywidgets
from ._type_alias import Figure, Axes
from .plot import BasePairFormatter


class _GenomeViewerWidget:
    def __init__(self, figure: Figure, axes: Sequence[Axes]):
        self.figure: Figure = figure
        self.axes: Axes = axes
        self._initial_xlim: tuple[float, float] = self.get_xlim()

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

        self._app = ipywidgets.AppLayout(
            center=center_widget, footer=footer_widget, pane_heights=[0, 1, "100px"]
        )

    @property
    def app(self) -> ipywidgets.AppLayout:
        return self._app

    def _ipython_display_(self) -> None:
        self._initial_xlim = self.get_xlim()
        self._update_display()
        display(self.app)

    def get_xlim(self) -> tuple[float, float]:
        return self.axes[0].get_xlim()

    def set_xlim(self, *args, **kw) -> tuple[float, float]:
        return self.axes[0].set_xlim(*args, **kw)

    def _zoom(self, scale: float, /) -> None:
        for ax in self.figure.axes:
            old_xlim = ax.get_xlim()
            center = sum(old_xlim) / 2
            radius = (old_xlim[1] - old_xlim[0]) / 2
            new_radius = radius * scale
            new_xlim = (center - new_radius, center + new_radius)
            ax.set_xlim(new_xlim)
        self._update_display()

    def _shift(self, distance: float, /) -> None:
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
        self.set_xlim(self._initial_xlim)
        self._update_display()

    def _update_xaxis_ticklabels(self) -> None:
        formatter = BasePairFormatter.from_limits(self.get_xlim())
        self.axes[-1].xaxis.set_major_formatter(formatter)

    def _update_region_text(self) -> None:
        start, end = self.get_xlim()
        self._region_text.value = f"{int(start):,} - {int(end):,}"

    def _update_display(self):
        self._update_region_text()
        self._update_xaxis_ticklabels()
        self._center_widget.clear_output(wait=True)
        with self._center_widget:
            display(self.figure)


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
        xlim: tuple[float, float] | None = None,
        use_tick_formatter: bool = True,
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
        self._widget: _GenomeViewerWidget | None = None
        self._use_tick_formatter:bool = use_tick_formatter
        if xlim is not None:
            self.set_xlim(xlim)

    @property
    def figure(self) -> Figure:
        return self._figure

    @property
    def axes(self) -> Sequence[Axes]:
        return self._axes

    @property
    def widget(self) -> _GenomeViewerWidget:
        if self._widget is None:
            widget = _GenomeViewerWidget(self.figure, self.axes)
            self._widget = widget
        return self._widget

    @property
    def xaxis(self) -> mpl.axis.XAxis:
        return self.axes[-1].xaxis

    def set_xlim(self, *args, **kw) -> tuple[float, float]:
        start, end = self.axes[0].set_xlim(*args, **kw)
        if self._use_tick_formatter:
            formatter = BasePairFormatter.from_limits(self.get_xlim())
            self.axes[-1].xaxis.set_major_formatter(formatter)
        return start, end

    def get_xlim(self) -> tuple[float, float]:
        return self.axes[0].get_xlim()

    def use_base_pair_formatter(self) -> None:
        formatter = BasePairFormatter.from_limits(self.get_xlim())
        self.axes[-1].xaxis.set_major_formatter(formatter)

    def set_xlabel(self, *args, **kw) -> mpl.text.Text:
        self.axes[-1].set_xlabel(*args, **kw)

    def set_title(self, *args, **kw) -> mpl.text.Text:
        self.axes[0].set_title(*args, **kw)

    def savefig(self, *args, dpi=300, bbox_inches="tight", **kw):
        return self.figure.savefig(*args, dpi=300, bbox_inches=bbox_inches, **kw)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(figure={self.figure!r})"

    def clear(self) -> None:
        """
        Clear all contents in each Axes.
        """
        for ax in self.axes:
            ax.cla()
