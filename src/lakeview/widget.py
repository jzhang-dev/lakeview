#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Optional, Sequence
from warnings import warn
import matplotlib as mpl
import matplotlib.pyplot as plt
from IPython.display import display
import ipywidgets
from .custom_types import Figure, Axes

# TODO: auto x tick formatter

class GenomeViewer:
    def __init__(
        self,
        nrows: int = 1,
        *,
        height_ratios: Optional[Sequence[float]] = None,
        figsize: tuple[float, float] = (10, 10),
    ) -> None:
        self._initial_xlim: Optional[tuple[float, float]] = None
        with plt.ioff():
            fig, axes = plt.subplots(
                nrows=nrows,
                ncols=1,
                figsize=figsize,
                sharex=True,
                sharey=False,
                squeeze=False,
                gridspec_kw=dict(height_ratios=height_ratios),
                constrained_layout=True,
            )
        self.figure: Figure = fig
        self.axes: list[Axes] = list(axes[:, 0])
        self._init_app()

    def _init_app(self) -> None:
        center_widget = ipywidgets.Output()

        zoom_in_button = ipywidgets.Button(description="Zoom in")
        zoom_out_button = ipywidgets.Button(description="Zoom out")
        zoom_in_button.on_click(lambda button: self.zoom(1 / 1.2))
        zoom_out_button.on_click(lambda button: self.zoom(1.2))

        shift_left_button = ipywidgets.Button(description="<<")
        shift_right_button = ipywidgets.Button(description=">>")
        shift_left_button.on_click(lambda __: self.shift(-0.15))
        shift_right_button.on_click(lambda __: self.shift(0.15))

        reset_button = ipywidgets.Button(description="Reset")
        reset_button.on_click(lambda _: self.reset_xlim())

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
        #self.update_display()
        self.app = ipywidgets.AppLayout(
            center=center_widget, footer=footer_widget, pane_heights=[0, 1, "100px"]
        )

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
        self.update_display()
        display(self.app)

    def zoom(self, scale, /):
        for ax in self.figure.axes:
            old_xlim = ax.get_xlim()
            center = sum(old_xlim) / 2
            radius = (old_xlim[1] - old_xlim[0]) / 2
            new_radius = radius * scale
            new_xlim = (center - new_radius, center + new_radius)
            ax.set_xlim(new_xlim)
        self.update_display()

    def shift(self, distance, /):
        for ax in self.figure.axes:
            old_xlim = ax.get_xlim()
            radius = (old_xlim[1] - old_xlim[0]) / 2
            new_xlim = (
                old_xlim[0] + distance * radius,
                old_xlim[1] + distance * radius,
            )
            ax.set_xlim(new_xlim)
        self.update_display()

    def goto(self, start: float, end: float) -> None:
        self.set_xlim(start, end)
        self.update_display()

    def _goto_region(self, region: str) -> None:
        start, end = region.split("-")
        if start.isnumeric() and end.isnumeric():
            self.goto(float(start), float(end))
        else:
            raise ValueError(f"Invalid region {region!r}")

    def reset_xlim(self) -> None:
        initial_xlim = self._initial_xlim
        if initial_xlim is not None:
            self.set_xlim(initial_xlim)
        self.update_display()

    def update_display(self):
        start, end = self.get_xlim()
        self._region_text.value = f"{int(start):,} - {int(end):,}"
        self._center_widget.clear_output(wait=True)
        with self._center_widget:
            display(self.figure)