#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass, field
from numbers import Real
from typing import List, Dict, Optional, Tuple, Callable
import matplotlib as mpl
import matplotlib.pyplot as plt
from IPython.display import display
import ipywidgets


# TODO: more buttons; auto x tick labels; prevent fig resizing
# TODO: fix first display problem; smart scrolling / resizing; show certain artists only when zooming in;

# TODO: enable region as parameter
# TODO: rewrite GV to emulate mpl api
@dataclass(repr=False)
class GenomeViewer:
    tracks: List[Callable] = field(default_factory=list)
    track_params: List[Dict] = field(default_factory=list)
    width: Real = 10
    height: Real = 12
    xlim: Optional[Tuple[Real, Real]] = None
    use_ipympl: bool = True

    def add_track(self, track, relative_height=1, ylim=None, ylabel=None, **kw):
        self.tracks.append(track)
        self.track_params.append(
            {
                "_relative_height": relative_height,
                "_ylim": ylim,
                "_ylabel": ylabel,
                **kw,
            }
        )
        self.figure = None
        self.app = None
        # TODO: use traitlets to monitor track changes

    def __post_init__(self):
        if len(self.tracks) != len(self.track_params):
            raise ValueError
        self.figure = None
        self.app = None

    def init_figure(self):
        height_ratios = [params["_relative_height"] for params in self.track_params]
        fig, axes = plt.subplots(
            len(self.tracks),
            1,
            squeeze=False,
            figsize=(self.width, self.height),
            gridspec_kw=dict(height_ratios=height_ratios),
            sharex=True,
            constrained_layout=True,
        )
        for ax, painter, kw in zip(axes[:, 0], self.tracks, self.track_params):
            ylim = kw.get("_ylim")
            ylabel = kw.get("_ylabel")
            kw = {k: v for k, v in kw.items() if not k.startswith("_")}
            painter(ax, **kw)
            if ylim is not None:
                ax.set_ylim(ylim)
            if ylim is not None:
                ax.set_ylabel(ylabel)
        self.figure = fig

    def init_app(self):
        if self.figure is None:
            self.init_figure()
        if self.use_ipympl:
            backend = mpl.get_backend()
            if backend == "module://ipympl.backend_nbagg":
                canvas = self.figure.canvas
                canvas.toolbar_visible = True
                canvas.toolbar_position = "right"
                canvas.header_visible = False
                canvas.footer_visible = False
                self.vlines = []
                canvas.mpl_connect(
                    "button_press_event", lambda event: self.on_click(event)
                )
                center_widget = canvas
            else:
                raise RuntimeError(
                    f"Please activate ipympl backend or set `use_ipympl` to `False`. Current backend: {backend}"
                )
        else:
            center_widget = ipywidgets.Output()

        zoom_in_button = ipywidgets.Button(description="Zoom in")
        zoom_out_button = ipywidgets.Button(description="Zoom out")
        zoom_in_button.on_click(lambda button: self.zoom(1 / 1.2))
        zoom_out_button.on_click(lambda button: self.zoom(1.2))

        shift_left_button = ipywidgets.Button(description="<<")
        shift_right_button = ipywidgets.Button(description=">>")
        shift_left_button.on_click(lambda __: self.shift(-0.15))
        shift_right_button.on_click(lambda __: self.shift(0.15))

        region_text = ipywidgets.Text(
            value="", placeholder="", description="", disabled=False,
        )
        self.region_text = region_text
        go_button = ipywidgets.Button(description="Go")
        go_button.on_click(lambda __: self.goto(self.region_text.value))

        footer_widget = ipywidgets.VBox(
            [
                ipywidgets.HBox([region_text, go_button]),
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

        self.center_widget = center_widget
        self.footer_widget = footer_widget
        self.update_center_widget()
        self.app = ipywidgets.AppLayout(
            center=center_widget, footer=footer_widget, pane_heights=[0, 1, "100px"]
        )
        self.reset_xlim()

    def _ipython_display_(self):
        if not self.tracks:
            return display(repr(self))
        if self.app is None:
            self.init_app()
        return display(self.app)

    def reset_xlim(self):
        if self.xlim is not None and self.figure is not None:
            self.figure.axes[0].set_xlim(self.xlim)

    def zoom(self, scale, /):
        for ax in self.figure.axes:
            old_xlim = ax.get_xlim()
            center = sum(old_xlim) / 2
            radius = (old_xlim[1] - old_xlim[0]) / 2
            new_radius = radius * scale
            new_xlim = (center - new_radius, center + new_radius)
            ax.set_xlim(new_xlim)
        self.update_center_widget()

    def shift(self, distance, /):
        for ax in self.figure.axes:
            old_xlim = ax.get_xlim()
            radius = (old_xlim[1] - old_xlim[0]) / 2
            new_xlim = (
                old_xlim[0] + distance * radius,
                old_xlim[1] + distance * radius,
            )
            ax.set_xlim(new_xlim)
        self.update_center_widget()

    def goto(self, region):
        start, end = region.split("-")
        start = float("".join(x for x in start if x not in (" ", "\t", ",")))
        end = float("".join(x for x in end if x not in (" ", "\t", ",")))
        for ax in self.figure.axes:
            ax.set_xlim(start, end)
        self.update_center_widget()

    def on_click(self, event):
        # TODO: use blitz to improve performance
        for line in self.vlines:
            line.set_visible(False)
        if event.inaxes:
            for ax in self.figure.axes:
                line = ax.axvline(event.xdata, color="k", lw=1, zorder=10)
                self.vlines.append(line)
        self.update_center_widget()

    def update_center_widget(self):
        start, end = self.figure.axes[0].get_xlim()
        self.region_text.value = f"{int(start):,} - {int(end):,}"
        if self.use_ipympl:
            self.figure.canvas.draw()
            self.figure.canvas.flush_events()
        else:
            self.center_widget.clear_output(wait=True)
            with self.center_widget:
                display(self.figure)

