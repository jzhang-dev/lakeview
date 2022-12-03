#!/usr/bin/env python
# coding: utf-8

from warnings import warn
from dataclasses import dataclass, field
from numbers import Real
from typing import List, Dict, Optional, Tuple, Callable
import matplotlib as mpl
import matplotlib.pyplot as plt
from IPython.display import display
import ipywidgets

# TODO: give up interactivity?
# TODO: more buttons; auto x tick labels; prevent fig resizing
# TODO: fix first display problem; smart scrolling / resizing; show certain artists only when zooming in;
# TODO: Fix interactivity; JS error
# TODO: soft clips should display as lines when not showing arrows


class GenomeViewer:
    def __init__(self, nrows=1, *, height_ratios=None, figsize=(10, 10)):
        self._background = None
        self._vlines = None
        self._initial_xlim = None
        self._interacted = False

        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=1,
            figsize=figsize,
            sharex=True,
            sharey=False,
            squeeze=False,
            gridspec_kw=dict(height_ratios=height_ratios),
            constrained_layout=True
        )
        self.figure = fig
        #self.blit_manager = _BlitManager(fig.canvas)
        self.axes = list(axes[:, 0])
        self._init_app()

    def _init_app(self):
        backend = mpl.get_backend()
        if backend == "module://ipympl.backend_nbagg":
            self.use_ipympl=True
            canvas = self.figure.canvas
            canvas.toolbar_visible = True
            canvas.toolbar_position = "right"
            canvas.header_visible = False
            canvas.footer_visible = False
            canvas.mpl_connect("button_press_event", lambda event: self.on_click(event))
            center_widget = canvas
        else:
            # warn(
            #     f"Backend ipympl not activated. Current backend: {backend}. To enable additional interactive functionality, please activate ipympl using `%matplotlib widget`."
            # )
            self.use_ipympl=False
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
        reset_button.on_click(lambda _: self.reset())

        region_text = ipywidgets.Text(
            value="", placeholder="", description="", disabled=False,
        )
        self.region_text = region_text
        go_button = ipywidgets.Button(description="Go")
        go_button.on_click(lambda __: self.goto(self.region_text.value))

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
        self.center_widget = center_widget
        self.footer_widget = footer_widget
        self.update(interactive=False)
        self.app = ipywidgets.AppLayout(
            center=center_widget, footer=footer_widget, pane_heights=[0, 1, "100px"]
        )

    def set_xlim(self, *args, **kw):
        self.axes[0].set_xlim(*args, **kw)

    def get_xlim(self):
        return self.axes[0].get_xlim()
    
    def set_xlabel(self, *args, **kw):
        self.axes[-1].set_xlabel(*args, **kw)

    def set_title(self, *args, **kw):
        self.axes[0].set_title(*args, **kw)

    def savefig(self, *args, **kw):
        if self._interacted:
            warn("The figure has been modified interactively before saving. This will make it harder for the output to be reproduced in the future.")
        return self.figure.savefig(*args, **kw)

    def _ipython_display_(self):
        self._initial_xlim = self.get_xlim()
        self.update(interactive=False)
        return display(self.app)


    def zoom(self, scale, /):
        for ax in self.figure.axes:
            old_xlim = ax.get_xlim()
            center = sum(old_xlim) / 2
            radius = (old_xlim[1] - old_xlim[0]) / 2
            new_radius = radius * scale
            new_xlim = (center - new_radius, center + new_radius)
            ax.set_xlim(new_xlim)
        self.update()

    def shift(self, distance, /):
        for ax in self.figure.axes:
            old_xlim = ax.get_xlim()
            radius = (old_xlim[1] - old_xlim[0]) / 2
            new_xlim = (
                old_xlim[0] + distance * radius,
                old_xlim[1] + distance * radius,
            )
            ax.set_xlim(new_xlim)
        self.update()

    def goto(self, region):
        start, end = region.split("-")
        start = float("".join(x for x in start if x not in (" ", "\t", ",")))
        end = float("".join(x for x in end if x not in (" ", "\t", ",")))
        for ax in self.figure.axes:
            ax.set_xlim(start, end)
        self.update()


    def show_vertical_line(self, x=None):
        vlines = self._vlines
        if vlines is None and x is not None:
            vlines = [ax.axvline(x, color="k", lw=1, zorder=10, animated=True) for ax in self.axes]
            self._vlines = vlines
        if x is not None:
            for l in vlines:
                l.set_xdata([x, x])
        
        bg = self._background
        if bg:
            self.figure.canvas.restore_region(bg)
            if vlines:
                for ax, l in zip(self.axes, vlines):
                    ax.draw_artist(l)
            self.figure.canvas.blit(self.figure.bbox)
            self.figure.canvas.flush_events()
        else:
            self.update()
        
        
    def hide_vertical_line(self):
        bg = self._background
        if bg:
            self.figure.canvas.restore_region(bg)
            self.figure.canvas.blit(self.figure.bbox)
            self.figure.canvas.flush_events()
        else:
            self.update()

    def on_click(self, event):
        # TODO: use blitz to improve performance
        if event.inaxes:
            self.show_vertical_line(event.xdata)
        else:
            self.hide_vertical_line()

    def reset(self):
        initial_xlim = self._initial_xlim
        if initial_xlim:
            self.set_xlim(initial_xlim)
        vlines = self._vlines
        if vlines:
            for line in vlines:
                line.set_visible(False)
        self.update()
        


    def update(self, *, interactive=True):
        if interactive:
            self._interacted = True
        start, end = self.get_xlim()
        self.region_text.value = f"{int(start):,} - {int(end):,}"
        if self.use_ipympl:
            self.figure.canvas.draw()
        else:
            self.center_widget.clear_output(wait=True)
            with self.center_widget:
                display(self.figure)
        self._background = self.figure.canvas.copy_from_bbox(self.figure.bbox)
        #self.show_vertical_line()

