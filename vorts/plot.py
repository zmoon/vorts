"""
Plotting routines
"""

from cycler import cycler
import matplotlib.pyplot as plt
# import numpy as np


# Tableau's newer version of tab10
_NEW_TAB10 = [
    "#4e79a7",
    "#f28e2b",
    "#e15759",
    "#76b7b2",
    "#59a14f",
    "#edc948",
    "#b07aa1",
    "#ff9da7",
    "#9c755f",
    "#bab0ac",
]

def plot_vorton_trajectories(ds, **kwargs):
    """Plot lines: one for each vorton's trajectory.

    **kwargs are passed on to `plt.subplots()`
    """
    # select vortons
    iv = ds.G != 0
    ds = ds.sel(v=iv)

    fig, ax = plt.subplots(**kwargs)

    color_cycle = cycler(color=_NEW_TAB10)
    ax.set_prop_cycle(color_cycle)

    # could be done with one plot command, but...
    nv = ds.v.size
    for i in range(nv):
        ts_i = ds.isel(v=i)
        x = ts_i.x
        y = ts_i.y
        l, = ax.plot(x, y, lw=0.5, alpha=0.5)
        # highlight starting position
        ax.plot(x[0], y[0], 'o', c=l.get_color())


    ax.set(
        xlabel="$x$",
        ylabel="$y$",
        title="Vortons",
    )

    ax.set_aspect("equal", "box")

    fig.tight_layout()


# note much shared with vorton traj plot
def plot_tracer_trajectories(ds, **kwargs):
    """Plot tracer trajectories."""
    # select tracers
    it = ds.G == 0  # tracers boolean
    ds = ds.sel(v=it)

    fig, ax = plt.subplots(**kwargs)

    nv = ds.v.size
    for i in range(nv):
        ts_i = ds.isel(v=i)
        x = ts_i.x
        y = ts_i.y
        ax.plot(x, y, c="0.5", lw=0.5, alpha=0.5)

    ax.set(
        xlabel="$x$",
        ylabel="$y$",
        title="Tracers",
    )

    ax.set_aspect("equal", "box")

    fig.tight_layout()



# TODO
# def plot_ps
