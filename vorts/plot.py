"""
Plotting routines
"""
import inspect
import functools
import operator
import warnings

from cycler import cycler
import matplotlib.pyplot as plt
import numpy as np


# Tableau's newer version of tab10
# https://www.tableau.com/about/blog/2016/7/colors-upgrade-tableau-10-56782
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

# TODO: PlotMethods object that could allow .plot_type(...) and ("plot_type", ...) and be attached to the Model? (like pandas/xarray)

# TODO: should plotters return `fig, ax`, or just `ax`? or something else? xarray returns the set of matplotlib artists

# TODO: routine to determine system rotation; plot trajectories with respect to this rotating ref frame

def plot_vorton_trajectories(ds, ax=None, **kwargs):
    """Plot lines: one for each vorton's trajectory.

    Parameters
    ----------
    ds : xarray.Dataset
        `hist` attribute of the model instance.
    **kwargs
        Passed on to [`plt.subplots()`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html)
    """
    # select vortons
    iv = ds.G != 0
    ds = ds.sel(v=iv)

    fig, ax = maybe_new_figure(ax=ax, **kwargs)

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
def plot_tracer_trajectories(ds, ax=None, **kwargs):
    """Plot tracer trajectories.

    Parameters
    ----------
    ds : xarray.Dataset
        `hist` attribute of the model instance.
    **kwargs
        Passed on to [`plt.subplots()`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html)
    """
    # select tracers
    it = ds.G == 0  # tracers boolean
    ds = ds.sel(v=it)

    fig, ax = maybe_new_figure(ax=ax, **kwargs)

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


def ps_data(ds, iv_ref=0, *, xtol=1e-2):
    """From full set of data, extract data corresponding to times for the Poincare Section.

    We find the times when the reference vorton is in a certain place
    or passing through a certain plane (TODO).

    Parameters
    ----------
    ds : xarray.Dataset
        `hist` attribute of the model instance.
    iv_ref : int
        Index of the vorton to use for reference.
    xtol : float
        Tolerance to use when searching for timesteps when the reference vorton is in its original position.

    Returns
    -------
    xarray.Dataset
        Only consisting of the timesteps when the reference vorton is approximately in the desired position.
    """
    # initial position
    r0_ref = ds.isel(t=0, v=iv_ref)
    x0 = r0_ref.x#.values
    # y0 = r0_ref.y#.values

    # xtol = 1e-2
    # ytol = 1e-2

    # TODO: need to be more careful. should make wrt. center-of-vort, ...
    # cond = (np.abs(ds.x - x0) <= xtol) & (np.abs(ds.y - y0) <= ytol) & (ds.v == iv_ref)
    # cond = (np.abs(ds.x - x0) <= xtol) & (ds.v == iv_ref)
    cond = (np.abs(ds.x - x0) <= xtol) & (ds.y > 0) & (ds.v == iv_ref)
    ds_ps_ref = ds.where(cond, drop=True)
    t_ps = ds_ps_ref.t  # ps times

    ds_ps = ds.sel(t=t_ps)

    return ds_ps


_allowed_subplots_kwargs = ("figsize", "linewidth", "frameon", "tight_layout", "constrained_layout")
_allowed_plot_kwargs = ("c", "color", "ms", "markersize", "alpha")


def plot_ps(ds, *,
    iv_ref=0,
    c="0.35",
    ms=0.2,
    alpha=0.5,
    cycle_comp="add",
    cycle_by="vorton",
    title="Poincaré section (tracers)",
    frame="none",
    ax=None,
    **kwargs
):
    """Poincaré section plot.

    Here using the data set of all data.

    Parameters
    ----------
    ds : xarray.Dataset
        Model output or output from `ps_data`.
    iv_ref : int
        Index of the vorton to use for reference (passed to `ps_data`).
    c : str or array_like
        Marker color (any of [these formats](https://matplotlib.org/stable/tutorials/colors/colors.html)).
        OR iterable of individual colors, which will be cycled.
    ms : float
        Marker size.
        OR iterable of individual sizes, which will be cycled.
    alpha : float
        Marker alpha.
        OR iterable of individual alpha values, which will be cycled.
    cycle_comp : str, {'add', 'multiply'}
        How to compose the property cyclers---[add](https://matplotlib.org/cycler/#addition)
        or [multiply](https://matplotlib.org/cycler/#integer-multiplication).
    cycle_by : str, {'vorton', 'time'}
        Whether certain times will have certain properties, or certain vortons will.
    title : str
        Plot title.
    frame : str {'none', 'only', 'default'}
        No frame (using `remove_frame`), frame only (using `frame_only`), or default style (do nothing).
    ax : matplotlib.axes.Axes, optional
        Optionally pass `ax` on which to plot. Otherwise a new figure will be created.
    **kwargs
        Passed on to `ps_data` (if applicable)
        or [`plt.subplots()`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html)
        or [`ax.plot()`](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.plot.html)

    See also
    --------
    ps_data
    """
    # Separate kwargs
    ps_data_kwargs = {k: kwargs.pop(k) for k in inspect.getfullargspec(ps_data).kwonlyargs if k in kwargs}
    subplots_kwargs = {k: kwargs.pop(k) for k in _allowed_subplots_kwargs if k in kwargs}
    plot_kwargs = {k: kwargs.pop(k) for k in _allowed_plot_kwargs if k in kwargs}

    # Subset data to approximate Poincare section
    ds = ps_data(ds, iv_ref, **ps_data_kwargs)

    # select tracers
    it = ds.G == 0
    ds = ds.sel(v=it)

    fig, ax = maybe_new_figure(ax=ax, **subplots_kwargs)

    # TODO: (optionally?) plot vorton initial positions / positions at reference time?

    # TODO: algo for choosing marker size and alpha (but also allow passing in)

    # Set up property cycling
    cyclers = [
        cycler(color=_to_list(c)),
        cycler(markersize=_to_list(ms)),
        cycler(alpha=_to_list(alpha)),
    ]
    if cycle_comp == "add":
        cyclers = _reconcile_cyclers_for_adding(cyclers)  # make all the same length
        cycle = sum(cyclers[1:], start=cyclers[0])
    elif cycle_comp == "multiply":
        cycle = functools.reduce(operator.mul, cyclers[1:], cyclers[0])
    else:
        raise ValueError("invalid value for `cycle_comp`: {cycle_comp!r}")
    ax.set_prop_cycle(cycle)

    # If not cycling properties, it should be more efficient to plot the flattened data.
    if len(cycle) > 1:
        if cycle_by == "vorton":
            x = ds.x.values  # (nt, nv), columns are vorton time series, each plotted separately
            y = ds.y.values
            if ds.v.size < len(cycle):
                warnings.warn(
                    f"cycle length {len(cycle)} is longer than the number of vortons+tracers ({ds.v.size}).",
                )
        elif cycle_by == "time":
            x = ds.x.values.T  # (nv, nt), columns are times, each plotted separately
            y = ds.y.values.T
            if ds.t.size < len(cycle):
                warnings.warn(
                    f"cycle length {len(cycle)} is longer than the number of times ({ds.t.size}).",
                )
        else:
            raise ValueError
    else:
        x = ds.x.values.ravel()
        y = ds.y.values.ravel()

    # Plot points. Other marker attributes are included in the cycler.
    ax.plot(x, y, ".", mew=0, **plot_kwargs)

    # Set labels, title, frame settings
    _fig_post(fig, ax, title=title, frame=frame)


def frame_only(ax=None, *, keep_ax_labels=True, keep_title=True):
    """Remove ticks and tick labels from `ax` (uses current by default)."""
    if ax is None:
        ax = plt.gca()

    ax.tick_params(
        axis="both",
        which="both",  # 'major', 'minor', or 'both'
        bottom=False,
        labelbottom=False,
        left=False,
        labelleft=False,
        top=False,
        labeltop=False,
        right=False,
        labelright=False,
    )

    if not keep_ax_labels:
        ax.set(
            xlabel="",
            ylabel="",
        )

    if not keep_title:
        ax.set_title("")


def remove_frame(ax=None, *, keep_title=True):
    """Remove ticks, tick labels, and frame (spines) from `ax` (uses current by default)."""
    if ax is None:
        ax = plt.gca()

    frame_only(ax, keep_ax_labels=False, keep_title=keep_title)

    ax.set(
        frame_on=False,
    )

def maybe_new_figure(ax=None, **kwargs):
    """Return `fig, ax`, both new if `ax` is `None`. `**kwargs` passed to `plt.subplots()`."""
    if ax is None:
        fig, ax = plt.subplots(**kwargs)
    else:
        fig = ax.get_figure()

    return fig, ax


def _is_iterable(x):
    """Check if `x` is iterable via duck typing."""
    try:
        _ = iter(x)
    except TypeError:
        return False
    else:
        return True


def _to_list(v):
    """Convert value(s) to list. If `v` is `str`, preserve it instead of returning a list of chars."""
    return list(v) if (
        _is_iterable(v) and not isinstance(v, (str,))
    ) else [v]


def _reconcile_cyclers_for_adding(cyclers):
    """Make all cyclers in list `cyclers` the same length."""
    max_len = max(len(c) for c in cyclers)
    ret = []
    for c in cyclers:
        if not max_len % len(c) == 0:
            raise ValueError(
                "all property set lengths must be even multiples of the longest length, "
                f"currently {max_len}"
            )
        ret.append(c * (max_len // len(c)))
    return ret


def _fig_post(fig, ax, *, title, frame, **kwargs):
    """Set axis labels, title, aspect equal, frame settings, ..."""
    ax.set_title(title)
    frame = frame.lower()
    if frame == "only":
        frame_only(ax=ax, **kwargs)
    elif frame in ("none", "remove"):
        remove_frame(ax=ax, **kwargs)
    elif frame == "default":
        ax.set(
            xlabel="$x$",
            ylabel="$y$",
        )
    else:
        raise ValueError("invalid value for `frame`: {frame!r}")

    ax.set_aspect("equal", "box")
    fig.set_tight_layout(True)
