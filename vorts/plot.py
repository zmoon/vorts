"""
Plotting routines
"""
import inspect
import functools
import operator
import warnings

import cycler
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

def plot_vorton_trajectories(ds, title="Vortons", ax=None, **kwargs):
    """Plot lines: one for each vorton's trajectory.

    Parameters
    ----------
    ds : xarray.Dataset
        `hist` attribute of the model instance.
    title: str
        Plot title.
    **kwargs
        Passed on to [`plt.subplots()`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html)
    """
    # Select vortons
    iv = ds.G != 0  # TODO: move this to the ds creation in model
    ds = ds.sel(v=iv)

    fig, ax = _maybe_new_fig(ax=ax, **kwargs)

    color_cycle = cycler.cycler(color=_NEW_TAB10)
    ax.set_prop_cycle(color_cycle)

    # could be done with one plot command, but...
    nv = ds.v.size
    for i in range(nv):
        ts_i = ds.isel(v=i)
        x = ts_i.x
        y = ts_i.y
        l, = ax.plot(x, y, lw=0.5, alpha=0.5)
        # Highlight starting position
        ax.plot(x[0], y[0], 'o', c=l.get_color())

    _fig_post(fig, ax, title=title, frame="default")


# Note much shared with vorton traj plot
def plot_tracer_trajectories(ds, title="Tracers", ax=None, **kwargs):
    """Plot tracer trajectories.

    Parameters
    ----------
    ds : xarray.Dataset
        `hist` attribute of the model instance.
    title : str
        Plot title.
    **kwargs
        Passed on to [`plt.subplots()`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html)
    """
    # Select tracers
    it = ds.G == 0  # tracers boolean
    ds = ds.sel(v=it)

    fig, ax = _maybe_new_fig(ax=ax, **kwargs)

    nv = ds.v.size
    for i in range(nv):
        ts_i = ds.isel(v=i)
        x = ts_i.x
        y = ts_i.y
        ax.plot(x, y, c="0.5", lw=0.5, alpha=0.5)

    _fig_post(fig, ax, title=title, frame="default")


def select_poincare_times(ds, iv_ref=0, *, xtol=1e-2, ytol=1e-2):
    """From a model output dataset, extract data corresponding to times to use for the Poincaré section.

    We find the times when the reference vorton (with index `iv_ref`) is in its original position.

    Parameters
    ----------
    ds : xarray.Dataset
        For example, the `hist` attribute of the model instance.
    iv_ref : int
        Index of the vorton to use for reference.
    xtol, ytol : float
        Tolerance in each direction to use when searching for timesteps
        when the reference vorton is in its original position.

    Returns
    -------
    xarray.Dataset
        A selection (`.sel`) of the original `ds` consisting only of the timesteps
        when the reference vorton is approximately in the desired position.
    """
    ds_ref = ds.isel(v=iv_ref)  # selected vorton
    ds_ref0 = ds_ref.isel(t=0)  # initial position
    is_ps = (
        (np.abs(ds_ref.x - ds_ref0.x) <= xtol) &
        (np.abs(ds_ref.y - ds_ref0.y) <= ytol)
    )
    return ds.isel(t=is_ps.values)  # have to use `.values` since `v` dim does not match original


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
    plot_vortons=False,
    vorton_colors=None,
    ax=None,
    **kwargs
):
    """Poincaré section plot.

    Parameters
    ----------
    ds : xarray.Dataset
        Model output or output from `select_poincare_times`.
    iv_ref : int
        Index of the vorton to use for reference (passed to `select_poincare_times`).
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
    frame : str, {'none', 'only', 'default'}
        No frame (using `remove_frame`), frame only (using `frame_only`), or default style (do nothing).
    plot_vortons : bool
        Whether to plot the initial vorton positions on top.
    vorton_colors
        Colors, of valid format (like `c`). OR single color.
        By default, cycles through new tab10.
    ax : matplotlib.axes.Axes, optional
        Optionally pass `ax` on which to plot. Otherwise a new figure will be created.
    **kwargs
        Passed on to `select_poincare_times` (if applicable)
        or [`plt.subplots()`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html)
        or [`ax.plot()`](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.plot.html)

    See also
    --------
    select_poincare_times
    """
    # TODO: algo for choosing marker size and alpha based on total number of points (but also allow passing in)

    # Separate kwargs
    select_poincare_times_kwargs = {k: kwargs.pop(k) for k in inspect.getfullargspec(select_poincare_times).kwonlyargs if k in kwargs}
    subplots_kwargs = {k: kwargs.pop(k) for k in _allowed_subplots_kwargs if k in kwargs}
    plot_kwargs = {k: kwargs.pop(k) for k in _allowed_plot_kwargs if k in kwargs}

    # Subset data to approximate Poincare section
    ds = select_poincare_times(ds, iv_ref, **select_poincare_times_kwargs)

    # Select tracers
    it = ds.G == 0
    ds_v0 = ds.sel(v=~it).isel(t=0)
    ds = ds.sel(v=it)

    fig, ax = _maybe_new_fig(ax=ax, **subplots_kwargs)

    # Set up property cycling
    cyclers = [
        cycler.cycler(color=_to_list(c)),
        cycler.cycler(markersize=_to_list(ms)),
        cycler.cycler(alpha=_to_list(alpha)),
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

    # Plot points. Other marker attributes are included in the property cycle.
    ax.plot(x, y, ".", mew=0, **plot_kwargs)

    # Plot vortons?
    if plot_vortons:
        if not vorton_colors:
            vorton_colors = _NEW_TAB10
        for x_v, y_v, c_v in zip(ds_v0.x.values, ds_v0.y.values, cycler.cycle(_to_list(vorton_colors))):
            ax.plot(x_v, y_v, "o", c=c_v, ms=10, alpha=1)

    # Set labels, title, frame settings
    _fig_post(fig, ax, title=title, frame=frame)


def frame_only(ax=None, *, keep_ax_labels=True, keep_title=True):
    """Remove ticks and tick labels from `ax` (uses current by default)."""
    if ax is None:
        ax = plt.gca()

    # Remove all ticks and tick labels
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

def _maybe_new_fig(ax=None, **kwargs):
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
