"""
Base class for models (classes which wrap functions that implement the integration).
"""
import abc
import copy
import warnings

import numpy as np
import xarray as xr

from .plot import plot_vorton_trajectories, plot_tracer_trajectories, plot_ps
from .vortons import Vortons, Tracers


def init_hist(
    n_vorton: int,
    n_time: int,  # in addition to t=0
    dt: float,
    #
    n_tracer=None,
    *,
    ds_attrs=None,
):
    """Create initial history `xr.Dataset`."""

    # if n_tracer is not None:
    if ds_attrs is None:
        ds_attrs = {}

    t = np.arange(0, n_time+1)*dt
    nt = t.size

    v = np.arange(0, n_vorton)
    nv = n_vorton

    def emp_v():
        return np.empty((nv,))

    def emp_tv():
        return np.empty((nt, nv))

    ds = xr.Dataset(
        coords={
            "t": ("t", t, {"long_name": "unitless time"}),
            "v": ("v", v, {"long_name": "vorton num"}),
        },
        data_vars={
            "G": (("v",), emp_v(), {"long_name": r"Vorton strength $\Gamma$ (circulation)"}),
            "x": (("t", "v"), emp_tv(), {"long_name": "Vorton x position"}),
            "y": (("t", "v"), emp_tv(), {"long_name": "Vorton y position"}),
        },
        attrs=ds_attrs,
    )

    return ds


class ModelBase(abc.ABC):
    """Base class for the models."""
    def __init__(
        self,
        vortons: Vortons = None,
        tracers: Tracers = None,
        *,
        dt=0.1,
        nt=1000,
    ):
        """Set up vortons, tracers, etc., initialize hist..."""
        # vortons (default to equilateral triangle)
        if vortons is None:
            vortons = Vortons.regular_polygon(3, G=1)
        self.vortons = vortons
        self.vortons0 = copy.deepcopy(vortons)  # store initial state (in hist as well)

        # tracers (leave if `None`)
        self.tracers = tracers
        self.tracers0 = copy.deepcopy(self.tracers)

        # sim settings/parameters for every model
        self.dt = float(dt)
        self.nt = int(nt)
        self.nv = self.vortons0.n
        self.n_tracers = self.tracers0.n if self.tracers0 is not None else 0
        self.n_vortons = self.nv

        # initialize hist (an `xr.Dataset`)
        v0 = self.vortons0.maybe_with_tracers(self.tracers0)
        # self.hist = init_hist(G, x0, y0, self.nv, self.nt, self.dt)
        self.hist = init_hist(self.nv + self.n_tracers, self.nt, self.dt)
        self.hist["G"].loc[:] = v0.G  # G doesn't change during the sim
        # TODO: having to set the initial values this way is a bit awkward
        t_hist = self.hist.t
        self.hist["x"].loc[dict(t=t_hist[t_hist == 0])] = v0.x
        self.hist["y"].loc[dict(t=t_hist[t_hist == 0])] = v0.y

        # initially, the model hasn't been run
        self._has_run = False

    @abc.abstractmethod
    def _run(self):
        """`_run` method should (1) integrate the system 0->(nt*dt) and (2) update `hist`."""
        ...

    def run(self):
        """Integrate and update hist."""
        if self._has_run:
            warnings.warn("Note that the model has already been run.")
        self._run()
        self._has_run = True
        # TODO: with hist having been updated (presumably), update vortons?

    # might be better to use _plot_methods dict of name: function
    # so that certain models could extend the options
    def plot(self, which="vortons", **kwargs):
        """Plot results stored in history data set `hist`.

        `**kwargs` passed through to the corresponding plotting function.
        """
        if not self._has_run:
            raise Exception("The model has not yet been run.")

        if which == "vortons":
            plot_vorton_trajectories(self.hist, **kwargs)

        elif which == "tracers":
            plot_tracer_trajectories(self.hist, **kwargs)

        elif which == "poincare":
            plot_ps(self.hist, **kwargs)

        else:
            raise NotImplementedError(f"which={which!r}")
