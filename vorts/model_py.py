"""
Driver for the Python version.
"""
import copy
# from glob import glob
# import os
# import subprocess
# from typing import List, Tuple
import warnings

#import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from .plot import plot_vorton_trajectories, plot_tracer_trajectories
from .py import (integrate_manual, integrate_scipy, MANUAL_STEPPERS, SCIPY_METHODS)
from .vortons import Vortons



# model_f could use it when collecting the data
# until the Fortran model writes nc files
# or could be moved to a different module
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
            "G": (("v",), emp_v(), {"long_name": "Vorton strength $\Gamma$ (circulation)"}),
            "x": (("t", "v"), emp_tv(), {"long_name": "Vorton x position"}),
            "y": (("t", "v"), emp_tv(), {"long_name": "Vorton y position"}),
        },
        attrs=ds_attrs,
    )

    return ds


class model_py:  # TODO: model base class?
    """Model in Python."""
    _manual_steppers = MANUAL_STEPPERS
    _scipy_methods = SCIPY_METHODS
    # _allowed_int_scheme_names = list(_manual_steppers) + list(_scipy_methods)
    _allowed_int_scheme_names = list(_scipy_methods)

    def __init__(
        self,
        vortons: Vortons = None,
        *,
        dt=0.1,
        nt=1000,
        int_scheme_name='scipy_RK45',
        **int_scheme_kwargs,
    ):
        """Create model with given settings.

        Parameters
        ----------
        vortons : Vortons
            default: equilateral triangle with all G=1

        dt : float
            time step for the output
            for the integrators, used as the constant integration time step
            or as the maximum integration time step
        nt : int
            number of time steps to run (not including t=0)
        int_scheme_name : str
            default: 'RK4_3' (handwritten basic RK4 stepper)

        **int_scheme_kwargs
            passed on to `integrate_manual` or `integrate_scipy`
            see their signatures
        """
        # vortons
        if vortons is None:
            vortons = Vortons.regular_polygon(3, G=1)
        self.vortons = vortons  # TODO: update these after the run!
        self.vortons0 = copy.deepcopy(self.vortons)  # store initial state

        # sim settings
        self.dt = float(dt)
        self.nt = int(nt)
        self.nv = self.vortons0.n
        self.int_scheme_name = int_scheme_name
        self.int_scheme_kwargs = int_scheme_kwargs

        # check `int_scheme_name`
        if self.int_scheme_name not in self._allowed_int_scheme_names:
            raise ValueError(
                f"{self.int_scheme_name!r} is not one of the allowed options for `int_scheme_name`:\n"
                f"{self._allowed_int_scheme_names}"
            )

        # self.l = 0  # time step index (doesn't seem to be used)

        # for adaptive time stepping calculations
        self.C_0 = self.vortons0.C()

        # initialize hist (xr.Dataset)
        # self.hist = init_hist(G, x0, y0, self.nv, self.nt, self.dt)
        self.hist = init_hist(self.nv, self.nt, self.dt)
        self.hist["G"].loc[:] = self.vortons0.G
        # TODO: having to set the initial values this way is a bit awkward
        t_hist = self.hist.t
        self.hist["x"].loc[dict(t=t_hist[t_hist == 0])] = self.vortons0.x
        self.hist["y"].loc[dict(t=t_hist[t_hist == 0])] = self.vortons0.y

        self._has_run = False


    def run(self):
        """Integrate from 0 to end (nt*dt)."""
        if self._has_run:
            warnings.warn("The model has already been run.")

        dt, nt = self.dt, self.nt
        # t_eval = np.arange(dt, (nt+1)*dt, dt)
        t_eval = np.arange(1, nt+1)*dt

        if "scipy" not in self.int_scheme_name:
            integrate_manual(
                self.vortons,
                self.C_0,
                t_eval,
                self.G_vals,
                stepper=self._manual_steppers[self.int_scheme_name],
                **self.int_scheme_kwargs
            )
        else:  # Scipy integrator
            data = integrate_scipy(
                self.vortons0.state_vec(),
                t_eval,
                self.vortons0.G_col,
                #
                method=self._scipy_methods[self.int_scheme_name],
                max_step=dt,
                **self.int_scheme_kwargs
            )
            # returned data has shape (2nv, nt), where n is number of vortons and nt number of time steps
            nv = self.nv
            t_hist = self.hist.t
            self.hist["x"].loc[dict(t=t_hist[t_hist > 0])] = data[:nv, :].T  # need to swap dims because t is first in hist
            self.hist["y"].loc[dict(t=t_hist[t_hist > 0])] = data[nv:, :].T

        self._has_run = True


    def plot(self, which="vortons", **kwargs):
        """Plot.

        **kwargs passed through to the corresponding plotting function.
        """
        if not self._has_run:
            raise Exception("The model has not yet been run.")

        if which == "vortons":
            plot_vorton_trajectories(self.hist, **kwargs)

        elif which == "tracers":
            plot_tracer_trajectories(self.hist, **kwargs)

        else:
            raise NotImplementedError(f"which={which!r}")
