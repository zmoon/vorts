"""
Driver for the Python version.
"""
# from typing import List, Tuple
# import warnings

#import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from ._model import ModelBase
from .py import (integrate_manual, integrate_scipy, MANUAL_STEPPERS, SCIPY_METHODS)
from .vortons import Vortons, Tracers


class Model_py(ModelBase):
    """Model in Python."""
    _manual_steppers = MANUAL_STEPPERS
    _scipy_methods = SCIPY_METHODS
    _allowed_int_scheme_names = list(_manual_steppers) + list(_scipy_methods)
    # _allowed_int_scheme_names = list(_scipy_methods)  # temporarily disable manual steppers methods

    def __init__(
        self,
        vortons: Vortons = None,
        tracers: Tracers = None,
        *,
        dt=0.1,
        nt=1000,
        # above are passed to base
        int_scheme_name='scipy_RK45',
        **int_scheme_kwargs,
    ):
        r"""

        Parameters
        ----------
        vortons : vorts.vortons.Vortons
            default: equilateral triangle with inscribing circle radius of $1$ and all $G=1$.

        tracers : vorts.vortons.Tracers
            default: no tracers

        dt : float
            Time step $\delta t$ for the output.
            Additionally, for the integrators, `dt` is used as the constant or maximum integration time step
            depending on the integration scheme.
        nt : int
            Number of time steps to run (not including $t=0$).

        int_scheme_name : str
            Time integration scheme name.

            default: `'scipy_RK45'` (the default method of [`scipy.integrate.solve_ivp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html))

            Other currently valid options are:

            * `'scipy_DOP853'` -- also from SciPy, higher-order RK method `'DOP853'` is written in Python like `'RK45'`

            * `'FT'` -- handwritten 1st-order forward Euler

            * `'RK4'` -- handwritten standard RK4

        **int_scheme_kwargs
            Passed on to `vorts.py.integ.integrate_manual` or `vorts.py.integ.integrate_scipy`.

        """
        # call base initialization
        super().__init__(vortons, tracers, dt=dt, nt=nt)

        # other inputs
        self.int_scheme_name = int_scheme_name
        self.int_scheme_kwargs = int_scheme_kwargs

        # check `int_scheme_name`
        if self.int_scheme_name not in self._allowed_int_scheme_names:
            raise ValueError(
                f"{self.int_scheme_name!r} is not one of the allowed options for `int_scheme_name`:\n"
                f"{self._allowed_int_scheme_names}"
            )

        # calculate initial C, used for adaptive time stepping tolerance checks
        self.C_0 = self.vortons0.C()

    # implement abstract method `_run`
    def _run(self):
        dt, nt = self.dt, self.nt
        # t_eval = np.arange(dt, (nt+1)*dt, dt)
        t_eval = np.arange(1, nt+1)*dt  # could start at 0?

        # manual (handwritten) integrators
        if "scipy" not in self.int_scheme_name:
            v0 = self.vortons0.maybe_with_tracers(self.tracers0)
            x0 = v0.x
            y0 = v0.y
            G = v0.G
            xhist, yhist = integrate_manual(
                G,
                x0,
                y0,
                #
                self.C_0,
                t_eval,
                stepper=self._manual_steppers[self.int_scheme_name],
                **self.int_scheme_kwargs
            )
            # returned data have shape (nv, nt)
            nv = v0.n
            t = self.hist.t
            self.hist["x"].loc[dict(t=t[t > 0])] = xhist.T
            self.hist["y"].loc[dict(t=t[t > 0])] = yhist.T

        # integration using SciPy
        else:
            v0 = self.vortons0.maybe_with_tracers(self.tracers0)
            y0 = v0.state_vec()
            G_col = v0.G_col
            data = integrate_scipy(
                y0,
                t_eval,
                G_col,
                #
                method=self._scipy_methods[self.int_scheme_name],
                max_step=dt,
                **self.int_scheme_kwargs
            )
            # returned data has shape (2nv, nt), where n is number of vortons and nt number of time steps
            nv = v0.n
            t = self.hist.t
            self.hist["x"].loc[dict(t=t[t > 0])] = data[:nv, :].T  # need to swap dims because t is first in hist
            self.hist["y"].loc[dict(t=t[t > 0])] = data[nv:, :].T
