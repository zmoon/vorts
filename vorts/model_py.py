"""
Driver for the Python version.
"""

from glob import glob
import os
import subprocess
from typing import List, Tuple
import warnings

#import matplotlib.pyplot as plt
import numpy as np

from .py import (integrate_manual, integrate_scipy, MANUAL_STEPPERS, SCIPY_METHODS, calc_C)
from .vortons import Vorton


class model_py:  # TODO: model base class?
    """Model in Python."""
    _manual_steppers = MANUAL_STEPPERS
    _scipy_methods = SCIPY_METHODS

    def __init__(
        self,
        G,
        xi,
        yi,
        *,
        dt=0.1,
        nt=1000,
        int_scheme_name='RK4_3',
        **int_scheme_kwargs,
    ):
        """Create model with given settings.

        **int_scheme_kwargs
            see signatures of `integrate_manual` and `integrate_scipy`
        """

        # vorton IC arrays
        self.G_vals = G
        self.xi_vals = xi
        self.yi_vals = yi

        # sim settings
        self.dt = float(dt)
        self.nt = nt
        self.int_scheme_name = int_scheme_name
        self.int_scheme_kwargs = int_scheme_kwargs

        self._allowed_int_scheme_names = list(self._manual_steppers) + list(self._scipy_methods)
        if self.int_scheme_name not in self._allowed_int_scheme_names:
            raise ValueError(
                f"{self.int_scheme_name!r} is not one of the allowed options for `int_scheme_name`:\n"
                f"{self._allowed_int_scheme_names}"
            )

        # create vortons
        self.vortons = []
        for Gxy in zip(self.G_vals, self.xi_vals, self.yi_vals):
            self.vortons.append(Vorton(*Gxy, nt))

        self.l = 0  # time step index

        # for adaptive time stepping calculations
        self.C_0 = calc_C(self.vortons, self.l)


    def run(self):
        """Integrate from 0 to end (nt*dt)."""

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
            integrate_scipy(
                self.vortons,
                t_eval,
                self.G_vals,
                method=self._scipy_methods[self.int_scheme_name],
                max_step=self.dt,
                **self.int_scheme_kwargs
            )
