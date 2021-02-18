"""
Model classes that act as drivers for the integration routines,
directing input and collecting output, etc.
"""
import abc
import copy
import glob
import os
import subprocess
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

from .plot import plot_vorton_trajectories, plot_tracer_trajectories, plot_ps
from .py import integrate_manual, integrate_scipy, MANUAL_STEPPERS, SCIPY_METHODS
from .vortons import Vortons, Tracers


def init_hist(  # TODO: could just be classmethod or staticmethod of ModelBase, maybe private
    n_vorton: int,
    n_time: int,  # in addition to t=0
    dt: float,
    #
    n_tracer=None,
    *,
    ds_attrs=None,
):
    """Create and return initial history `xr.Dataset`, with coordinates `'t'` (time) and `'v'` (vorton)."""

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
            "x": (("t", "v"), emp_tv(), {"long_name": "Vorton $x$ position"}),
            "y": (("t", "v"), emp_tv(), {"long_name": "Vorton $y4 position"}),
        },
        attrs=ds_attrs,
    )

    return ds


class ModelBase(abc.ABC):
    """Abstract base class for the models.
    Provides concrete methods `ModelBase.run` and `ModelBase.plot`, which work properly if the
    inheriting class defines a `_run` method that integrates the system and updates the `hist`.
    """
    def __init__(
        self,
        vortons: Vortons = None,
        tracers: Tracers = None,
        *,
        dt=0.1,
        nt=1000,
    ):
        """Set up vortons, tracers, etc., initialize history dataset using `init_hist`,..."""
        # vortons (default to equilateral triangle)
        if vortons is None:
            vortons = Vortons.regular_polygon(3, G=1)
        self.vortons = vortons
        """`vorts.vortons.Vortons` for the model instance."""
        self.vortons0 = copy.deepcopy(vortons)  # store initial state (in hist as well)
        """A copy of the initial `vorts.vortons.Vortons`, corresponding to the model input."""

        # tracers (leave if `None`)
        self.tracers = tracers
        """`vorts.vortons.Tracers` for the model instance, or `None` if no tracers."""
        self.tracers0 = copy.deepcopy(self.tracers)
        """A copy of the initial `vorts.vortons.Tracers` for the model instance, or `None` if no tracers."""

        # sim settings/parameters for every model
        self.dt = float(dt)
        r"""Time step $\delta t$ for the model output."""
        self.nt = int(nt)
        """The number of time steps to run for, such that `nt*dt` is the last time step."""
        self.nv = self.vortons0.n
        """Alias for `ModelBase.n_vortons`."""
        self.n_tracers = self.tracers0.n if self.tracers0 is not None else 0
        """The number of tracers."""
        self.n_vortons = self.nv
        """The number of vortons."""

        # initialize hist (an `xr.Dataset`)
        v0 = self.vortons0.maybe_with_tracers(self.tracers0)  # TODO: this maybe should be private
        # self.hist = init_hist(G, x0, y0, self.nv, self.nt, self.dt)
        self.hist = init_hist(self.nv + self.n_tracers, self.nt, self.dt)
        """An `xr.Dataset` with coordinates `'t'` (time) and `'v'` (vorton),
        created by `init_hist`.
        """
        self.hist["G"].loc[:] = v0.G  # G doesn't change during the sim
        # TODO: having to set the initial values this way is a bit awkward
        t_hist = self.hist.t
        self.hist["x"].loc[dict(t=t_hist[t_hist == 0])] = v0.x
        self.hist["y"].loc[dict(t=t_hist[t_hist == 0])] = v0.y

        # initially, the model hasn't been run
        self._has_run = False

    @abc.abstractmethod
    def _run(self):
        """The `_run` method in concrete classes should:

        1. Integrate the system from t=0 -> t=(nt*dt).

        2. Update `hist`.
        """
        ...

    def run(self):
        """Integrate and update the history dataset `hist`."""
        if self._has_run:
            warnings.warn("Note that the model has already been run.")
        self._run()
        self._has_run = True
        # TODO: with hist having been updated (presumably), update vortons?

        return self

    # might be better to use _plot_methods dict of name: function
    # so that certain models could extend the options
    def plot(self, which="vortons", **kwargs):
        """Plot results stored in the history data set `hist`.

        `**kwargs` are passed through to the corresponding plotting function.
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


FORT_BASE_DIR = Path(__file__).parent / "f" #.absolute()
# ^ the one that `bin`, `in`, `out`, `src` are in
assert (FORT_BASE_DIR / "src").exists()  # make sure that this is the right spot


def fort_bool(b: bool):
    """Convert Python boolean to a string of the corresponding Fortran logical (`.true.`/`.false.`)."""
    return ".true." if b else ".false."


class Model_f(ModelBase):
    """Thin wrapper for functionality of the Fortran model, whose source code is in `vorts/f/src/`.

    .. note::
       In this implementation we communicate with the Fortran program via text files.
    """
    # _allowed_int_scheme_names = ("FT", "RK4")

    def __init__(
        self,
        vortons: Vortons = None,
        tracers: Tracers = None,
        *,
        dt=0.1,
        nt=1000,
        # above are passed to base
        int_scheme_name='RK4',
        #
        write_vortons=True,  # maybe should put `flag` or something in these names
        write_tracers=False,
        write_ps=False,
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

            options: `'RK4'` (standard RK4; default), `'FT'` (1st-order forward Euler)

        """
        # call base initialization
        super().__init__(vortons, tracers, dt=dt, nt=nt)

        # other inputs
        self.int_scheme_name = int_scheme_name  # {'FT', 'RK4'}

        # output option flags
        self.f_write_out_vortons = write_vortons
        self.f_write_out_tracers = write_tracers
        self.f_write_out_ps = write_ps

        # executing the model
        self.vorts_exe_path = FORT_BASE_DIR / 'bin/vorts.exe'
        self.oe = ''  # we will store standard output and error here

        # write the text input files to directory `vorts/f/in`
        self.create_inputs()

    def create_inputs(self):
        """
        Create input files for the Fortran model
          describing the initial condition
          and the simulation settings.
        """
        # write vorton system initial state
        mat = self.vortons0.state_mat_full()  # needs to be rows of G, xi, yi
        np.savetxt(FORT_BASE_DIR / 'in/vorts_in.txt', mat,
                   delimiter=' ', fmt='%.16f', header='Gamma xi yi')

        # write tracers initial state (positions only)
        mat = self.tracers0.state_mat
        np.savetxt(FORT_BASE_DIR / 'in/tracers_in.txt', mat,
                   delimiter=' ', fmt='%.16f', header='xi, yi')

        # write model options
        mat = [
            self.dt,
            self.nt,
            self.int_scheme_name,
            fort_bool(self.f_write_out_vortons),
            fort_bool(self.f_write_out_tracers),
            fort_bool(self.f_write_out_ps),
        ]
        np.savetxt(FORT_BASE_DIR / 'in/vorts_sim_in.txt', mat,
                   delimiter=' ', fmt='%s')

    def _maybe_try_compile(self):
        """Try to run `make` if the executable is missing."""
        if not self.vorts_exe_path.exists():
            cwd = os.getcwd()
            os.chdir(FORT_BASE_DIR / "src")
            print(f"{self.vorts_exe_path!r} doesn't exist, attempting to `make`.\n")
            try:
                subprocess.run("make")
            except Exception as e:
                raise Exception(
                    f"Attempted `make` failed with exception (see above). "
                    "The Fortran code must be compiled before running!"
                ) from e
            finally:
                os.chdir(cwd)

    # implement abstract method `_run`
    def _run(self):
        """Invoke the Fortran model's executable and load the results."""
        self._maybe_try_compile()
        # exe_abs = str(self.vorts_exe_path)
        exe_rel = str(self.vorts_exe_path.relative_to(FORT_BASE_DIR))
        cmd = exe_rel
        # print(cmd)

        # invoke the Fortran model's executable
        cwd = os.getcwd()
        os.chdir(FORT_BASE_DIR)
        for f in glob.glob('./out/*'):  # non-hidden files
            os.remove(f)
        self.oe = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        os.chdir(cwd)
        # ^ hack for now, but could instead pass FORT_BASE_DIR into the Fortran program
        #   using vorts_sim_in.txt

        # load results from the text file outputs into `self.hist` (an `xr.Dataset`)
        self._load_results()

    def _load_results(self):
        """Load results from a run of the Fortran model."""
        sf = 0  # there may be blank line at end of the files?

        G = self.hist.G
        is_t = G == 0
        is_v = ~ is_t

        if self.f_write_out_vortons:
            vortons_file = FORT_BASE_DIR / 'out/vortons.csv'
            data = np.genfromtxt(vortons_file, delimiter=',', skip_header=1, skip_footer=sf)
            nrows = data.shape[0]
            i1 = np.arange(0, nrows-1, 2)
            i2 = np.arange(1, nrows, 2)
            self.hist["x"].loc[dict(v=is_v)] = data[i1, :].T  # need to swap dims because t is first in hist
            self.hist["y"].loc[dict(v=is_v)] = data[i2, :].T

        if self.f_write_out_tracers:
            tracers_file = FORT_BASE_DIR / 'out/tracers.csv'
            data = np.genfromtxt(tracers_file, delimiter=',', skip_header=1, skip_footer=sf)
            nrows = data.shape[0]
            i1 = np.arange(0, nrows-1, 2)
            i2 = np.arange(1, nrows, 2)
            self.hist["x"].loc[dict(v=is_t)] = data[i1, :].T
            self.hist["y"].loc[dict(v=is_t)] = data[i2, :].T

        # note: the ps code of the Fortran model only works for a specific case
        # (initial equi tri with the second point having x=0 and y>0)
        if self.f_write_out_ps:
            ps_file = FORT_BASE_DIR / "out/ps.txt"
            data = np.genfromtxt(ps_file, skip_header=1, skip_footer=sf)
            self.ps = data
