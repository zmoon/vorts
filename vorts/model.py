"""
Model classes that act as drivers for the integration routines,
directing input and collecting output, etc.
"""

import abc
import contextlib
import copy
import glob
import os
import platform
import subprocess
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

from .plot import _PlotMethods
from .py import MANUAL_STEPPERS, SCIPY_METHODS, integrate_manual, integrate_scipy
from .vortons import Tracers, Vortons


class ModelBase(abc.ABC):
    """Abstract base class for the models.

    Provides concrete methods `ModelBase.run` and `ModelBase.plot`, which work properly if the
    inheriting class defines a `_run` method that integrates the system and updates `ModelBase.hist`,
    populating it with the run's output data.
    """

    def __init__(
        self,
        vortons: Vortons = None,
        tracers: Tracers = None,
        *,
        dt=0.1,
        nt=1000,
    ):
        """
        * `vortons` default to equilateral triangle (`vorts.vortons.Vortons.regular_polygon` with `3, G=1`)
        * `tracers` default to `None` (no tracers used)
        """
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
        """The number of time steps to run for, such that `t=nt*dt` is the last time simulated
        and we have a total of `nt+1` times when including the initial state.
        """
        self.n_tracers = self.tracers0.n if self.tracers0 is not None else 0
        """The number of tracers."""
        self.n_vortons = self.vortons0.n
        """The number of vortons (tracers not included!)."""
        self.n_points = self.n_vortons + self.n_tracers
        """The number of vortons + tracers."""
        self.n_timesteps = self.nt
        """Alias for `ModelBase.nt`."""

        # combine vortons and tracers (used by some models to feed combined initial states to integrators)
        self._vt0 = self.vortons0._maybe_add_tracers(self.tracers0)

        # initially no history
        self.hist = None
        """An [`xr.Dataset`](http://xarray.pydata.org/en/stable/generated/xarray.Dataset.html)
        with coordinates `'t'` (time) and `'v'` (vorton index, including any tracers).
        Equal to `None` before the model has run.
        """

        # initially, the model hasn't been run
        self._has_run = False

        # attach _PlotMethods instance
        self.plot = _PlotMethods(self)

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

    def _res_to_xr(self, xhist, yhist):
        """Take full trajectory histories `xhist` and `yhist`
        in $(t, v)$ dim order and create a model results `xr.dataset`.
        """
        vt0 = self._vt0

        G = vt0.G
        t = np.arange(0, self.nt + 1) * self.dt
        v = np.arange(0, vt0.n)

        is_t = G == 0

        ds = xr.Dataset(
            coords={
                "t": ("t", t, {"long_name": "Unitless elapsed time"}),
                "v": (
                    "v",
                    v,
                    {
                        "long_name": "Vorton index",
                        "description": (
                            r"This includes true vortons (with nonzero $\Gamma$, coming first) "
                            r"and may also include tracers ($\Gamma = 0$, coming after)."
                        ),
                    },
                ),
            },
            data_vars={
                "G": (("v",), G, {"long_name": r"Vorton strength $\Gamma$ (circulation)"}),
                "x": (("t", "v"), xhist, {"long_name": "Vorton $x$ position"}),
                "y": (("t", "v"), yhist, {"long_name": "Vorton $y$ position"}),
                "is_t": (
                    ("v",),
                    is_t,
                    {
                        "long_name": "Tracer mask",
                        "description": (
                            r"Vortons have nonzero $\Gamma$ values. Tracers have no $\Gamma$."
                        ),
                    },
                ),
            },
            attrs={
                "model_class": self.__class__.__name__,
                "int_scheme_name": getattr(self, "int_scheme_name", "?"),
            },
        )
        return ds


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
        int_scheme_name="RK4",
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

            default: `'RK4'` -- handwritten standard RK4, using Numba to calculate the tendencies

            Other currently valid options are:

            * `'scipy_*'` methods -- use [`scipy.integrate.solve_ivp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)

                where the `*` can be `RK45`, `DOP854`, `Radau`, `BDF`, `LSODA`

            * `'FT'` -- handwritten 1st-order forward Euler

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
        t_eval = np.arange(0, nt + 1) * dt
        v0 = self._vt0

        # manual (handwritten) integrators
        if "scipy" not in self.int_scheme_name:
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
                **self.int_scheme_kwargs,
            )
            # returned data have shape (nv, nt)
            self.hist = self._res_to_xr(xhist.T, yhist.T)

        # integration using SciPy
        else:
            y0 = v0.xy.T.flatten()  # needs to be a 1-d array!
            G_col = v0.G_col
            data = integrate_scipy(
                y0,
                t_eval,
                G_col,
                #
                method=self._scipy_methods[self.int_scheme_name],
                max_step=dt,
                **self.int_scheme_kwargs,
            )
            # returned data has shape (2nv, nt), where n is number of vortons and nt number of time steps
            nv = v0.n
            self.hist = self._res_to_xr(data[:nv, :].T, data[nv:, :].T)


FORT_BASE_DIR = Path(__file__).parent / "f"  # .absolute()
# ^ the one that `bin`, `in`, `out`, `src` are in


def fort_bool(b: bool):
    """Convert Python boolean to a string of the corresponding Fortran logical (`.true.`/`.false.`)."""
    return ".true." if b else ".false."


@contextlib.contextmanager
def _out_and_back(p):
    """Context manager: change working directory to path `p` but return to original cwd after."""
    # ref: https://stackoverflow.com/a/170174
    cwd = os.getcwd()
    os.chdir(p)
    try:
        yield
    finally:
        os.chdir(cwd)


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
        int_scheme_name="RK4",
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
        assert (FORT_BASE_DIR / "src").exists(), "Fortran source directory not found!"
        exe_name = "vorts.exe" if platform.system() == "Windows" else "vorts"
        self.vorts_exe_path = FORT_BASE_DIR / "bin" / exe_name
        self.oe = ""  # we will store standard output and error here

        # create the needed subdirs if they don't exist (pip install)
        for sub in ("bin", "in", "out"):
            p = FORT_BASE_DIR / sub
            if not p.exists():
                p.mkdir()

        # write the text input files to directory `vorts/f/in`
        self.create_inputs()

    def create_inputs(self):
        """
        Create input files for the Fortran model
          describing the initial condition
          and the simulation settings.
        """
        # # write vorton system initial state
        # mat = self.vortons0.state_mat_full()  # needs to be rows of G, xi, yi
        # np.savetxt(FORT_BASE_DIR / 'in/vorts_in.txt', mat,
        #            delimiter=' ', fmt='%.16f', header='Gamma xi yi')

        # # write tracers initial state (positions only)
        # mat = self.tracers0.state_mat
        # np.savetxt(FORT_BASE_DIR / 'in/tracers_in.txt', mat,
        #            delimiter=' ', fmt='%.16f', header='xi yi')

        # Write combined vortons + tracers, with vortons first
        mat = self._vt0.state_mat_full()
        np.savetxt(
            FORT_BASE_DIR / "in/vortons.txt", mat, delimiter=" ", fmt="%.16f", header="Gamma xi yi"
        )

        # write model options
        mat = [
            self.dt,
            self.nt,
            self.int_scheme_name,
            fort_bool(self.f_write_out_vortons),
            fort_bool(self.f_write_out_tracers),
            fort_bool(self.f_write_out_ps),
        ]
        np.savetxt(FORT_BASE_DIR / "in/settings.txt", mat, delimiter=" ", fmt="%s")

    def _maybe_try_compile(self):
        """Try to run `make` if the executable is missing."""
        if not self.vorts_exe_path.exists():
            with _out_and_back(FORT_BASE_DIR / "src"):
                print(f"{self.vorts_exe_path!r} doesn't exist, attempting to `make`.\n")
                try:
                    subprocess.run("make")
                except Exception as e:
                    raise Exception(
                        "Attempted `make` failed with exception (see above). "
                        "The Fortran code must be compiled before running!"
                    ) from e

    # implement abstract method `_run`
    def _run(self):
        """Invoke the Fortran model's executable and load the results."""
        self._maybe_try_compile()

        # Invoke the Fortran model's executable
        # Note: could instead pass FORT_BASE_DIR into the Fortran program so could run from anywhere
        cmd = str(self.vorts_exe_path.relative_to(FORT_BASE_DIR))
        with _out_and_back(FORT_BASE_DIR):
            for f in glob.glob("./out/*"):  # non-hidden files
                os.remove(f)
            self.oe = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            # Note: could instead use `error stop` in the Fortran to get non-zero exit code and `CalledProcessError`
            s_oe = str(self.oe, "utf-8")
            if s_oe.startswith("STOP"):
                raise Exception(f"the model stopped early, with message:\n{str(s_oe)}")

        # Load results from the text file outputs into `self.hist` (an `xr.Dataset`)
        self._load_results()

    def _load_results(self):
        """Load results from a run of the Fortran model."""
        sf = 0  # there may be blank line at end of the files?

        nt = self.nt  # time steps taken
        nv = self._vt0.n
        G = self._vt0.G
        is_t = G == 0
        is_v = ~is_t

        # TODO: was trying avoid pre-allocating empty array, but for now gonna do it...
        self.hist = self._res_to_xr(np.empty((nt + 1, nv)), np.empty((nt + 1, nv)))

        if self.f_write_out_vortons:
            vortons_file = FORT_BASE_DIR / "out/vortons.csv"
            data = np.genfromtxt(vortons_file, delimiter=",", skip_header=1, skip_footer=sf)
            nrows = data.shape[0]
            i1 = np.arange(0, nrows - 1, 2)
            i2 = np.arange(1, nrows, 2)
            # need to swap dims because t is first in hist
            self.hist["x"].loc[dict(v=is_v)] = data[i1, :].T
            self.hist["y"].loc[dict(v=is_v)] = data[i2, :].T

        if self.f_write_out_tracers:
            tracers_file = FORT_BASE_DIR / "out/tracers.csv"
            data = np.genfromtxt(tracers_file, delimiter=",", skip_header=1, skip_footer=sf)
            nrows = data.shape[0]
            i1 = np.arange(0, nrows - 1, 2)
            i2 = np.arange(1, nrows, 2)
            self.hist["x"].loc[dict(v=is_t)] = data[i1, :].T
            self.hist["y"].loc[dict(v=is_t)] = data[i2, :].T

        # note: the ps code of the Fortran model only works for a specific case
        # (initial equi tri with the second point having x=0 and y>0)
        if self.f_write_out_ps:
            ps_file = FORT_BASE_DIR / "out/ps.txt"
            data = np.genfromtxt(ps_file, skip_header=1, skip_footer=sf)
            self.ps = data


JL_BASE_DIR = Path(__file__).parent / "jl"


class Model_jl(ModelBase):
    """Interface to the Julia model, whose source code is in `vorts/jl/vorts.jl`.

    .. note::
       In this implementation we communicate with Julia using
       [pyjulia](https://github.com/JuliaPy/pyjulia)
    """

    def __init__(
        self,
        vortons: Vortons = None,
        tracers: Tracers = None,
        *,
        dt=0.1,
        nt=1000,
        # above are passed to base
        int_scheme_name="Tsit5",
        use_sysimage=False,
        run_method="pyjulia",
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
            Time integration scheme name. Any from the
            [ODE solvers list](https://diffeq.sciml.ai/stable/solvers/ode_solve/)
            are valid.
        use_sysimage : bool
            Whether to (try to) use the `sys_vorts.so` created by the `compile.jl` script.
        run_method : str, {'pyjulia', 'diffeqpy'}
            Currently, we can run using pyjulia alone or using diffeqpy.
            Both use the tendency function defined in `vorts.jl`.
        """
        # Call base initialization
        super().__init__(vortons, tracers, dt=dt, nt=nt)

        # Other inputs
        self.int_scheme_name = int_scheme_name
        self.use_sysimage = use_sysimage
        self.run_method = run_method

    @staticmethod
    def _load_jl_sysimage():
        # Use our sysimage with DifferentialEquations (and our tend fn) pre-compiled
        # (has to be done before importing any Julia modules)
        from julia import Julia

        jl = Julia(sysimage=str(JL_BASE_DIR / "sys_vorts.so"))
        # TODO: might want to also support `Julia(compiled_modules=False)` in some way
        return jl

    # pyjulia alone
    def _run_pyjulia(self):
        from julia.core import JuliaError

        if self.use_sysimage:
            self._load_jl_sysimage()

        # For some reason it doesn't work the first time on Windows (can't find PyCall)
        # but will work if we just do it again...
        try:
            from julia import Base, Main, Pkg
        except JuliaError as e:
            warnings.warn(f"Julia imports failed the first time with message:\n{e}")
            from julia import Base, Main, Pkg

        # Load our code
        # Base.println("Initial env status:")
        # Pkg.status()
        Pkg.activate(JL_BASE_DIR.as_posix())
        Pkg.status()
        Main.include((JL_BASE_DIR / "vorts.jl").as_posix())
        Base.println("vorts.jl has been loaded!")

        # Run by calling `integrate` from `vorts.jl`
        x, y = Main.integrate(
            self._vt0.state_mat,
            self._vt0.G,
            self.dt,
            self.nt,
            int_scheme_name=self.int_scheme_name,
            ret_for="python",
        )

        # Set output dataset
        self.hist = self._res_to_xr(x, y)

    # pyjulia + diffeqpy
    def _run_diffeqpy(self):
        from julia.core import JuliaError

        if self.use_sysimage:
            self._load_jl_sysimage()

        try:
            from julia import Base, Main, Pkg
        except JuliaError as e:
            warnings.warn(f"Julia imports failed the first time with message:\n{e}")
            from julia import Base, Main, Pkg

        Pkg.activate(JL_BASE_DIR.as_posix())
        Pkg.status()

        # `ode` from the readme doesn't work yet in stable (would have to install from GH)
        # This one would also take two tries in Win if we hadn't already done the above import
        from diffeqpy import de

        # Add the code with the tendency function
        Main.include((JL_BASE_DIR / "vorts.jl").as_posix())
        Base.println("vorts.jl has been loaded!")

        # Set up our inputs
        nt = self.nt
        dt = self.dt
        vt0 = self._vt0
        r0 = vt0.state_mat.T  # note transpose!
        tspan = (0, dt * nt)
        p = Main.TendParams(vt0.G)
        solver = getattr(de, self.int_scheme_name)()

        prob = de.ODEProblem(Main.tend_b, r0, tspan, p)
        sol = de.solve(prob, solver, saveat=np.arange(0, nt + 1) * dt)

        # Exract data and set the output dataset
        x = np.empty((self.nt + 1, vt0.n))
        y = np.empty_like(x)
        assert len(sol.u) == nt + 1
        assert sol.u[0].shape == (2, vt0.n)
        for i in range(len(sol.u)):
            x[i, :] = sol.u[i][0, :]
            y[i, :] = sol.u[i][1, :]

        self.hist = self._res_to_xr(x, y)

    def _run(self):
        if self.run_method == "pyjulia":
            self._run_pyjulia()
        elif self.run_method == "diffeqpy":
            self._run_diffeqpy()
        else:
            raise ValueError("`run_method` must be 'pyjulia' or 'diffeqpy'.")
