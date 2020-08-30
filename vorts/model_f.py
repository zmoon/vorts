"""
Python driver for the Fortran version.
"""
import copy
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

from .plot import plot_vorton_trajectories, plot_tracer_trajectories
from .model_py import init_hist
from .vortons import Vortons


FORT_BASE_DIR = Path(__file__).parent / "f" #.absolute()
# ^ the one that `bin`, `in`, `out`, `src` are in
assert (FORT_BASE_DIR / "src").exists()  # make sure that this is the right spot


def fort_bool(b: bool):
    if b:
        return ".true."
    else:
        return ".false."


class Model_f:
    """Thin wrapper for functionality of the Fortran model in `src/`.

    In this implementation we communicate with the Fortran program via text files.
    """
    # _allowed_int_scheme_names = ("FT", "RK4")

    def __init__(
        self,
        vortons: Vortons = None,
        *,
        dt=0.1,
        nt=1000,
        int_scheme_name='RK4',
        #
        write_vortons=True,  # maybe should put `flag` or something in these names
        write_tracers=False,
        write_ps=False,
    ):
        """Initialize and create inputs for the Fortran model."""
        # vortons
        # ! this part should be kept the same as in Model_py
        if vortons is None:
            vortons = Vortons.regular_polygon(3, G=1)
        self.vortons = vortons
        self.vortons0 = copy.deepcopy(self.vortons)

        # tracer initial positions
        # if xit is None or yit is None:
        #     xit, yit = [], []
        # self.xit_vals = np.asarray(xit)
        # self.yit_vals = np.asarray(yit)

        # sim settings
        self.dt = dt
        self.nt = int(nt)
        self.nv = self.vortons0.n
        self.int_scheme_name = int_scheme_name  # {'FT', 'RK4'}

        # executing the model
        self.vorts_exe_path = FORT_BASE_DIR / 'bin/vorts.exe'
        if not self.vorts_exe_path.exists():
            raise Exception(
                f"{self.vorts_exe_path!r} doesn't exist. "
                "The Fortran code must first be compiled to produce this file."
            )
        self.oe = ''  # standard output and error

        # output option flags
        self.write_vortons = write_vortons
        self.write_tracers = write_tracers
        self.write_ps = write_ps

        # initialize hist
        # ! should be same as in Model_py
        self.hist = init_hist(self.nv, self.nt, self.dt)
        self.hist["G"].loc[:] = self.vortons0.G
        t_hist = self.hist.t
        self.hist["x"].loc[dict(t=t_hist[t_hist == 0])] = self.vortons0.x
        self.hist["y"].loc[dict(t=t_hist[t_hist == 0])] = self.vortons0.y

        # write the text input files to directory `vorts/f/in`
        self.create_inputs()

        self._has_run = False


    def create_inputs(self):
        """
        Create input files for the Fotran model
          describing the initial condition
          and the simulation settings.
        """
        # write vorton system initial state
        mat = self.vortons0.state_mat_full()  # needs to be rows of G, xi, yi
        np.savetxt(FORT_BASE_DIR / 'in/vorts_in.txt', mat,
                   delimiter=' ', fmt='%.16f', header='Gamma xi yi')

        # write tracers initial state (positions only)
        # mat = np.vstack((self.xit_vals.flat, self.yit_vals.flat)).T
        # np.savetxt(FORT_BASE_DIR / 'in/tracers_in.txt', mat,
        #            delimiter=' ', fmt='%.3f', header='xi, yi')

        # write model options
        mat = [
            self.dt,
            self.nt,
            self.int_scheme_name,
            fort_bool(self.write_vortons),
            fort_bool(self.write_tracers),
            fort_bool(self.write_ps),
        ]
        np.savetxt(FORT_BASE_DIR / 'in/vorts_sim_in.txt', mat,
                   delimiter=' ', fmt='%s')


    def run(self):
        """Invoke the Fortran model's executable and load the results."""
        # exe_abs = str(self.vorts_exe_path)
        exe_rel = str(self.vorts_exe_path.relative_to(FORT_BASE_DIR))
        cmd = exe_rel

        # try:
        cwd = os.getcwd()
        os.chdir(FORT_BASE_DIR)
        os.system('rm ./out/*')
        # print(cmd)
        self.oe = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        os.chdir(cwd)
        # ^ hack for now, but could instead pass FORT_BASE_DIR into the Fortran program
        #   using vorts_sim_in.txt

        self._has_run = True

        # load results from the text file outputs into `self.hist` (an `xr.Dataset`)
        self.load_results()


    def load_results(self):
        """Load results from a run of the Fortran model."""

        sf = 0  # there may be blank line at end of the files?

        #> version for new output files follows
        t_hist = self.hist.t

        if self.write_vortons:
            vortons_file = FORT_BASE_DIR / 'out/vortons.csv'
            data = np.genfromtxt(vortons_file, delimiter=',', skip_header=1, skip_footer=sf)
            nrows = data.shape[0]
            i1 = np.arange(0, nrows-1, 2)
            i2 = np.arange(1, nrows, 2)
            # self.hist["x"].loc[dict(t=t_hist[t_hist > 0])] = data[i1, :].T  # need to swap dims because t is first in hist
            # self.hist["y"].loc[dict(t=t_hist[t_hist > 0])] = data[i2, :].T
            self.hist["x"].loc[:] = data[i1, :].T  # need to swap dims because t is first in hist
            self.hist["y"].loc[:] = data[i2, :].T

        if self.write_tracers:
            tracers_file = FORT_BASE_DIR / 'out/tracers.csv'
            data = np.genfromtxt(tracers_file, delimiter=',', skip_header=1, skip_footer=sf)
            nrows = data.shape[0]
            i1 = np.arange(0, nrows-1, 2)
            i2 = np.arange(1, nrows, 2)

        if self.write_ps:
            ps_file = FORT_BASE_DIR / "out/ps.txt"
            data = np.genfromtxt(ps_file, skip_header=1, skip_footer=sf)
            self.ps = data


    def plot(self, which="vortons", **kwargs):
        """Plot.

        **kwargs passed through to the corresponding plotting function.
        """
        # ! should be the same as in Model_py
        if not self._has_run:
            raise Exception("The model has not yet been run.")

        if which == "vortons":
            plot_vorton_trajectories(self.hist, **kwargs)

        elif which == "tracers":
            plot_tracer_trajectories(self.hist, **kwargs)

        else:
            raise NotImplementedError(f"which={which!r}")
