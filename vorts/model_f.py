"""
Python driver for the Fortran version.
"""
import os
from pathlib import Path
import subprocess

import numpy as np

from ._model import ModelBase
from .vortons import Vortons, Tracers


FORT_BASE_DIR = Path(__file__).parent / "f" #.absolute()
# ^ the one that `bin`, `in`, `out`, `src` are in
assert (FORT_BASE_DIR / "src").exists()  # make sure that this is the right spot


def fort_bool(b: bool):
    """Convert Python boolean to a string of the Fortran form."""
    # TODO: shorten
    if b:
        return ".true."
    else:
        return ".false."


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
        if not self.vorts_exe_path.exists():
            raise Exception(
                f"{self.vorts_exe_path!r} doesn't exist. "
                "The Fortran code must first be compiled to produce this file."
            )
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

    # implement abstract method `_run`
    def _run(self):
        """Invoke the Fortran model's executable and load the results."""
        # exe_abs = str(self.vorts_exe_path)
        exe_rel = str(self.vorts_exe_path.relative_to(FORT_BASE_DIR))
        cmd = exe_rel

        # invoke the Fortran model's executable
        cwd = os.getcwd()
        os.chdir(FORT_BASE_DIR)
        os.system('rm ./out/*')
        # print(cmd)
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
