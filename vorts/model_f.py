"""
Python driver for the Fortran version.
"""
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

from .vortons import Vorton


FORT_BASE_DIR = Path(__file__).parent / "f" #.absolute()
# ^ the one that `bin`, `in`, `out`, `src` are in
assert (FORT_BASE_DIR / "src").exists()  # make sure that this is the right spot


def fort_bool(b: bool):
    if b:
        return ".true."
    else:
        return ".false."


class model_f:
    """Thin wrapper for functionality of the Fortran model in `src/`.

    In this implementation we communicate with the Fortran program via text files.
    """

    def __init__(
        self,
        G,
        xi,
        yi,
        xit=None,
        yit=None,
        dt=0.1,
        nt=1000,
        int_scheme_name='RK4',
        #
        write_vortons=True,
        write_tracers=False,
        write_ps=False,
    ):
        """Initialize and create inputs for the Fortran model."""
        # vorton IC arrays
        self.G_vals = np.asarray(G)
        self.xi_vals = np.asarray(xi)
        self.yi_vals = np.asarray(yi)

        # tracer initial positions
        if xit is None or yit is None:
            xit, yit = [], []
        self.xit_vals = np.asarray(xit)
        self.yit_vals = np.asarray(yit)

        # sim settings
        self.dt = dt
        self.nt = int(nt)
        self.int_scheme_name = int_scheme_name  # {'FT', 'RK4'}

        # executing the model
        self.vorts_exe_path = FORT_BASE_DIR / 'bin/vorts.exe'
        if not self.vorts_exe_path.exists():
            raise Exception(
                f"{self.vorts_exe_path!r} doesn't exist. "
                "The Fortran code must first be compiled to produce this file."
            )
        self.oe = ''  # standard output and error

        # output data
        self.vortons = []
        self.tracers = []
        self.write_vortons = write_vortons
        self.write_tracers = write_tracers
        self.write_ps = write_ps

        self.create_inputs()


    def create_inputs(self):
        """
        Create input files for the Fotran model
          describing the initial condition
          and the simulation settings.
        """

        mat = np.vstack((self.G_vals, self.xi_vals, self.yi_vals)).T
        np.savetxt(FORT_BASE_DIR / 'in/vorts_in.txt', mat,
                   delimiter=' ', fmt='%.16f', header='Gamma xi yi')

        mat = np.vstack((self.xit_vals.flat, self.yit_vals.flat)).T
        np.savetxt(FORT_BASE_DIR / 'in/tracers_in.txt', mat,
                   delimiter=' ', fmt='%.3f', header='xi, yi')

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

        self.load_results()


    def load_results(self):
        """Load results from a run of the Fortran model."""

        sf = 0  # there may be blank line at end of the files?

        #> version for new output files follows

        if self.write_vortons:
            vortons_file = FORT_BASE_DIR / 'out/vortons.csv'
            data = np.genfromtxt(vortons_file, delimiter=',', skip_header=1, skip_footer=sf)
            self.vortons = []
            for i in range(0, data.shape[0]-1, 2):
                ihalf = int(i/2)
                self.vortons.append(Vorton(self.G_vals[ihalf], data[i,:], data[i+1,:], self.nt))

        if self.write_tracers:
            tracers_file = FORT_BASE_DIR / 'out/tracers.csv'
            data = np.genfromtxt(tracers_file, delimiter=',', skip_header=1, skip_footer=sf)
            self.tracers = []
            for i in range(0, data.shape[0]-1, 2):
                self.tracers.append(Vorton(0, data[i,:], data[i+1,:], self.nt))

        if self.write_ps:
            ps_file = FORT_BASE_DIR / "out/ps.txt"
            data = np.genfromtxt(ps_file, skip_header=1, skip_footer=sf)
            self.ps = data
