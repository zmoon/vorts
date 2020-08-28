"""
Python driver for the Fortran version.
"""
import os
import subprocess

import numpy as np

from .vortons import Vorton


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
        int_scheme_name='RK4'
    ):
        """ """

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
        self.nt = nt
        self.int_scheme_name = int_scheme_name  # {'FT', 'RK4'}

        # executing the model
        self.vorts_exe_path = './bin/vorts.exe'
        self.oe = ''  # standard output and error

        # output data
        self.vortons = []
        self.tracers = []

        self.create_inputs()


    def create_inputs(self):
        """
        Create input files for the Fotran model
          describing the initial condition
          and the simulation settings.
        """

        mat = np.vstack((self.G_vals, self.xi_vals, self.yi_vals)).T
        np.savetxt('./in/vorts_in.txt', mat,
                   delimiter=' ', fmt='%.16f', header='Gamma xi yi')

        mat = np.vstack((self.xit_vals.flat, self.yit_vals.flat)).T
        np.savetxt('./in/tracers_in.txt', mat,
                   delimiter=' ', fmt='%.3f', header='xi, yi')

        mat = [self.dt, self.nt, self.int_scheme_name]
        np.savetxt('./in/vorts_sim_in.txt', mat,
                   delimiter=' ', fmt='%s')


    def run(self):
        """Invoke the Fortran model's executable."""

        os.system('rm ./out/*')

        #os.system(vorts_exe)
        self.oe = subprocess.check_output(self.vorts_exe_path, stderr=subprocess.STDOUT)

        #self.load_results()


    def load_results(self):
        """Load results from a run of the Fortran model."""

        sf = 0  # there may be blank line at end of the files?

        #> version for new output files follows

        vortons_file = './out/vortons.csv'
        data = np.genfromtxt(vortons_file, delimiter=',', skip_header=1, skip_footer=sf)

        self.vortons = []
        for i in range(0, data.shape[0]-1, 2):
            self.vortons.append(Vorton(self.G_vals[i/2], data[i,:], data[i+1,:], self.nt))

        tracers_file = './out/tracers.csv'
        data = np.genfromtxt(tracers_file, delimiter=',', skip_header=1, skip_footer=sf)

        self.tracers = []
        for i in range(0, data.shape[0]-1, 2):
            self.tracers.append(Vorton(0, data[i,:], data[i+1,:], self.nt))