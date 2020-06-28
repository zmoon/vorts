"""
Run Python and Fortran versions of the N-vortex model
"""

from glob import glob
import os
import subprocess
from typing import List, Tuple
import warnings

#import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


class Vorton():
    def __init__(self, G, x, y, nt):
        """ 
        x, y inputs can be either a single point (IC)
        or an array 
        so that same class can be used for catching/organizing the Fortran model output
        """
#        self.x = xi
#        self.y = yi
        self.G = G

        try: 
            self.xhist = np.zeros(nt+1)
            self.yhist = np.zeros(nt+1)
            self.xhist[0] = x
            self.yhist[0] = y
            
        except ValueError:  # ValueError: setting array with seq; TypeError: ints and floats have no len()
            self.xhist = x
            self.yhist = y
            
            if self.xhist.size != nt+1:
                warnings.warn(
                    f"Fortran model output should be size {nt+1:d} but is {self.xhist.size:d}"
                )



class model_f():
    """Thin wrapper for functionality of the Fortran model in `src/`.
    
    In this implementation we communicate with the Fortran program via text files.
    """
    
    def __init__(self, G, xi, yi, 
                 xit=[], yit=[], 
                 dt=0.1, nt=1000, int_scheme_name='RK4'):
        """ """
        
        # vorton IC arrays
        self.G_vals = np.array(G)
        self.xi_vals = np.array(xi)
        self.yi_vals = np.array(yi)
        
        # tracer initial positions
        self.xit_vals = np.array(xit)
        self.yit_vals = np.array(yit)
        
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
        "Invoke the Fortran model's executable."
        
        os.system('rm ./out/*')
        
        #os.system(vorts_exe)
        self.oe = subprocess.check_output(self.vorts_exe_path, stderr=subprocess.STDOUT)
        
        #self.load_results()
        
    
    def load_results(self):
        "Load results from a run of the Fortran model."

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



class model_py():
    """Model in Python."""
    
    def __init__(self, G, xi, yi, 
                 dt=0.1, nt=1000, int_scheme_name='RK4_3', 
                 adapt_tstep=False):
        """Create model with given settings."""

        # vorton IC arrays
        self.G_vals = G
        self.xi_vals = xi
        self.yi_vals = yi
        
        # sim settings
        self.dt = float(dt)
        self.nt = nt
        self.int_scheme_name = int_scheme_name  # {'FT', 'RK4'}
        self.adapt_tstep = adapt_tstep
        
        # set int scheme method from the options
        int_schemes = {
            'FT': FT_step, 
            'FT_2': FT_2_step,
            'RK4': RK4_step, 
            'RK4_2': RK4_2_step, 
            'RK4_3': RK4_3_step,
        }
        self.int_scheme = int_schemes[int_scheme_name]
        
        # create vortons
        self.vortons = []
        for Gxy in zip(self.G_vals, self.xi_vals, self.yi_vals):
            self.vortons.append(Vorton(*Gxy, nt))
            
        self.l = 0  # time step index
        
        # for adaptive time stepping calculations
        self.C_0 = calc_C(self.vortons, self.l)
        self.C_err_reltol = 1e-9
        self.dt_min = 1e-4
            
            
    def run(self):
        """ """
        
        C_lm1 = self.C_0
        
        for self.l in tqdm(range(1, self.nt+1)):  # hists at time l=0 are already set in Vorton class init
            
            Gs = self.G_vals
            x_lm1s = [v.xhist[self.l-1] for v in self.vortons]  # positions at time lm1 (the previous step)
            y_lm1s = [v.yhist[self.l-1] for v in self.vortons]
            
            if self.adapt_tstep:
                
                C_relerr = 100
                dt = self.dt
                while C_relerr > self.C_err_reltol and dt >= self.dt_min:
                    
                    # print(f"tstep: {self.l:d}, dt: {dt:.1e}")
                    
                    state0 = list(zip(Gs, x_lm1s, y_lm1s))  # positions before taking this step
                    # `list(zip(` in order to take length etc.

                    for ll in range(int(np.round(self.dt/dt))):
                    
                        xnews, ynews = self.int_scheme(state0, dt=dt)
                        
                        state0 = list(zip(Gs, xnews, ynews))  # new state0
                        
                        
                    # decrease dt for next round (if necessary)
                    dt /= 2
                    
                    # update hists
                    for i, v in enumerate(self.vortons):
                        v.xhist[self.l] = xnews[i]
                        v.yhist[self.l] = ynews[i]
                        
                    # calculate C error
                    C_l = calc_C(self.vortons, self.l)
#                    C_relerr = np.abs((C_l - self.C_0) / self.C_0)  # may want to compare to C at lm1 instead
                    C_relerr = np.abs((C_l - C_lm1) / C_lm1)
                
                
                print(f"tstep: {self.l:d}, min dt used: {dt:.1e}")

            
            else:
                
                state0 = list(zip(Gs, x_lm1s, y_lm1s))  # positions before taking this step
                
                xnews, ynews = self.int_scheme(state0, dt=self.dt)
                
                for i, v in enumerate(self.vortons):
                    v.xhist[self.l] = xnews[i]
                    v.yhist[self.l] = ynews[i]
                    
            C_lm1 = calc_C(self.vortons, self.l)
  
            
            
def calc_lsqd(x1, y1, x2, y2):
    "Calculate intervortical distance $l^2$."
    
    return (x1-x2)**2 + (y1-y2)**2


def calc_C(vortons: List[Vorton], l: int):
    """Calculate $C$ at time $l$.

    $C$ is supposed to be a conserved quantity in this system 
    - Chamecki (2005) eq. 15
    
    We use deparature from $C_0$ (initial value of $C$) in the adaptive stepping
    to see if we need to go back and step with smaller dt.
    """
    
    C = 0
    for i, j in zip(*np.triu_indices(len(vortons), 1)):
        
        xi, yi = vortons[i].xhist[l], vortons[i].yhist[l]
        xj, yj = vortons[j].xhist[l], vortons[j].yhist[l]
        
        lij_sqd = (xi-xj)**2 + (yi-yj)**2
        
        Gi = vortons[i].G
        Gj = vortons[j].G
        
        C += Gi * Gj * lij_sqd
        
    return C


def calc_xtend(xa: float, ya: float, others: List[Tuple]):
    """Calculate x-tend for one vorton."""
    # note: could use matrix math instead of for-loop

    dxdt = 0
    for Gb, xb, yb in others:
                
        lab_sqd = calc_lsqd(xa, ya, xb, yb)
        
        if lab_sqd:
            dxdt += -1/(2*np.pi) * Gb * (ya-yb) / lab_sqd
        
    return dxdt
        
        
def calc_ytend(xa: float, ya: float, others: List[Tuple]):
    """Calculate y-tend for one vorton."""
    
    dydt = 0
    for Gb, xb, yb in others:
    
        lab_sqd = calc_lsqd(xa, ya, xb, yb)
        
        if lab_sqd:
            dydt += 1/(2*np.pi) * Gb * (xa-xb) / lab_sqd
        
    return dydt


def calc_xtend_all(G0: List[float], x0: List[float], y0: List[float]):
    """Calculate x-tend for all vortons."""
    
    dxdt = np.zeros(x0.size)  # pre-allocate
    G = G0
    
    for i, y0i in enumerate(y0):
        
        x0i = x0[i]
        
        dxdti = 0
        for j, y0j in enumerate(y0):

            x0j = x0[j]
            Gj = G[j]
            
            if i == j:
                pass

            else:
                lij = calc_lsqd(x0i, y0i, x0j, y0j)
                dxdti += -1/(2*np.pi) * Gj * (y0i-y0j) / lij
                
        dxdt[i] = dxdti
        
    return dxdt
        
        
def calc_ytend_all(G0: List[float], x0: List[float], y0: List[float]):
    """Calculate y-tend for all vortons."""
    
    dydt = np.zeros(x0.size)
    G = G0
    
    for i, y0i in enumerate(y0):
        
        x0i = x0[i]
        
        dydti = 0
        for j, y0j in enumerate(y0):

            x0j = x0[j]
            Gj = G[j]
            
            if i == j:
                pass

            else:
                lij = calc_lsqd(x0i, y0i, x0j, y0j)
                dydti += 1/(2*np.pi) * Gj * (x0i-x0j) / lij
                
        dydt[i] = dydti
        
    return dydt


def FT_step(state0: List[Tuple], dt: float):
    """Step using 1st-O forward-in-time. 
    
    Calculate tendencies / integrate one vorton by one. 
    This doesn't affect the result for FT, but for RK4 it does...
    """
    # pre-allocate
    xnews = np.zeros(len(state0))
    ynews = np.zeros_like(xnews)
    
    for i in range(len(state0)):

        _, xa, ya = state0[i]
        
        others = state0[:i] + state0[i+1:] 
    
        dxdt = calc_xtend(xa, ya, others)
        dydt = calc_ytend(xa, ya, others)
        
        xnew = xa + dxdt*dt
        ynew = ya + dydt*dt
        
        xnews[i] = xnew
        ynews[i] = ynew
        
    return xnews, ynews


def FT_2_step(state0: List[Tuple], dt: float):
    """Step using 1st-O forward-in-time. 
    
    Calculate tendencies / integrate all vortons at once.
    """
    
    G, x0, y0 = zip(*state0)

    G = np.array(G)
    x0 = np.array(x0)
    y0 = np.array(y0)
    
    dxdt = calc_xtend_all(G, x0, y0)
    dydt = calc_ytend_all(G, x0, y0)
    
    xnews = x0 + dxdt*dt
    ynews = y0 + dydt*dt
        
    return xnews, ynews
        
        
def RK4_step(state0: List[Tuple], dt: float):
    """Step using RK4. 

    One vorton at a time. 
    **This doesn't work for RK4 since we have sub-steps where the tendencies due to 
    all other vortons need to be up to date or we introduce error.**
    """
    
    xnews = np.zeros(len(state0))
    ynews = np.zeros_like(xnews)
    
    for i in range(len(state0)):

        _, xa, ya = state0[i]
        
        others = state0[:i] + state0[i+1:]
        
        x0 = xa
        y0 = ya
    
        x1 = x0
        y1 = y0
        k1y = calc_ytend(x1, y1, others)
        k1x = calc_xtend(x1, y1, others)
        
        x2 = x0 + dt/2*k1x
        y2 = y0 + dt/2*k1y
        k2x = calc_xtend(x2, y2, others)
        k2y = calc_ytend(x2, y2, others)

        x3 = x0 + dt/2*k2x
        y3 = y0 + dt/2*k2y
        k3x = calc_xtend(x3, y3, others)
        k3y = calc_ytend(x3, y3, others)
        
        x4 = x0 + dt/1*k3x
        y4 = y0 + dt/1*k3y
        k4x = calc_xtend(x4, y4, others)
        k4y = calc_ytend(x4, y4, others)
        
        xnew = x0 + dt/6*(k1x + 2*k2x + 2*k3x + k4x) 
        ynew = y0 + dt/6*(k1y + 2*k2y + 2*k3y + k4y)
        
        xnews[i] = xnew
        ynews[i] = ynew
        
    return xnews, ynews


def RK4_2_step(state0: List[Tuple], dt: float):
    """Step using RK4.

    Whole system at once -- matrix math version.
    """
    
    G, x0, y0 = zip(*state0)
    G = np.array(G)
    G.shape = (G.size, 1)  # column vector
    
    # TODO: could pull these out of here and append `_mat` to names?
    def calc_lsqd(xarr, yarr):
        """ avoid dividing by 0 by adding I """
        
        return (xarr - xarr.T)**2 + (yarr - yarr.T)**2 + 1e-13*np.eye(len(state0))
    
    
    def calc_xtend(x, y):
        
        xarr, yarr = np.meshgrid(x, y)
        
        ydiff = yarr - yarr.T
        
        lsqd = calc_lsqd(xarr, yarr)
        
        return -1/(2*np.pi) * np.sum(G.T * ydiff / lsqd, axis=1)

    
    def calc_ytend(x, y):
        
        xarr, yarr = np.meshgrid(x, y)

        xdiff = xarr - xarr.T
        
        lsqd = calc_lsqd(xarr, yarr)
        
        return 1/(2*np.pi) * np.sum(G * xdiff / lsqd, axis=0)  


    def calc_tend(x, y):
        """Calculate both x- and y-tend."""
        # a more-optimized calculation trying to reduce mem usage / repetition
        # but still is slower than RK4_3 (at least for nt=2000, 3 vortons)
        # (seems to overtake RK4_3 in performance as increase N (vortons))
        
        xarr, yarr = np.meshgrid(x, y)

        xdiff = xarr - xarr.T
        ydiff = yarr - yarr.T
        
        lsqd = calc_lsqd(xarr, yarr)
        
        return (
            -1/(2*np.pi) * np.sum(G.T * ydiff / lsqd, axis=1),  # x-tend
            1/(2*np.pi) * np.sum(G * xdiff / lsqd, axis=0)  # y-tend
        )


    # do that RK4
        
    x1 = x0
    y1 = y0
    # k1x = calc_xtend(x1, y1)
    # k1y = calc_ytend(x1, y1)
    k1x, k1y = calc_tend(x1, y1)
                    
    x2 = x0 + dt/2*k1x
    y2 = y0 + dt/2*k1y
    # k2x = calc_xtend(x2, y2)
    # k2y = calc_ytend(x2, y2)
    k2x, k2y = calc_tend(x2, y2)

    x3 = x0 + dt/2*k2x
    y3 = y0 + dt/2*k2y
    # k3x = calc_xtend(x3, y3)
    # k3y = calc_ytend(x3, y3)
    k3x, k3y = calc_tend(x3, y3)
    
    x4 = x0 + dt/1*k3x
    y4 = y0 + dt/1*k3y
    # k4x = calc_xtend(x4, y4)
    # k4y = calc_ytend(x4, y4)
    k4x, k4y = calc_tend(x3, y3)
    
    xnews = x0 + dt/6*(k1x + 2*k2x + 2*k3x + k4x) 
    ynews = y0 + dt/6*(k1y + 2*k2y + 2*k3y + k4y)
        
    return xnews, ynews


def RK4_3_step(state0: List[Tuple], dt: float):
    """Step using RK4.

    Whole system at once -- version without using matrix math (cf. RK4_2)
    """

    G, x0, y0 = zip(*state0)

    x0 = np.array(x0)
    y0 = np.array(y0)
            
    
    # do that RK4
        
    x1 = x0
    y1 = y0
    k1x = calc_xtend_all(G, x1, y1)
    k1y = calc_ytend_all(G, x1, y1)
                    
    x2 = x0 + dt/2*k1x
    y2 = y0 + dt/2*k1y
    k2x = calc_xtend_all(G, x2, y2)
    k2y = calc_ytend_all(G, x2, y2)

    x3 = x0 + dt/2*k2x
    y3 = y0 + dt/2*k2y
    k3x = calc_xtend_all(G, x3, y3)
    k3y = calc_ytend_all(G, x3, y3)
    
    x4 = x0 + dt/1*k3x
    y4 = y0 + dt/1*k3y
    k4x = calc_xtend_all(G, x4, y4)
    k4y = calc_ytend_all(G, x4, y4)
    
    xnews = x0 + dt/6*(k1x + 2*k2x + 2*k3x + k4x) 
    ynews = y0 + dt/6*(k1y + 2*k2y + 2*k3y + k4y)
        
    return xnews, ynews        

