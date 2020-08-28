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
        """Create vorton.

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


class model_py():
    """Model in Python."""

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

        # TODO: make module of class variable after moving steppers away
        self.manual_steppers = {
            'FT': FT_step,
            'FT_2': FT_2_step,
            'RK4': RK4_step,
            'RK4_2': RK4_2_step,
            'RK4_3': RK4_3_step,
        }
        self.scipy_methods = {  # could create this programatically
            "scipy_RK45": "RK45",
            "scipy_DOP853": "DOP853",
        }
        self._allowed_int_scheme_names = list(self.manual_steppers.keys()) + list(self.scipy_methods.keys())
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
                stepper=self.manual_steppers[self.int_scheme_name],
                **self.int_scheme_kwargs
            )
        else:  # Scipy integrator
            integrate_scipy(
                self.vortons,
                t_eval,
                self.G_vals,
                method=self.scipy_methods[self.int_scheme_name],
                max_step=self.dt,
                **self.int_scheme_kwargs
            )


def integrate_scipy(
    vortons: List[Vorton],
    t_eval,
    G_vals,
    *,
    # solve_ivp kwargs
    method="RK45",
    max_step: float,
    **options,
):
    """Integrate using `scipy.integrate.solve_ivp`.
    
    **options
        passed through to `solve_ivp`
    """
    from scipy.integrate import solve_ivp

    def calc_lsqd(xarr, yarr):
        # avoid dividing by 0 by adding I
        # taken from RK4_2_step
        return (xarr - xarr.T)**2 + (yarr - yarr.T)**2 + 1e-13*np.eye(*xarr.shape)

    def fun(t, y, G):
        """Calculate both x- and y-tend."""
        # based on calc_tend in RK4_2_step
        # assume columns of `y` are G, x, y

        # TODO: try pulling G out instead (or pass in additional args), since it has no tend
        # G, xpos, ypos = y.T  # rows are unpacked, so need to transpose

        m = G.size  # number of vortons
        xpos, ypos = y[:m], y[m:]

        xarr, yarr = np.meshgrid(xpos, ypos)

        xdiff = xarr - xarr.T
        ydiff = yarr - yarr.T

        lsqd = calc_lsqd(xarr, yarr)

        return np.concatenate((
            -1/(2*np.pi) * np.sum(G.T * ydiff / lsqd, axis=1),  # x-tend
            1/(2*np.pi) * np.sum(G * xdiff / lsqd, axis=0)  # y-tend
        ))

    t_span = 0, t_eval[-1]

    # return of fun(t, y) must have shape (n,) or (n, k)
    # but y0 must be 1-dimensional (n,)
    Gs = np.array(G_vals)[:, np.newaxis]  # column vector (so it can be transposed)

    x0s = np.array([v.xhist[0] for v in vortons])
    y0s = np.array([v.yhist[0] for v in vortons])
    y0 = np.concatenate((x0s, y0s))

    res = solve_ivp(
        fun, 
        t_span=t_span, 
        y0=y0,
        t_eval=t_eval,
        #
        method=method, 
        vectorized=True,  # not sure what impact this has...
        args=(Gs,),  # additional arguments after `t, y` used in fun, jac, etc.
        max_step=max_step,
        **options,
    )

    # update hists
    m = Gs.size  # number of vortons
    for i, v in enumerate(vortons):
        v.xhist[1:] = res.y[i, :]
        v.yhist[1:] = res.y[i+m, :]


def integrate_manual(
    vortons: List[Vorton],
    C_0: float,  # could be computed instead (so could be optional arg)
    t_eval,
    G_vals,
    stepper,
    *,
    adapt_tstep=False,
    use_tqdm=True,
    C_err_reltol: float = 1.0e-9,
    dt_min: float = 1.0e-4,

):
    """Integration routine for use with my handwritten FT/RK4 steppers.

    Optional naive adaptive time-stepping.

    t_eval : array_like
        times to store
        not including 0, which has already been stored

    """
    if adapt_tstep:
        use_tqdm = False  # override this for now

    # initial previous C val is C_0
    C_lm1 = C_0

    # check for non-constant dt
    dt_eval = np.diff(t_eval)
    dt0 = dt_eval[0]
    # if not np.all(dt_eval == dt0):  # fails for some reason even when used np.arange
    if not np.allclose(dt_eval, dt0):
        print(dt_eval, dt0, dt_eval==dt0)
        raise NotImplementedError("non-constant dt in `t_eval`")

    # optionally use tqdm
    iter_l = range(1, len(t_eval)+1)
    if use_tqdm:
        iter_l = tqdm(iter_l)

    # iterate over time index `l`
    for l in iter_l:  # hists at time l=0 are already set in Vorton class init

        Gs = G_vals
        x_lm1s = [v.xhist[l-1] for v in vortons]  # positions at time lm1 (the previous step)
        y_lm1s = [v.yhist[l-1] for v in vortons]
        state00 = list(zip(Gs, x_lm1s, y_lm1s))  # state before taking this step
        # `list(zip(` in order to take length etc.
        if adapt_tstep:
            C_relerr = 100
            dt = dt0
            while C_relerr > C_err_reltol and dt >= dt_min:
                state0 = state00.copy()  # before sub-stepping
                for ll in range(int(np.round(dt0/dt))):
                    # step
                    xnews, ynews = stepper(state0, dt=dt)
                    state0 = list(zip(Gs, xnews, ynews))  # new state0 for potential next step

                # decrease dt for potential next round
                dt /= 2

                # update hists (calc_C currently uses the hists and needs this)
                for i, v in enumerate(vortons):
                    v.xhist[l] = xnews[i]
                    v.yhist[l] = ynews[i]

                # calculate C error
                C_l = calc_C(vortons, l)
                C_relerr = np.abs((C_l - C_lm1) / C_lm1)

            print(f"tstep: {l:d}, min dt used: {dt:.1e}")

        else:  # no adaptive time stepping

            # step
            xnews, ynews = stepper(state00, dt=dt0)

            # store
            for i, v in enumerate(vortons):
                v.xhist[l] = xnews[i]
                v.yhist[l] = ynews[i]

        # calculate new previous C val
        C_lm1 = calc_C(vortons, l)




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
    This doesn't affect the result for FT, but for RK4 it does (since sub-steps are used)...
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
    G = np.asarray(G)[:, np.newaxis]  # column vector (so can transpose)

    # TODO: could pull these out of here and append `_mat` to names?
    def calc_lsqd(xdiff, ydiff):
        """Avoid dividing by lsqd=0 by adding I"""
        return (xdiff)**2 + (ydiff)**2 + 1e-10*np.eye(len(state0))


    def calc_tend(x, y):
        """Calculate both x- and y-tend."""
        # a more-optimized calculation trying to reduce mem usage / repetition
        # but still is slower than RK4_3 (at least for nt=2000, 3 vortons)
        # (seems to overtake RK4_3 in performance as increase N (vortons))

        xarr, yarr = np.meshgrid(x, y)

        xdiff = xarr - xarr.T
        ydiff = yarr - yarr.T

        lsqd = calc_lsqd(xdiff, ydiff)

        return (
            -1/(2*np.pi) * np.sum(G.T * ydiff / lsqd, axis=1),  # x-tend
            1/(2*np.pi) * np.sum(G * xdiff / lsqd, axis=0)  # y-tend
        )

    # TODO: x, y together in matrix?
    # RK4
    x1 = x0
    y1 = y0
    k1x, k1y = calc_tend(x1, y1)

    x2 = x0 + dt/2*k1x
    y2 = y0 + dt/2*k1y
    k2x, k2y = calc_tend(x2, y2)

    x3 = x0 + dt/2*k2x
    y3 = y0 + dt/2*k2y
    k3x, k3y = calc_tend(x3, y3)

    x4 = x0 + dt/1*k3x
    y4 = y0 + dt/1*k3y
    k4x, k4y = calc_tend(x4, y4)

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

    # RK4
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
