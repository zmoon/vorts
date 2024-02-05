"""
Integration routines and time steppers in Python.

In addition to the sub-par handwritten RK and FT schemes,
SciPy RK routines (which are written in Python as well) are also available.
"""

import numba
import numpy as np
from tqdm import tqdm

from ..vortons import Vorton  # noqa: F401

_TEND_PRE = 1 / (2 * np.pi)


def integrate_scipy(
    y0,
    t_eval,
    G_col,
    *,
    # `solve_ivp` kwargs
    method: str = "RK45",
    max_step: float,
    **options,
):
    r"""Integrate using [`scipy.integrate.solve_ivp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html),
    to which `**options` are passed, along with `t_eval`, `method`, and `max_step`.

    Parameters
    ----------
    y0 : array_like
        Vorton state ($x$ and $y$ only) as a single column.
    G_col : array_like
        Array of $\Gamma$ values as a column vector, e.g., from `vorts.vortons.Vortons.G_col`.
    """
    from scipy.integrate import solve_ivp

    def calc_lsqd(xarr, yarr):
        # avoid dividing by 0 by adding I
        # taken from RK4_2_step
        # fmt: off
        return (xarr - xarr.T)**2 + (yarr - yarr.T)**2 + 1e-13*np.eye(*xarr.shape)
        # fmt: on

    def fun(t, y, G):
        """Calculate both x- and y-tend."""
        # based on calc_tend in RK4_2_step

        # unpack
        m = G.size  # number of vortons
        xpos, ypos = y[:m], y[m:]

        xarr, yarr = np.meshgrid(xpos, ypos)

        xdiff = xarr - xarr.T
        ydiff = yarr - yarr.T

        lsqd = calc_lsqd(xarr, yarr)

        return np.concatenate(
            (
                -1 / (2 * np.pi) * np.sum(G.T * ydiff / lsqd, axis=1),  # x-tend
                1 / (2 * np.pi) * np.sum(G * xdiff / lsqd, axis=0),  # y-tend
            )
        )

    t_span = 0, t_eval[-1]

    # return of fun(t, y) must have shape (n,) or (n, k)
    # but y0 must be 1-dimensional (n,)
    res = solve_ivp(
        fun,
        t_span=t_span,
        y0=y0,
        t_eval=t_eval,
        #
        method=method,
        vectorized=method in ("RK45", "DOP853"),  # implicit ones don't work with `vectorized=True`
        args=(G_col,),  # additional arguments after `t, y` used in fun, jac, etc.
        max_step=max_step,
        **options,
    )

    return res.y  # return data only


def integrate_manual(
    G,
    x0,
    y0,
    C_0: float,  # TODO: could be computed instead (so could be optional arg)
    t_eval,
    stepper,  # f(G, x, y)
    *,
    adapt_tstep: bool = False,
    use_tqdm: bool = True,
    C_err_reltol: float = 1.0e-9,
    dt_min: float = 1.0e-4,
):
    r"""Integration routine for use with my handwritten FT/RK4 steppers.

    Optional naive adaptive time-stepping by setting `adapt_tstep=True`.

    Parameters
    ----------
    G, x0, y0 : array_like
        Vectors of $\Gamma$ and initial positions ($x$, $y$) for the system of vortons.
    C_0: float
        Initial of value of $C$ (see description in `vorts.vortons.Vortons.C`)
        to compare to if doing adaptive time stepping.
    t_eval : array_like
        Times to store (not including $t=0$, which has already been stored).
    stepper
        A time stepping function that returns new positions after one step,
        e.g., one from the `vorts.py.integ.MANUAL_STEPPERS` dict.
    adapt_tstep: bool
        Whether to apply adaptive time-stepping.
    use_tqdm: bool, str
        If `True`, or `'console'`, the console version of tqdm is used.
        Pass `'notebook'` to activate the Jupyter widget version of tqdm.
    """
    if adapt_tstep:
        use_tqdm = False  # override this for now

    nt = t_eval.size - 1  # number of integration steps
    nv = x0.size  # number of vortons

    # initial previous C val is C_0
    C_lm1 = C_0

    # check for non-constant dt
    dt_eval = np.diff(t_eval)
    dt0 = dt_eval[0]
    # if not np.all(dt_eval == dt0):  # fails for some reason even when used np.arange
    if not np.allclose(dt_eval, dt0):
        # print(dt_eval, dt0, dt_eval==dt0)
        raise NotImplementedError("non-constant dt in `t_eval`")

    # optionally use tqdm
    iter_l = range(1, nt + 1)  # note starting at 1, not 0
    if use_tqdm is True or use_tqdm == "console":
        iter_l = tqdm(iter_l)
    elif use_tqdm == "notebook":
        from tqdm import tqdm_notebook

        iter_l = tqdm_notebook(iter_l)

    # pre-allocate return arrays
    xhist = np.empty((nv, nt + 1))
    yhist = np.empty_like(xhist)
    xhist[:, 0] = x0
    yhist[:, 0] = y0

    # iterate over time index `l`
    x_lm1 = x0  # l-1 starts at 0
    y_lm1 = y0
    for l in iter_l:  # noqa: E741
        # adaptive time stepping
        if adapt_tstep:
            C_relerr = 100
            dt = dt0
            while C_relerr > C_err_reltol and dt >= dt_min:
                # sub-steps
                x_l = x_lm1.copy()
                y_l = y_lm1.copy()
                for _ in range(int(np.round(dt0 / dt))):
                    # step
                    x_l, y_l = stepper(G, x_l, y_l, dt=dt)

                # decrease dt for potential next round
                dt /= 2

                # calculate C error
                C_l = calc_C(G, x_l, y_l)
                C_relerr = np.abs((C_l - C_lm1) / C_lm1)

            # print(f"tstep: {l:d}, min dt used: {dt:.1e}")

            # new C_lm1 for next main time step
            C_lm1 = C_l

        # no adaptive time stepping
        else:
            # step
            x_l, y_l = stepper(G, x_lm1, y_lm1, dt=dt0)

        # store
        xhist[:, l] = x_l
        yhist[:, l] = y_l

        # new lm1 terms
        x_lm1 = x_l
        y_lm1 = y_l

    return xhist, yhist


# TODO: try numba (especially for the looping ones) (wrapping in njit could be a model option)

# TODO: take advantage of the fact that tracers are only impacted by the vortons
# (calculations involved in their tends can be reduced)

# TODO: streamline calculating x/y tends together (reducing mem usage, etc.)


def calc_lsqd_xy(x1, y1, x2, y2):
    """Calculate intervortical distance $l^2$, using positions to compute `dx` and `dy`."""
    return (x1 - x2) ** 2 + (y1 - y2) ** 2


def calc_lsqd_diff(dx, dy):
    """Calculate intervortical distance $l^2$, passing in already-computed `dx` and `dy`."""
    return dx**2 + dy**2


def calc_C(G, x, y):
    r"""Calculate $C$ at time $l$ (see equation and info in `vorts.vortons.Vortons.C`).

    We use deparature from $C_0$ (initial value of $C$) in the adaptive stepping
    to see if we need to go back and step with smaller $\delta t$.
    """
    nv = x.size  # number of vortons

    C = 0
    for i, j in zip(*np.triu_indices(nv, 1)):
        xi, yi = x[i], y[i]
        xj, yj = x[j], y[j]

        lij_sqd = (xi - xj) ** 2 + (yi - yj) ** 2

        Gi = G[i]
        Gj = G[j]

        C += Gi * Gj * lij_sqd

    return C


def calc_tend_vec(G, x, y):
    """Calculate both $x$- and $y$-tend written in a vectorized way."""
    # a more-optimized calculation trying to reduce mem usage / repetition
    # but still is slower than RK4_3 (at least for nt=2000, 3 vortons)
    # (seems to overtake RK4_3 in performance as increase N (vortons))

    # note: meshgrid not allowed in numba
    xarr, yarr = np.meshgrid(x, y)

    dx = xarr - xarr.T
    dy = yarr - yarr.T

    # avoid dividing by lsqd=0 by adding I
    # fmt: off
    lsqd = calc_lsqd_diff(dx, dy) + 1e-10*np.eye(*xarr.shape)
    # // lsqd = xdiff**2 + ydiff**2

    return (
        -1/(2*np.pi) * np.sum(G.T * dy / lsqd, axis=1),  # x-tend
        1/(2*np.pi) * np.sum(G * dx / lsqd, axis=0)  # y-tend
    )
    # fmt: on


@numba.njit
def calc_tend_vec_premesh(G, X, Y):
    """Calculate both $x$- and $y$-tend vectorized, but they must be passed in in meshgrid form.

    Note that Numba [doesn't support](https://numba.pydata.org/numba-doc/dev/reference/numpysupported.html) `numpy.meshgrid`.

    .. note::
       This function is wrapped with `numba.njit`.
    """
    # with pre-calculated meshgrid so can use Numba
    dx = X - X.T
    dy = Y - Y.T

    # avoid dividing by lsqd=0 by adding I
    # fmt: off
    lsqd = dx**2 + dy**2 + 1e-10*np.eye(*X.shape)

    return (
        -1/(2*np.pi) * np.sum(G.T * dy / lsqd, axis=1),  # x-tend
        1/(2*np.pi) * np.sum(G * dx / lsqd, axis=0)  # y-tend
    )
    # fmt: on


@numba.njit
def calc_tend(G, x, y):
    """Calculate tendencies for each vorton (or tracer).
    Using explicit loops intead of vector(ized) operations.

    .. note::
       This function is wrapped with `numba.njit`.
    """
    # pre-allocate
    dxdt = np.zeros(x.size)
    dydt = np.zeros_like(dxdt)

    # TODO: try filtering out tracers (`np.delete` ?) here to use in the 2nd loop

    # calculate x and y tendencies for each i vorton
    for i, (xi, yi) in enumerate(zip(x, y)):
        # add up contributions
        for j, (Gj, xj, yj) in enumerate(zip(G, x, y)):
            # only other vortons contribute (not tracers)
            if i != j and Gj != 0:
                dx = xi - xj
                dy = yi - yj
                lij_sqd = dx**2 + dy**2
                # add the contributions of j to i's tends
                dxdt[i] += -_TEND_PRE * Gj * dy / lij_sqd
                dydt[i] += _TEND_PRE * Gj * dx / lij_sqd

    return dxdt, dydt


def calc_tend_one(xi, yi, Gn, xn, yn):
    """Calculate tendencies for one position (`xi`, `yi`) based on others (`Gn`, `xn`, `yn`).
    Using explicit loops intead of vector(ized) operations.
    """
    dxdt, dydt = 0, 0
    for Gj, xj, yj in zip(Gn, xn, yn):
        if Gj != 0:  # tracers don't contribute
            dx = xi - xj
            dy = yi - yj
            lij_sqd = dx**2 + dy**2
            # add contribution of j to i's tends
            dxdt += -_TEND_PRE * Gj * dy / lij_sqd
            dydt += _TEND_PRE * Gj * dx / lij_sqd

    return dxdt, dydt


def FT_step_1b1(G, x0, y0, dt: float):
    """Step using 1st-O forward-in-time.

    Calculate tendencies / integrate one vorton by one.
    This doesn't affect the result for FT, but for RK4 it does (since sub-steps are used)...
    """
    # pre-allocate
    nv = x0.size
    xnew = np.zeros(nv)
    ynew = np.zeros_like(xnew)

    iall = np.arange(nv)
    for i in iall:
        # point i
        xi, yi = x0[i], y0[i]
        # others
        iothers = np.delete(iall, i)  # iothers is a new array
        Gn = G[iothers]
        xn = x0[iothers]
        yn = y0[iothers]
        # calc tends
        dxdt, dydt = calc_tend_one(xi, yi, Gn, xn, yn)
        # increment and store
        xnew[i] = xi + dxdt * dt
        ynew[i] = yi + dydt * dt

    return xnew, ynew


def FT_step(G, x0, y0, dt: float, *, tend_fn=calc_tend):
    """Step using 1st-O forward-in-time.

    Calculate tendencies / integrate all vortons at once.
    """
    # calc tends
    dxdt, dydt = tend_fn(G, x0, y0)
    # increment
    xnew = x0 + dxdt * dt
    ynew = y0 + dydt * dt

    return xnew, ynew


def RK4_step_1b1(G, x0_all, y0_all, dt: float):
    """Step using RK4.

    One vorton at a time.

    .. warning::
       This doesn't work for RK4 since we have sub-steps where the tendencies due to
       all other vortons need to be up to date or we introduce error.

       So don't use this if you want to get a good solution.
    """
    # pre-allocate
    nv = x0_all.size
    xnew = np.zeros(nv)
    ynew = np.zeros_like(xnew)

    iall = np.arange(nv)
    for i in iall:
        # point i
        x0, y0 = x0_all[i], y0_all[i]
        # others
        iothers = np.delete(iall, i)  # iothers is a new array
        Gn = G[iothers]
        xn = x0_all[iothers]
        yn = y0_all[iothers]

        # RK4
        # fmt: off
        x1 = x0
        y1 = y0
        k1x, k1y = calc_tend_one(x1, y1, Gn, xn, yn)

        x2 = x0 + dt/2*k1x
        y2 = y0 + dt/2*k1y
        k2x, k2y = calc_tend_one(x2, y2, Gn, xn, yn)

        x3 = x0 + dt/2*k2x
        y3 = y0 + dt/2*k2y
        k3x, k3y = calc_tend_one(x3, y3, Gn, xn, yn)

        x4 = x0 + dt/1*k3x
        y4 = y0 + dt/1*k3y
        k4x, k4y = calc_tend_one(x4, y4, Gn, xn, yn)

        xnew[i] = x0 + dt/6*(k1x + 2*k2x + 2*k3x + k4x)
        ynew[i] = y0 + dt/6*(k1y + 2*k2y + 2*k3y + k4y)
        # fmt: on

    return xnew, ynew


def RK4_step(G, x0, y0, dt: float, *, tend_fn=calc_tend):
    """Step using RK4.

    Whole system at once -- matrix math version.
    """
    # RK4
    x1 = x0
    y1 = y0
    k1x, k1y = tend_fn(G, x1, y1)

    # fmt: off
    x2 = x0 + dt/2*k1x
    y2 = y0 + dt/2*k1y
    k2x, k2y = tend_fn(G, x2, y2)

    x3 = x0 + dt/2*k2x
    y3 = y0 + dt/2*k2y
    k3x, k3y = tend_fn(G, x3, y3)

    x4 = x0 + dt/1*k3x
    y4 = y0 + dt/1*k3y
    k4x, k4y = tend_fn(G, x4, y4)

    xnew = x0 + dt/6*(k1x + 2*k2x + 2*k3x + k4x)
    ynew = y0 + dt/6*(k1y + 2*k2y + 2*k3y + k4y)
    # fmt: on

    return xnew, ynew


# these two dicts are used by model_py to select integration method
MANUAL_STEPPERS = {
    "FT": FT_step,
    "FT_1b1": FT_step_1b1,
    "RK4": RK4_step,
    "RK4_1b1": RK4_step_1b1,
}
"""Time steppers (return new positions after stepping forward in time once)."""

# for selecting tend fn when creating model, like stepper with `int_scheme_name`
# TODO: implement selecting these when creating model
TEND_FNS = {
    "double-loop": calc_tend,
    "vec": calc_tend_vec,
    "vec_premesh": calc_tend_vec_premesh,
}

# TODO: method of registering user custom stepper/tend fns

# keys are used to select integration method when creating model
# values are used as `method` for `scipy.integrate.solve_ivp`'s
SCIPY_METHODS = {
    f"scipy_{method}": method for method in ["RK45", "DOP853", "Radau", "BDF", "LSODA"]
}
