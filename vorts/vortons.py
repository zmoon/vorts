"""
Vorton class and Vortons container class (TODO).
"""
import warnings

import numpy as np
import xarray as xr

class Vorton0:
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


# new Vorton
# doesn't know its history
# just current state
from typing import NamedTuple

class Vorton(NamedTuple):
    G: float
    x: float
    y: float



# class Tracer(Vorton):
    # """Tracer -- a vorton with G=0 (no power)."""


class Vortons:
    """Collection of Vortons."""
    def __init__(self, G, x, y):
        """Create vorton collection.

        Parameters
        ----------
        G, x, y : array_like (n_vortons,)
            G: Gamma (strength)
            x: x position
            y: y position

        """
        self.G = np.asarray(G)

        # the state matrix has shape (n_vortons, n_pos_dims) (G excluded since time-invariant)
        self.state_mat = np.column_stack((x, y))

        # create initial corresponding Vorton objects
        self._update_vortons()


    # these 2 don't really need to be property?
    # maybe shouldn't be, to emphasize that state_mat is the real data
    @property
    def state_vec(self):
        """Return flattened state matrix (G not included).

        Needed to feed to `scipy.integrate.solve_ivp`,
        which requires a 1-d array for the `y0` input.
        """
        return self.state_mat.T.flatten()

    @property
    def state_mat_full(self):
        """Return full state matrix: G and positions."""
        return np.column_stack((self.G, self.state_mat))

    @property
    def G_col(self):
        """G as a column vector."""
        return self.G[:, np.newaxis]

    @property
    def x(self):
        # slice indexing should give just a view into `self.state_mat`
        # thus `self.x.base` will return the state mat
        # i.e., `vs.x.base is vs.state_mat`
        return self.state_mat[:,0]

    @property
    def y(self):
        return self.state_mat[:,1]

    def __repr__(self):
        # n_vortons should be too many, so let's show all
        s_vorts = "\n".join(f"  {v}" for v in self._vortons)
        return f"Vortons(\n{s_vorts}\n)"

    def _update_vortons(self):
        self._vortons = [Vorton(G, x, y) for G, x, y in self.state_mat_full]


    def C(self):
        """Calculate $C$.

        $C$ is supposed to be a conserved quantity in this system.
        - Chamecki (2005) eq. 15
        """
        n_vortons = self.G.size
        G = self.G
        C = 0
        for i, j in zip(*np.triu_indices(n_vortons, 1)):  # all combinations without repetition

            xi, yi = self.x[i], self.y[i]
            xj, yj = self.x[j], self.y[j]

            lij_sqd = (xi-xj)**2 + (yi-yj)**2

            Gi, Gj = G[i], G[j]

            C += Gi * Gj * lij_sqd

        return C


    def plot(self):
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()

        # plot vorton positions
        l, = ax.plot(self.x, self.y, "o", ms=7)

        # plot center-of-mass
        x_cm, y_cm = self.cm()
        s_cm = f"({x_cm:.4g}, {y_cm:.4g})"
        ax.plot(x_cm, y_cm, "*", ms=13, c="gold", label=f"center-of-mass\n{s_cm}")

        ax.set(
            title=f"$C = {self.C():.4g}$",
            xlabel="$x$",
            ylabel="$y$",
        )
        ax.set_aspect("equal", "box")
        ax.legend()
        ax.grid(True)
        fig.tight_layout()

        return l


    def cm(self):
        """Compute center-of-mass using Gamma as mass."""
        G = self.G_col
        G_tot = G.sum()
        x = self.state_mat  # x, y (columns)

        # vector center-of-mass
        x_cm = (G * x).sum(axis=0) / G_tot  # sum along vortons dim

        return x_cm



    # TODO: indexing dunder methods


# class Tracers(Vortons):



# should probably be in model_py
# except that model_f could also use it when collecting the data
# until the Fortran model writes nc files
def init_hist(
    n_vorton: int,
    n_time: int,  # in addition to t=0
    dt: float,
    #
    n_tracer=None,
    *,
    ds_attrs=None,
):
    """Create initial history `xr.Dataset`."""

    # if n_tracer is not None:
    if ds_attrs is None:
        ds_attrs = {}

    t = np.arange(0, n_time+1)*dt
    n_t = t.size

    ds = xr.Dataset(
        coords={
            "t": ("t", t, {"long_name": "non-dimensional? time"}),
        },
        data_vars={
            "x": (("t"), np.empty((n_t,)), {"long_name": "Vorton x position"}),
            "y": (("t"), np.empty((n_t,)), {"long_name": "Vorton y position"}),
        },
        attrs=ds_attrs,
    )

    return ds


if __name__ == "__main__":

    vs = Vortons([1, 1, 1], [-0.666, 0, 0.666], [-0.4, 0.8, -0.4])

    vs.plot()
