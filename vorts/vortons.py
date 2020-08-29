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


# could exchange x,y for r at some point, to open 3-d option more easily
class Vortons:
    """Collection of Vortons."""
    def __init__(self, G, x, y):
        """Create vorton collection.

        Parameters
        ----------
        G, x, y : array_like (n_vortons,)
            G: Gamma (strength of the circulation, with sign to indicate direction)
                In fluid dynamics, circulation $\Gamma$ is the line integral of velocity
                or flux of vorticity vectors through a surface (here the xy-plane).
            x: x position
            y: y position

        """
        self.G = np.asarray(G)

        # the state matrix has shape (n_vortons, n_pos_dims) (G excluded since time-invariant)
        x = np.asarray(x, dtype=np.float)
        y = np.asarray(y, dtype=np.float)
        self.state_mat = np.column_stack((x, y))

        assert self.G.ndim == 1 and self.state_mat.shape[1] == 2

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
        # TODO: should this call `self._update_vortons`? so as to not be out-of-date if state changes?

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
        """Plot the vortons.
        (Only their current positions, which are all this container knows about.)
        """
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()

        # plot vorton positions
        c_Gp = "cadetblue"
        c_Gm = "salmon"
        G = self.G
        Gp, Gm = G > 0, G < 0
        x, y = self.x, self.y
        ax.plot(x[Gp], y[Gp], "o", ms=7, c=c_Gp, label="$\Gamma > 0$")
        ax.plot(x[Gm], y[Gm], "o", ms=7, c=c_Gm, label="$\Gamma < 0$")

        # plot center-of-mass
        x_cm, y_cm = self.cm()
        s_cm = f"({x_cm:.4g}, {y_cm:.4g})"
        ax.plot(x_cm, y_cm, "*", ms=13, c="gold", label=f"center-of-mass\n{s_cm}")

        # 2nd mom
        x_cm2, y_cm2 = self.mom(2)
        s_cm2 = f"({x_cm2:.4g}, {y_cm2:.4g})"
        ax.plot(x_cm2, y_cm2, "*", ms=13, c="0.4", label=f"mom2\n{s_cm2}")

        # 3nd mom
        # TODO: helper fn to DRY this
        x_cm3, y_cm3 = self.mom(3)
        s_cm3 = f"({x_cm3:.4g}, {y_cm3:.4g})"
        ax.plot(x_cm3, y_cm3, "*", ms=13, c="0.55", label=f"mom3\n{s_cm3}")

        ax.set(
            title=f"$C = {self.C():.4g}$",
            xlabel="$x$",
            ylabel="$y$",
        )
        ax.set_aspect("equal", "box")
        fig.legend()
        ax.grid(True)
        fig.tight_layout()

        # return


    def mom(self, n, *, abs_G=False, center=False):
        """Compute `n`-th moment.

        Parameters
        ----------
        n : int
            which moment
            https://en.wikipedia.org/wiki/Moment_(mathematics)
        abs_G : bool, optional (default False)
            whether to take the absolute value of G values
        center : bool, optional (default True)
            True: evaluate moment wrt. center-of-mass
            False: evaluate moment wrt. (0, 0)
        """
        # seems like a moment but that might not be the correct terminology...
        G = self.G_col
        if abs_G:
            G = np.abs(G)
        G_tot = G.sum()

        x = self.state_mat  # x, y (columns)

        c = self.cm() if center else 0

        x_mom = (G * (x-c)**n).sum(axis=0) / G_tot  # sum along vortons dim, giving a position
        # ^ maybe this should be x - x_cm here...

        return x_mom


    def cm(self):
        """Compute center-of-mass using Gamma as mass."""
        # TODO: what impact should sign of G have on cm? mass is always pos. but G can be neg.
        return self.mom(1, abs_G=True, center=False)


    def center_coords(self, inplace=False):
        """Make center-of-mass (0, 0)."""
        xy_cm = self.cm()
        x_cm, y_cm = xy_cm
        if not inplace:
            return Vortons(self.G, self.x-x_cm, self.y-y_cm)
        else:
            self.state_mat -= x_cm


    @staticmethod
    def regular_polygon(n, *, G=None, **kwargs):
        """Create Vortons with positions corresponding to regular polygon.

        Parameters
        ----------
        n : int
            polygon order
        G : int, array-like, optional
            Gamma value(s) to use
            single value or array of values
            default: 1.0

        `**kwargs` are passed on to `vortons.regular_polygon_vertices`.
        See signature there.
        """
        if G is None:
            G = 1.0  # default
        G = np.asarray(G)
        if G.size == 1:  # single G provided, or using the default
            G = np.full((n,), G)  # TODO: could also the constructor to accept single G
        if G.size != n:
            raise ValueError(f"`G` must have size `n` or 1, but is {G.size!r}")

        xy = regular_polygon_vertices(n, **kwargs).T  # x, y cols-> rows (for unpacking)

        return Vortons(G, *xy)


    # TODO: indexing dunder methods



    # TODO: class method to take List[Vorton] and return a Vortons?



# class Tracers(Vortons):



def rotmat_2d(ang_deg):
    """Return rotation matrix for rotation `ang_deg` in degrees.
    For left-multiplication of a column position vector.

    `scipy.spatial.transform.Rotation` can be used for 3-d.
    """
    ang = np.deg2rad(ang_deg)
    c, s = np.cos(ang), np.sin(ang)
    R = np.array([
        [c, -s],
        [s, c]
    ])
    return R



def regular_polygon_vertices(n, *, c=(0, 0), r_c=1):
    """Regular polygon vertices.

    Parameters
    ----------
    n : int
        order (number of sides/vertices)
    c : 2-tuple / array-like
        center coordinate of the inscribing circle
    r_c : float, int
        radius of the inscribing circle
    """
    c = np.asarray(c)

    # initial vertex
    vert0 = np.r_[0, r_c]

    # rotation matrix -- left-multiplies a column position vector to give rotated position
    rotmat = rotmat_2d(360/n)

    verts = np.full((n, 2), vert0, dtype=np.float)
    # successive rotations
    for i in range(1, n):
        verts[i, :] = (rotmat @ verts[i-1, :][:, np.newaxis]).squeeze()

    return verts + c




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
            "t": ("t", t, {"long_name": "unitless time"}),
        },
        data_vars={
            "x": (("t"), np.empty((n_t,)), {"long_name": "Vorton x position"}),
            "y": (("t"), np.empty((n_t,)), {"long_name": "Vorton y position"}),
        },
        attrs=ds_attrs,
    )

    return ds


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    plt.close("all")

    vs = Vortons([1, 1], [0, 1], [0, 0])
    vs.plot()

    # G sum here is 0, messing up the mom's...
    Vortons([1, -1], [0, 1], [0, 0]).plot()

    Vortons.regular_polygon(3).plot()

    Vortons.regular_polygon(10, c=(1, 0), r_c=0.5).plot()
