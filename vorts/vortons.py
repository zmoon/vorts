"""
Vorton class and Vortons container class (TODO).
"""
import warnings

import numpy as np
import xarray as xr

# class Vorton0:
#     def __init__(self, G, x, y, nt):
#         """Create vorton.

#         x, y inputs can be either a single point (IC)
#         or an array
#         so that same class can be used for catching/organizing the Fortran model output
#         """
# #        self.x = xi
# #        self.y = yi
#         self.G = G

#         try:
#             self.xhist = np.zeros(nt+1)
#             self.yhist = np.zeros(nt+1)
#             self.xhist[0] = x
#             self.yhist[0] = y

#         except ValueError:  # ValueError: setting array with seq; TypeError: ints and floats have no len()
#             self.xhist = x
#             self.yhist = y

#             if self.xhist.size != nt+1:
#                 warnings.warn(
#                     f"Fortran model output should be size {nt+1:d} but is {self.xhist.size:d}"
#                 )


# new Vorton
# doesn't know its history
# just current state
from typing import NamedTuple

class Vorton(NamedTuple):
    G: float
    x: float
    y: float


# TODO: PointVortices ABC that implements adding, has position state_mat, etc.
#       Vortons and Tracers could both be based on it

# not sure if this is necessary...
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
        # TODO: should we check G nonzero? since tracers are treated separately

        # the state matrix has shape (n_vortons, n_pos_dims) (G excluded since time-invariant)
        x = np.asarray(x, dtype=np.float)
        y = np.asarray(y, dtype=np.float)
        self.state_mat = np.column_stack((x, y))

        assert self.G.ndim == 1 and self.state_mat.shape[1] == 2
        assert self.G.size == self.state_mat.shape[0]  # n_vortons

        # create initial corresponding Vorton objects
        self._update_vortons()



    # these 2 don't really need to be property?
    # maybe shouldn't be, to emphasize that state_mat is the real data
    # @property
    def state_vec(self):
        """Return flattened state matrix (G not included).

        Needed to feed to `scipy.integrate.solve_ivp`,
        which requires a 1-d array for the `y0` input.
        """
        return self.state_mat.T.flatten()

    # @property
    def state_mat_full(self):
        """Return full state matrix: G and positions."""
        return np.column_stack((self.G, self.state_mat))

    # seems to return a view into self.G, so ok to be property
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

    @property
    def n(self):
        """Number of vortons."""
        # numpy.ndarray size lookups are esentially free
        return self.G.size  # will have to change if want to allow single G at some point

    def __repr__(self):
        # n_vortons should be too many, so let's show all
        s_vorts = "\n".join(f"  {v}" for v in self._vortons)
        return f"Vortons(\n{s_vorts}\n)"
        # TODO: should this call `self._update_vortons`? so as to not be out-of-date if state changes?

    def _update_vortons(self):
        self._vortons = [Vorton(G, x, y) for G, x, y in self.state_mat_full()]


    def C(self):
        """Calculate $C$.

        $$
        C = \sum_{\alpha, \beta = 1; \alpha \neq \beta}^{N}
            \Gamma_{\alpha} \Gamma_{\beta} l_{\alpha \beta}^{2}
        $$

        $C$ is supposed to be a conserved quantity in this system.
        - Chamecki (2005) eq. 15, which references Aref (1979)
        """
        n_vortons = self.n
        G = self.G
        C = 0
        for i, j in zip(*np.triu_indices(n_vortons, 1)):  # all combinations without repetition

            xi, yi = self.x[i], self.y[i]
            xj, yj = self.x[j], self.y[j]

            lij_sqd = (xi-xj)**2 + (yi-yj)**2

            Gi, Gj = G[i], G[j]

            C += Gi * Gj * lij_sqd

        return C


    def H(self):
        """Calculate $H$, the Hamiltonian of the system.

        $$
        H = -\frac{1}{4 \pi} \sum_{\alpha, \beta = 1; \alpha \neq \beta}^{N}
            \Gamma_{\alpha} \Gamma_{\beta}
            \ln | r_{\alpha} - r_{\beta} |
        $$
        """
        nv = self.n
        G = self.G
        r = self.state_mat  # vorton positions
        H = 0
        for a, b in zip(*np.triu_indices(nv, 1)):
            ra, rb = r[a], r[b]
            Ga, Gb = G[a], G[b]
            H += -1/(4*np.pi) * Ga * Gb * np.log(np.linalg.norm(ra - rb))

        return H


    def I(self):
        """Calculate $I$, the angular impulse of the system.

        $$
        I = \sum_{\alpha = 1}^{N} \Gamma_{\alpha} | r_{\alpha} |^2
        $$
        """
        G = self.G
        # r = self.state_mat
        x = self.x
        y = self.y

        # r_hat_sqd =

        return (G * (x**2 + y**2)).sum()


    # TODO: P and Q (coordinates of the center-of-vorticity)


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


    # Chamecki notes suggest this should be called "center of vorticity" or "linear impulse"
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


    def add_tracers(self, n, *, method="randu"):
        """Add `n` passive tracers (vortons with G=0)."""

        if method == "randu":
            d = 3.0  # displacement off center to select initial tracer coords from
            Gt = np.zeros((n,))
            xit = np.random.uniform(-d, d, (n,))  # TODO: optional centering choice (or cm?)
            yit = np.random.uniform(-d, d, (n,))

        else:
            raise NotImplementedError(f"method={method!r}")

        # TODO: spiral, circles, grid, randn, etc.

        xyt = np.column_stack((xit, yit))
        self.state_mat = np.append(self.state_mat, xyt, axis=0)
        self.G = np.append(self.G, Gt)



    # TODO: indexing dunder methods

    # TODO: class method to take List[Vorton] and return a Vortons?




def rotmat_2d(ang_deg):
    """Return rotation matrix for rotation `ang_deg` in degrees.
    For left-multiplication of a column position vector.

    Note: `scipy.spatial.transform.Rotation` can be used for 3-d rotations.
    """
    ang = np.deg2rad(ang_deg)
    c, s = np.cos(ang), np.sin(ang)
    R = np.array([
        [c, -s],
        [s, c]
    ])
    return R


def rotate_2d(x, *, ang_deg=None, rotmat=None):
    """Rotate vector `x` by `ang_deg` degrees.

    Either `ang_deg` or `rotmat` must be provided.

    Parameters
    ----------
    x : array-like (1-d)
        the vector to be rotated
    ang_deg : int, float
        degrees to rotate `x` about the origin
        positive -> counter-clockwise
    rotmat : array, shape (2, 2), optional
        rotation matrix -- left-multiplies a column position vector to give rotated position

    Optionally can pass `rotmat` instead to avoid computing it multiple times.
    """
    x = np.asarray(x)
    if ang_deg and rotmat:
        raise Exception("Only one of `ang_deg` and `rotmat` should be specified.")

    assert x.ndim == 1  # need a true vector

    if ang_deg:
        rotmat = rotmat_2d(ang_deg)
    else:
        if rotmat is None:
            raise Exception("If `ang_deg` is not provided, `rotmat` must be.")

    return (rotmat @ x[:, np.newaxis]).squeeze()


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
        verts[i, :] = rotate_2d(verts[i-1, :], rotmat=rotmat)

    return verts + c



if __name__ == "__main__":
    import matplotlib.pyplot as plt

    plt.close("all")

    vs = Vortons([1, 1], [0, 1], [0, 0])
    vs.add_tracers(10)
    vs.plot()

    # G sum here is 0, messing up the mom's...
    Vortons([1, -1], [0, 1], [0, 0]).plot()

    Vortons.regular_polygon(3).plot()

    Vortons.regular_polygon(10, c=(1, 0), r_c=0.5).plot()
