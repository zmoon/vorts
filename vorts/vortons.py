"""
`Vorton`/`Tracer` classes and `Vortons`/`Tracers` container classes.
"""
import functools
import inspect
from typing import NamedTuple
import warnings

import makefun
import numpy as np
import xarray as xr

from . import plot


# TODO: add global snippets dict but still allow `snippets` to override
def _add_snippets(func=None, *, snippets=None):
    """Decorator for adding snippets to a docstring. This function
    uses ``%(name)s`` substitution rather than `str.format` substitution so
    that the `snippets` keys can be invalid variable names.

    Based on [this one](https://github.com/lukelbd/proplot/blob/master/proplot/internals/docstring.py),
    but snippets passed as an argument instead of using a global dict.
    """
    if func is None:
        return functools.partial(_add_snippets, snippets=snippets)

    snippets = snippets if snippets else {}

    func.__doc__ = inspect.getdoc(func)
    if func.__doc__:
        func.__doc__ %= {key: value.strip() for key, value in snippets.items()}

    return func


class Tracers:
    """Collection of `Tracer`s."""
    def __init__(self, x, y):
        """
        Parameters
        ----------
        x, y : array_like
            shape: `(n_vortons,)`

            Tracer initial $x$ and $y$ positions.
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)

        assert x.shape == y.shape and x.ndim == 1

        # TODO: change to private attr, `_xy` or somesuch
        self.state_mat = np.column_stack((x, y))


    def __repr__(self):
        # unlike Vortons, might have many tracers
        # so don't need to show all in the repr
        n = self.n
        return f"Tracers(n={n})"

    @property
    def n(self):
        return self.state_mat.shape[0]

    @property
    def x(self):
        return self.state_mat[:,0]

    @property
    def y(self):
        return self.state_mat[:,1]

    def state_vec(self):
        return self.state_mat.T.flatten()

    def state_mat_full(self):
        """Full state mat for tracers doesn't include G."""
        warnings.warn("Note that `state_mat_full` for tracers is the same as `state_mat` (no G).")
        return self.state_mat

    # @classmethod
    # @_add_snippets(snippets=dict(params=_points_randu_params))
    # @makefun.with_signature(inspect.signature(points_randu))
    # def randu(cls, *args, **kwargs):
    #     """Return `Tracers` created by sampling from normal distributions, using `points_randu`.

    #     Parameters
    #     ----------
    #     %(params)s
    #     """
    #     return cls(*points_randu(*args, **kwargs).T)

    # @staticmethod
    # @functools.wraps(points_randu)
    # def randu(*args, **kwargs):
    #     xy = points_randu(*args, **kwargs).T
    #     return Tracers(*xy)

    def plot(self, *, connect=False, ax=None):
        """Plot tracers, connected if `connect=True`."""
        import matplotlib.pyplot as plt

        fig, ax = plot.maybe_new_figure(ax)

        x, y = self.x, self.y
        fmt = "-o" if connect else "o"
        ax.plot(x, y, fmt, c="0.5", ms=4, label="tracers")

        ax.set(
            xlabel="$x$",
            ylabel="$y$",
        )
        ax.set_aspect("equal", "box")
        fig.legend()
        ax.grid(True)
        fig.set_tight_layout(True)


def _extract_params_block(f):
    """Extract params block from `f`'s docstring."""
    lines = inspect.getdoc(f).splitlines()

    # Beginning of block, skipping header
    a = lines.index("Parameters") + 2

    # Find end of the block
    b = -1
    for i, line in enumerate(lines[a:]):
        if line.startswith("--"):  # new block header
            b = a + i - 2
            break

    return "\n".join(lines[a:b+1])


def _add_to_tracers(points_method=None, *, short=None):
    """Decorator for adding points fns to `Tracers`."""
    if points_method is None:
        return functools.partial(_add_to_tracers, short=short)

    short = short if short else ""
    params = _extract_params_block(points_method)

    @staticmethod
    @makefun.with_signature(inspect.signature(points_method))
    @_add_snippets(snippets=dict(params=params, short=short))
    def f(*args, **kwargs):
        """%(short)s

        Parameters
        ----------
        %(params)s

        Returns
        -------
        Tracers
        """
        return Tracers(*points_method(*args, **kwargs).T)

    # Add method to `Tracers`, removing the `points_` part of the name
    setattr(Tracers, f"{points_method.__name__[7:]}", f)

    return points_method


_points_returns = """
numpy.ndarray
    2-d array with first column $x$ and second column $y$.
""".strip()


@_add_to_tracers(short="Create `Tracers` by sampling from uniform random distributions using `points_randu`.")
@_add_snippets(snippets=dict(returns=_points_returns))
def points_randu(n, *, c=(0, 0), dx=2, dy=2):
    """Sample from 2-d uniform.

    Parameters
    ----------
    n : int
        Number of points.
    c : array_like
        Coordinates of the center ($x_c$, $y_c$).
    dx, dy : float
        $x$ positions will be sampled from $[$`-dx`, `dx`$)$, and $y$ similarly.

    Returns
    -------
    %(returns)s
    """
    c = np.asarray(c)
    x = np.random.uniform(-dx, dx, (n,))
    y = np.random.uniform(-dy, dy, (n,))
    return np.column_stack((x, y)) + c


@_add_to_tracers(short="Create spiral arrangement of `Tracers` using `points_spiral`.")
@_add_snippets(snippets=dict(returns=_points_returns))
def points_spiral(n, *, c=(0, 0), rmin=0, rmax=2, revs=3):
    """Create spiral of points.

    Parameters
    ----------
    n : int
        Number of points.
    c : array_like
        Coordinates of the center ($x_c$, $y_c$).
    rmin : float
        Minimum radius (distance from the center for the innermost point).
    rmax : float
        Maximum radius (distance from the center for the outermost point).
    revs : float
        Total number of revolutions in the spiral.

    Returns
    -------
    %(returns)s
    """
    # TODO: option for linear distance between consecutive points
    c = np.asarray(c)

    rad = np.linspace(rmin, rmax, n)  # radius

    deg_tot = revs*360
    rotmat = rotmat_2d(deg_tot/n)
    rhat = np.full((n, 2), (0, 1), dtype=float)  # rhat: unit vectors
    for i in range(1, n):
        rhat[i, :] = rotate_2d(rhat[i-1, :], rotmat=rotmat)
    # TODO: here would be simpler to do polar coords first then convert to x,y

    return rad[:, np.newaxis] * rhat + c


@_add_to_tracers(short="Create `Tracers` by sampling from normal distributions using `points_randn`.")
@_add_snippets(snippets=dict(returns=_points_returns))
def points_randn(n, *, mu_x=0, mu_y=0, sig_x=1, sig_y=1, c=(0, 0)):
    """Sample from normal distribution.

    Parameters
    ----------
    n : int
        Number of points.
    mu_x, mu_y : float
        Mean/center of the distribution in each direction.
    sig_x, sig_y : float
        Standard deviation of the distribution in each direction.
    c : array_like
        Coordinates of the center ($x_c$, $y_c$).

    Returns
    -------
    %(returns)s
    """
    c = np.asarray(c)
    x = np.random.normal(mu_x, sig_x, (n,))
    y = np.random.normal(mu_y, sig_y, (n,))
    return np.column_stack((x, y)) + c


# TODO: sample from any scipy dist, optionally different for x and y


@_add_to_tracers(short="Create gridded arrangement of `Tracers` using `points_grid`.")
@_add_snippets(snippets=dict(returns=_points_returns))
def points_grid(nx, ny, *, xbounds=(-2, 2), ybounds=(-2, 2), c=(0, 0)):
    """Points on a grid.

    Parameters
    ----------
    nx, ny : int
        Number of points in the grid in each direction.
    xbounds, ybounds : array_like
        Inclusive bounds in each direction (lower, upper).
    c : array_like
        Coordinates of the center ($x_c$, $y_c$).

    Returns
    -------
    %(returns)s
    """
    c = np.asarray(c)
    x = np.linspace(*xbounds, nx)
    y = np.linspace(*ybounds, ny)
    X, Y = np.meshgrid(x, y)
    return np.column_stack((X.ravel(), Y.ravel())) + c


@_add_to_tracers(short="Create concentric circle arrangement of `Tracers` using `points_circles`.")
@_add_snippets(snippets=dict(returns=_points_returns))
def points_circles(ns=(10, 20, 34, 50), rs=(0.5, 1, 1.5, 2), *, c=(0, 0)):
    """Concentric circles.

    Parameters
    ----------
    ns : array_like
        Number of points in each circle.
    rs : array_like
        Radii of each circle (one for each value of `ns`).
    c : array_like
        Coordinates of the center ($x_c$, $y_c$).

    Returns
    -------
    %(returns)s
    """
    c = np.asarray(c)
    x = []
    y = []
    for n, r in zip(ns, rs):
        dtheta = 360 / n
        thetas = np.deg2rad(np.linspace(0, 360-dtheta, n))
        x = np.append(x, r*np.cos(thetas))
        y = np.append(y, r*np.sin(thetas))

    return np.column_stack((x, y)) + c


def rotmat_2d(ang_deg):  # TODO: could lru_cache?
    """Return rotation matrix for rotation `ang_deg` in degrees.
    For left-multiplication of a column position vector.

    .. note::
       [`scipy.spatial.transform.Rotation`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html)
       can be used for 3-d rotations.
    """
    ang = np.deg2rad(ang_deg)
    c, s = np.cos(ang), np.sin(ang)
    R = np.array([
        [c, -s],
        [s, c]
    ])
    return R


def rotate_2d(x, *, ang_deg=None, rotmat=None):
    r"""Rotate vector `x` by `ang_deg` degrees.

    .. important::
       Either `ang_deg` or `rotmat` can be provided to specify the degree of rotation, but not both.

       If `ang_deg` is used, the rotation matrix will be computed with `rotmat_2d`, so
       you can pass `rotmat` instead to avoid computing it multiple times.

    Parameters
    ----------
    x : array_like
        The vector to be rotated.
    ang_deg : int, float
        Degrees by which to rotate `x` about the origin.

        positive $\to$ counter-clockwise rotation
    rotmat : array_like
        shape: `(2, 2)`

        Rotation matrix -- left-multiplies a column position vector to give rotated position.

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


def vertices_regular_polygon(n, *, c=(0, 0), r_c=1):
    """Regular polygon vertices.

    Parameters
    ----------
    n : int
        Polygon order (number of sides/vertices).
    c : array_like
        Coordinates of the center of the inscribing circle ($x_c$, $y_c$).
    r_c : float, int
        Radius $r_c$ of the inscribing circle.
    """
    c = np.asarray(c)

    # initial vertex
    vert0 = np.r_[0, r_c]

    # rotation matrix -- left-multiplies a column position vector to give rotated position
    rotmat = rotmat_2d(360/n)

    verts = np.full((n, 2), vert0, dtype=float)
    # successive rotations
    for i in range(1, n):
        verts[i, :] = rotate_2d(verts[i-1, :], rotmat=rotmat)

    return verts + c


def vertices_isos_triangle(*, theta_deg=None, Lambda=None):
    r"""Isosceles triangle vertices.
    With fixed top point $(0, 1)$ and fixed left & right $y=-0.5$.

    .. important::
       Either `theta_deg` or `Lambda` can be used to specify the angle, but not both.

    Parameters
    ----------
    theta_deg : float
        Value of the two angles $\theta$ between the horizontal base and connections to the top point at $(0,1)$
        in degrees.

        $\theta = 72^{\circ} \to \Lambda_c$ (equal to $1/\sqrt{2}$)

        $\theta = 60^{\circ} \to$ equilateral triangle (can also create with `vertices_regular_polygon`,
        which gives control over size and location)

    Lambda : float
        $\Lambda \in (0, 1]$. Related to $\theta$ by $\theta = \pi / (\Lambda^2 + 2)$

        $\Lambda = 1 \to$ equilateral triangle

    """
    if (theta_deg is not None and Lambda is not None) or (theta_deg is None and Lambda is None):
        raise Exception("Specify either `theta_deg` or `Lambda` (not both).")

    if Lambda:
        assert Lambda > 0 and Lambda <= 1
        theta_deg = 180 / (Lambda**2 + 2)

    theta = np.deg2rad(theta_deg)

    xb = 1.5/np.tan(theta)  # one half of x base

    xi = [-xb,  0,  xb]
    yi = [-0.5, 1, -0.5]

    Lambda = np.sqrt( (180-2*theta_deg) / float(theta_deg) )  # Marcelo eqns 17--19

    return np.column_stack((xi, yi))


class Vorton(NamedTuple):
    """A vorton that knows its current state (position and strength)."""
    G: float
    r"""$\Gamma$, the strength of the circulation, with sign to indicate direction.

    See also
    --------
    Vortons : For a more detailed description.
    """
    x: float
    """$x$ position"""
    y: float
    """$y$ position"""


class Tracer(NamedTuple):
    r"""Tracer -- a vorton with $\Gamma=0$ (no circulation/mass) that knows its current position."""
    x: float
    """$x$ position"""
    y: float
    """$y$ position"""


# TODO: PointVortices ABC that implements adding, has position state_mat, n, x, y, state_vec, etc.
#       Vortons and Tracers could both be based on it
#       also should add xy (state_mat for both) and xy_vec (state_vec)


# could exchange x,y for r at some point, to open 3-d option more easily
class Vortons:
    """Collection of `Vorton`s."""
    def __init__(self, G, x, y):
        r"""

        Parameters
        ----------
        G, x, y : array_like
            shape: `(n_vortons,)`

            `G`: $\Gamma$s ("G" for [Gamma](https://en.wikipedia.org/wiki/Gamma)).

            $\Gamma$ represents the strength of the circulation, with sign to indicate direction.
            In fluid dynamics, circulation $\Gamma$ is the line integral of velocity
            or flux of vorticity vectors through a surface (here the $xy$-plane).

            `x`: $x$ positions

            `y`: $y$ positions

        """
        self.G = np.asarray(G)
        r"""Array of vorton strengths ($\Gamma$)."""
        if np.any(self.G == 0):
            warnings.warn(
                "Tracers should be in a `Tracers` instance. "
                "The ability to add them here may be removed in the future."
            )

        # the state matrix has shape (n_vortons, n_pos_dims) (G excluded since time-invariant)
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        self.state_mat = np.column_stack((x, y))
        """2-d array of $(x, y)$ coordinates -- each row is the coordinate of one vorton."""

        assert self.G.ndim == 1 and self.state_mat.shape[1] == 2
        assert self.G.size == self.state_mat.shape[0]  # n_vortons

        # create initial corresponding Vorton objects
        self._update_vortons()


    # these 2 don't really need to be property?
    # maybe shouldn't be, to emphasize that state_mat is the real data
    # @property
    def state_vec(self):
        """Return flattened state matrix (`Vortons.state_mat`; `Vortons.G` not included).

        If using `vorts.py.integ.integrate_scipy`,
        this is needed to feed to [`scipy.integrate.solve_ivp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html),
        which requires a 1-d array for its `y0` input.
        """
        return self.state_mat.T.flatten()  # TODO: change to ravel, to return view when possible

    # @property
    def state_mat_full(self):
        """Return full state matrix: (`Vortons.G`, `Vortons.x`, `Vortons.y`) as 3 columns."""
        return np.column_stack((self.G, self.state_mat))

    # seems to return a view into self.G, so ok to be property
    @property
    def G_col(self):
        """`Vortons.G` as a column vector."""
        return self.G[:, np.newaxis]

    @property
    def x(self):
        """Array of $x$ positions (a view into `Vortons.state_mat`)."""
        # slice indexing should give just a view into `self.state_mat`
        # thus `self.x.base` will return the state mat
        # i.e., `vs.x.base is vs.state_mat`
        return self.state_mat[:,0]

    @property
    def y(self):
        """Array of $y$ positions (a view into `Vortons.state_mat`)."""
        return self.state_mat[:,1]

    @property
    def n(self):
        """Number of vortons."""
        # numpy.ndarray size lookups are esentially free
        return self.G.size  # will have to change if want to allow single G at some point

    def __repr__(self):
        # n_vortons shouldn't be too many, so let's show all
        s_vorts = "\n".join(f"  {v}" for v in self._vortons)
        return f"Vortons(\n{s_vorts}\n)"
        # TODO: should this call `self._update_vortons`? so as to not be out-of-date if state changes?

    def _update_vortons(self):
        self._vortons = [Vorton(G, x, y) for G, x, y in self.state_mat_full()]


    def C(self):
        r"""Calculate $C$.

        $$
        C = \sum_{\alpha, \beta = 1; \alpha \neq \beta}^{N}
            \Gamma_{\alpha} \Gamma_{\beta} l_{\alpha \beta}^{2}
        $$

        $C$ is supposed to be a conserved quantity in this system.
        -- Chamecki (2005) eq. 15, which references Aref (1979)
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
        r"""Calculate $H$, the Hamiltonian of the system.

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
        r"""Calculate $I$, the angular impulse of the system.

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


    # TODO: results are not right for equi tri... need to check formulas
    def theta(self):
        r"""Calculate $\theta$, the action angles??

        Chamecki eq. 19
        """
        N = self.n
        I = self.I()
        H = self.H()

        return (2/(N-1))**(N*(N-1)/2) * I**(N*(N-1)) * np.exp(4*np.pi*H)


    def plot(self):
        """Plot the vortons.
        (Only their current positions, which are all `Vortons` knows about.)
        """
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()

        # plot vorton positions
        c_Gp = "cadetblue"
        c_Gm = "salmon"
        G = self.G
        Gp, Gm = G > 0, G < 0
        x, y = self.x, self.y
        ax.plot(x[Gp], y[Gp], "o", ms=7, c=c_Gp, label=r"$\Gamma > 0$")
        ax.plot(x[Gm], y[Gm], "o", ms=7, c=c_Gm, label=r"$\Gamma < 0$")

        # plot center of mass
        x_cm, y_cm = self.cm()
        s_cm = f"({x_cm:.4g}, {y_cm:.4g})"
        ax.plot(x_cm, y_cm, "*", ms=13, c="gold", label=f"center of mass\n{s_cm}")

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
        ax.set_aspect("equal", "datalim")
        fig.legend()
        ax.grid(True)
        fig.tight_layout()

        # return


    def mom(self, n, *, abs_G=False, center=False):
        r"""Compute `n`-th moment.

        Parameters
        ----------
        n : int
            Which [moment](https://en.wikipedia.org/wiki/Moment_(mathematics)) to calculate.
        abs_G : bool
            Whether to take the absolute value of the $\Gamma$ values (false by default).
        center : bool
            `True`: evaluate moment wrt. center of mass from `Vortons.cm`

            `False`: evaluate moment wrt. $(0, 0)$
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
    def center_of_mass(self):
        r"""Compute [center of mass](https://en.wikipedia.org/wiki/Center_of_mass#A_system_of_particles)
        using $\Gamma$ (`Vortons.G`) as mass.
        Equivalent to `Vortons.mom` with `n=1`, `abs_G=True` (currently), `center=False`.
        """
        # TODO: what impact should sign of G have on cm? mass is always pos. but G can be neg.
        return self.mom(1, abs_G=True, center=False)

    def cm(self):
        """Alias for `Vortons.center_of_mass`."""
        return self.center_of_mass()

    def center_coords(self, inplace=False):
        """Make $(0, 0)$ the center of mass."""
        xy_cm = self.cm()
        x_cm, y_cm = xy_cm
        if not inplace:
            return Vortons(self.G, self.x-x_cm, self.y-y_cm)
        else:
            self.state_mat -= x_cm


    @classmethod
    def regular_polygon(cls, n, *, G=None, **kwargs):
        r"""Create Vortons with positions corresponding to the vertices of a regular polygon.

        Parameters
        ----------
        n : int
            Polygon order.
        G : int, array_like, optional
            $\Gamma$ value(s) to use.

            Single value or size-$n$ array-like of values.

            default: 1.0

        **kwargs
            Passed on to `vertices_regular_polygon`.
        """
        G = _maybe_fill_G(G, n)
        xy = vertices_regular_polygon(n, **kwargs).T  # x, y cols-> rows (for unpacking)
        return cls(G, *xy)


    @classmethod
    def isos_triangle(cls, *, G=None, **kwargs):
        r"""Create Vortons with isosceles triangle vertices.

        Parameters
        ----------
        G : int, array_like, optional
            $\Gamma$ value(s) to use.

            Single value or size-3 array-like of values.

            default: 1.0

        `**kwargs`
            Passed on to `vertices_isos_triangle`.
        """
        G = _maybe_fill_G(G, 3)
        xy = vertices_isos_triangle(**kwargs).T
        return cls(G, *xy)


    def maybe_with_tracers(self, tracers: Tracers = None):
        """Return new `Vortons` with the tracers (maybe) included.

        .. caution::
           Temporary? hack to get full (combined) `state_vec` for the whole system.

        If `Tracers` is `None`, just returns `self`.
        """
        if tracers is None:
            return self

        # TODO: add option to extend existing Vortons instead of creating a new. `inplace`?

        G_v = self.G
        xy_v = self.state_mat
        xy_t = tracers.state_mat
        G_t = np.zeros((tracers.n,))

        G = np.append(G_v, G_t)
        x, y = np.append(xy_v, xy_t, axis=0).T  # unpacking arrays goes by rows

        return Vortons(G, x, y)

    # TODO: indexing dunder methods

    # TODO: overload addition and such

    # TODO: class method to take List[Vorton] and return a Vortons?



def _maybe_fill_G(G, n):
    if G is None:  # this first part maybe shouldn't be here? or kwarg for default G val?
        G = 1.0
    G = np.asarray(G)
    if G.size == 1:  # single G provided, or using the default
        G = np.full((n,), G)  # TODO: could also the constructor to accept single G
    if G.size != n:
        raise ValueError(f"`G` must have size `n` or 1, but is {G.size!r}")

    return G


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    plt.close("all")

    vs = Vortons([1, 1], [0, 1], [0, 0])
    vs.plot()

    # G sum here is 0, messing up the mom's...
    Vortons([1, -1], [0, 1], [0, 0]).plot()

    Vortons.regular_polygon(3).plot()

    Vortons.regular_polygon(10, c=(1, 0), r_c=0.5).plot()

    Vortons.isos_triangle(theta_deg=72).plot()

    Vortons.isos_triangle(Lambda=0.49).plot()

    ts = Tracers.randu(50)

    Tracers.spiral(100).plot()

    Tracers.spiral(200, c=(1, 0), revs=10).plot(connect=True)
