"""
`Vorton`/`Tracer` classes and `Vortons`/`Tracers` container classes.

The module name is "vortons" because most of the focus is on the vorton collection.

`Vortons` and `Tracers` can combined using `+`.
Currently the result has the class of the one on the left in the addition.
>>> import vorts
>>> ts = vorts.Tracers(0, 0)
>>> ts + ts
Tracers(
  Tracer(x=0.0, y=0.0)
  Tracer(x=0.0, y=0.0)
)
>>> vorts.Vortons.regular_polygon(3) + ts
Vortons(
  Vorton(G=1.0, x=0.0, y=1.0)
  Vorton(G=1.0, x=-0.8660254037844387, y=-0.4999999999999998)
  Vorton(G=1.0, x=0.8660254037844384, y=-0.5000000000000003)
  Vorton(G=0.0, x=0.0, y=0.0)
)

`Vortons`/`Tracers` can also be transformed (creating a new object).
>>> ts = vorts.Tracers(2, 1)
>>> ts
Tracers(
  Tracer(x=2.0, y=1.0)
)
>>> ts + (1, 1)  # translate
Tracers(
  Tracer(x=3.0, y=2.0)
)
>>> 2 * ts  # scale
Tracers(
  Tracer(x=4.0, y=2.0)
)
>>> ts.rotate(90)  # rotate
Tracers(
  Tracer(x=-0.9999999999999999, y=2.0)
)
"""
import abc
import functools
import inspect
import warnings
from typing import NamedTuple

import makefun
import numpy as np

from .plot import _maybe_new_fig

_SNIPPETS = {}


def _add_snippets(func=None, *, snippets=None):
    """Decorator for adding snippets to a docstring. This function
    uses ``%(name)s`` substitution rather than `str.format` substitution so
    that the `snippets` keys can be invalid variable names.

    Based on [this one](https://github.com/lukelbd/proplot/blob/master/proplot/internals/docstring.py),
    but snippets passed as an argument instead of using a global dict.
    """
    if func is None:
        return functools.partial(_add_snippets, snippets=snippets)

    if snippets is None:
        snippets = {}

    snippets = {**_SNIPPETS, **snippets}

    func.__doc__ = inspect.getdoc(func)
    if func.__doc__:
        func.__doc__ %= {key: value.strip() for key, value in snippets.items()}

    return func


# class PointBase(NamedTupleMeta):
#     """Point base class with $x$ and $y$."""
#     x: float
#     """$x$ position"""
#     y: float
#     """$y$ position"""


class Vorton(NamedTuple):
    """A vorton that knows its current state (position and strength).

    See also
    --------
    Vortons : For a more detailed description.
    """

    G: float
    r"""$\Gamma$, the strength of the circulation, with sign to indicate direction."""
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


class PointsBase(abc.ABC):
    """Points base class with $x$ and $y$."""

    def __init__(self, x, y):
        """
        Parameters
        ----------
        x, y : array_like
            shape: `(n_points,)`

            Initial $x$ and $y$ positions.
        """
        x = np.atleast_1d(np.asarray(x, dtype=float))  # should be view if `x` is ndarray
        y = np.atleast_1d(np.asarray(y, dtype=float))
        assert x.shape == y.shape and x.ndim == 1

        self._xy = np.column_stack((x, y))

        self._points = None

    @abc.abstractmethod
    def _update_points(self):
        """Update `_points` list of corresponding point objects."""
        ...

    def __repr__(self):
        self._update_points()  # ensure consistency
        n_show = min(len(self._points), 10)
        s_points = "\n".join(f"  {v}" for v in self._points[:n_show])
        if n_show < self.n:
            s_points += "\n  ..."
        return f"{self.__class__.__name__}(\n{s_points}\n)"

    @property
    def n(self):
        """Number of points."""
        return self._xy.shape[0]

    @property
    def x(self):
        """Array of $x$ positions (should be a view)."""
        return self._xy[:, 0]

    @property
    def y(self):
        """Array of $y$ positions (should be a view)."""
        return self._xy[:, 1]

    @property
    def xy(self):
        """2-d array of $(x, y)$ coordinates -- each row is the coordinate of one point.
        This is the data array on which the others depend.
        """
        return self._xy

    @xy.setter
    def xy(self, xy):
        warnings.warn("The coordinates are not intended to be modified this way. Doing nothing.")
        # Elements can still be modified though! And through the other views to `_xy` as well.

    @abc.abstractmethod
    def state_mat_full(self):
        """Full state matrix (could be same as `xy` but should return a copy).
        Columns should be in the same order as the class init positional parameters.
        """
        ...

    @abc.abstractmethod
    def plot(self):
        """Plot state."""
        ...

    def __add__(self, other):
        if hasattr(other, "xy"):  # other points collection
            xy = np.append(self.xy, other.xy, axis=0)
        else:  # vector for translation?
            try:
                xyp = np.asarray(other)
                assert xyp.shape == (2,)
            except (TypeError, AssertionError) as e:
                raise TypeError(f"{other!r} is unsuitable for adding to {type(self)}.") from e
            else:
                xy = self.xy + xyp

        return self.__class__(*xy.T)

    # def __iadd__(self, other):

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            xy = self.xy * other
            return self.__class__(*xy.T)
        else:
            raise TypeError(f"Multiplication by {type(other)} is unsupported.")

    def __rmul__(self, other):
        return self.__mul__(other)

    def rotate(self, theta, *, units="deg", inplace=False):
        """Rotate coordinates about the origin by angle `theta` (units `'rad'` or `'deg'`)."""
        if inplace:
            raise NotImplementedError
        if units not in ("rad", "deg"):
            raise ValueError
        theta_deg = theta if units == "deg" else np.rad2deg(theta)
        rotmat = rotmat_2d(theta_deg)
        xy = self.xy.copy()
        for i in range(xy.shape[0]):  # over rows
            xy[i] = rotate_2d(xy[i], rotmat=rotmat)

        return self.__class__(*xy.T)


class Tracers(PointsBase):
    """Collection of `Tracer`s."""

    def __init__(self, x, y):
        """
        Parameters
        ----------
        x, y : array_like
            shape: `(n_tracers,)`

            Tracer initial $x$ and $y$ positions.
        """
        super().__init__(x=x, y=y)

    def _update_points(self):
        self._points = [Tracer(x, y) for x, y in self._xy]

    @property
    def tracers(self):
        """List of `Tracer` instances corresponding to the coordinates.

        .. warning::
           Modifying this will not update the `Tracers` data.
        """
        self._update_points()  # ensure consistency
        return self._points

    def state_mat_full(self):
        """Full state mat for tracers doesn't include G."""
        # warnings.warn("Note that `state_mat_full` for tracers is the same as `state_mat` (no G).")
        return self._xy.copy()

    def plot(self, *, connect=False, adjustable="box", ax=None, **kwargs):
        """Plot tracers, with points connected if `connect=True`."""
        fig, ax = _maybe_new_fig(ax=ax, **kwargs)

        x, y = self.x, self.y
        fmt = "-o" if connect else "o"
        ax.plot(x, y, fmt, c="0.5", ms=4, label="tracers")

        ax.set(
            xlabel="$x$",
            ylabel="$y$",
        )
        ax.set_aspect("equal", adjustable)
        fig.legend()
        ax.grid(True)
        fig.set_tight_layout(True)


# Note: could exchange x,y for r at some point, to open 3-d option more easily
class Vortons(PointsBase):
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
        super().__init__(x=x, y=y)

        self.G = np.atleast_1d(np.asarray(G, dtype=float))
        r"""Array of vorton strengths ($\Gamma$)."""
        # if np.any(self.G == 0):
        #     warnings.warn(
        #         "Tracers should be in a `Tracers` instance. "
        #         "The ability to add them here may be removed in the future."
        #     )
        assert self.G.ndim == 1 and self.G.size == self.n  # n_vortons

        # the state matrix has shape (n_vortons, n_pos_dims) (G excluded since time-invariant)
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        self.state_mat = np.column_stack((x, y))
        """2-d array of $(x, y)$ coordinates -- each row is the coordinate of one vorton."""

    def _update_points(self):
        self._points = [Vorton(G, x, y) for G, x, y in self.state_mat_full()]

    def vortons(self):
        """List of corresponding `Vorton` objects."""
        self._update_points()
        return self._points

    def state_mat_full(self):
        """Return full state matrix: ($G$, $x$, $y$ / `Vortons.G`, `Vortons.x`, `Vortons.y`) as 3 columns."""
        return np.column_stack((self.G, self.xy))

    # Seems to return a view into self.G, so ok to be property
    @property
    def G_col(self):
        """`Vortons.G` as a column vector."""
        return self.G[:, np.newaxis]

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

            lij_sqd = (xi - xj) ** 2 + (yi - yj) ** 2

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
            H += -1 / (4 * np.pi) * Ga * Gb * np.log(np.linalg.norm(ra - rb))

        return H

    def I(self):  # noqa: 743,741
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
        I = self.I()  # noqa: 741
        H = self.H()

        # fmt: off
        return (2/(N-1))**(N*(N-1)/2) * I**(N*(N-1)) * np.exp(4*np.pi*H)
        # fmt: on

    def plot(self, *, ax=None, adjustable="datalim", **kwargs):
        """Plot the vortons.
        (Only their current positions, which are all `Vortons` knows about.)
        """
        fig, ax = _maybe_new_fig(ax=ax, **kwargs)

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
        ax.plot(x_cm, y_cm, "o", ms=13, c="gold", label=f"center of mass\n{s_cm}")

        # 2nd mom
        x_cm2, y_cm2 = self.moment(2)
        s_cm2 = f"({x_cm2:.4g}, {y_cm2:.4g})"
        ax.plot(x_cm2, y_cm2, "*", ms=13, c="0.4", label=f"mom2\n{s_cm2}")

        # 3nd mom
        # TODO: helper fn to DRY this
        x_cm3, y_cm3 = self.moment(3)
        s_cm3 = f"({x_cm3:.4g}, {y_cm3:.4g})"
        ax.plot(x_cm3, y_cm3, "*", ms=13, c="0.55", label=f"mom3\n{s_cm3}")

        ax.set(
            title=f"$C = {self.C():.4g}$",
            xlabel="$x$",
            ylabel="$y$",
        )
        ax.set_aspect("equal", adjustable)
        fig.legend()
        ax.grid(True)
        fig.tight_layout()

    def moment(self, n, *, abs_G=False, center=False):
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

        x_mom = (G * (x - c) ** n).sum(axis=0) / G_tot  # sum along vortons dim, giving a position
        # ^ maybe this should be x - x_cm here...

        return x_mom

    # Chamecki notes suggest this should be called "center of vorticity" or "linear impulse"
    def center_of_mass(self):
        r"""Compute [center of mass](https://en.wikipedia.org/wiki/Center_of_mass#A_system_of_particles)
        using $\Gamma$ (`Vortons.G`) as mass.
        Equivalent to `Vortons.moment` with `n=1`, `abs_G=True` (currently), `center=False`.
        """
        # TODO: what impact should sign of G have on cm? mass is always pos. but G can be neg.
        return self.moment(1, abs_G=True, center=False)

    def cm(self):
        """Alias for `Vortons.center_of_mass`."""
        return self.center_of_mass()

    def center_coords(self, inplace=False):
        """Make $(0, 0)$ the center of mass."""
        xy_cm = self.cm()
        x_cm, y_cm = xy_cm
        if not inplace:
            return Vortons(self.G, self.x - x_cm, self.y - y_cm)
        else:
            self.state_mat -= x_cm

    def _add_vortons(self, vortons, inplace=False):
        if inplace:
            raise NotImplementedError
        Gxy = np.append(self.state_mat_full(), vortons.state_mat_full(), axis=0)
        return self.__class__(*Gxy.T)

    def _maybe_add_tracers(self, tracers, inplace=False):
        if tracers is None:
            return self
        if inplace:
            raise NotImplementedError
        G = np.append(self.G, np.zeros((tracers.n,)))
        x, y = np.append(self.xy, tracers.xy, axis=0).T
        return self.__class__(G, x, y)

    def to_tracers(self):
        """Return `Tracers` instance corresponding to the vorton positions."""
        return Tracers(self.x, self.y)

    # Overriding base class so can treat tracers and vortons differently (and due to G)
    def __add__(self, other):
        if isinstance(other, self.__class__):
            return self._add_vortons(other)
        elif isinstance(other, Tracers):
            return self._maybe_add_tracers(other)
        else:  # try translation
            xy = (self.to_tracers() + other).xy
            return self.__class__(self.G, *xy.T)

    # def __iadd__

    # Overriding base class due to G
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            xy = self.xy * other
            return self.__class__(self.G, *xy.T)
        else:  # keep message in sync with `PointsBase.__mul__`
            raise TypeError(f"Multiplication by {type(other)} is unsupported.")

    # Overriding base class due to G
    def rotate(self, theta, *, units="deg", inplace=False):
        """Rotate coordinates about the origin by angle `theta` (units `'rad'` or `'deg'`)."""
        xy = self.to_tracers().rotate(theta, units=units, inplace=inplace).xy
        return self.__class__(self.G, *xy.T)

    # TODO: indexing dunder methods

    # TODO: class method to take List[Vorton] and return a Vortons?


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

    return "\n".join(lines[a : b + 1])


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


_G_params = r"""
G : float, array_like, optional
    $\Gamma$ value(s) to use.

    Single value or array-like vector of values.

    default: 1
""".strip()


def _add_to_vortons(points_method=None, *, short=None):
    """Decorator for adding points fns to `Vortons`."""
    if points_method is None:
        return functools.partial(_add_to_vortons, short=short)

    short = short if short else ""
    params = _extract_params_block(points_method)

    # Add `G` as first param in the signature
    sig0 = inspect.signature(points_method)
    sig = makefun.add_signature_parameters(
        sig0, last=inspect.Parameter("G", kind=inspect.Parameter.KEYWORD_ONLY, default=1)
    )

    @staticmethod
    @makefun.with_signature(sig)
    @_add_snippets(snippets=dict(params=params, short=short, G_params=_G_params))
    def f(*args, **kwargs):
        """%(short)s

        Parameters
        ----------
        %(params)s
        %(G_params)s

        Returns
        -------
        Vortons
        """
        G = kwargs.pop("G")
        xy = points_method(*args, **kwargs)
        n = xy.shape[0]
        G_ = _maybe_fill_G(G, n)
        return Vortons(G_, *xy.T)

    # Add method to `Vortons`
    name0 = points_method.__name__
    name = name0[name0.index("_") + 1 :]
    setattr(Vortons, name, f)

    return points_method


_points_returns = """
numpy.ndarray
    2-d array with first column $x$ and second column $y$.
""".strip()


@_add_to_tracers(
    short="Create `Tracers` by sampling from uniform random distributions using `points_randu`."
)
@_add_snippets(snippets=dict(returns=_points_returns))
def points_randu(n, *, dx=2, dy=2):
    """Sample from 2-d uniform.

    Parameters
    ----------
    n : int
        Number of points.
    dx, dy : float
        $x$ positions will be sampled from $[$`-dx`, `dx`$)$, and $y$ similarly.

    Returns
    -------
    %(returns)s
    """
    x = np.random.uniform(-dx, dx, (n,))
    y = np.random.uniform(-dy, dy, (n,))
    return np.column_stack((x, y))


@_add_to_tracers(short="Create spiral arrangement of `Tracers` using `points_spiral`.")
@_add_snippets(snippets=dict(returns=_points_returns))
def points_spiral(n, *, rmin=0, rmax=2, revs=3, kind="Archimedean", spacing="linear"):
    r"""Create spiral of points.

    Parameters
    ----------
    n : int
        Number of points.
    rmin : float
        Minimum radius (distance from the center for the innermost point).
        Normally should be 0 (not really a spiral without the 0 point).
    rmax : float
        Maximum radius (distance from the center for the outermost point).
    revs : float
        Total number of revolutions in the spiral.
    kind : str, {'Archimedean', "Fermat's", 'logarithmic'}
        Type of spiral.
    spacing : str, {'linear', 'log', 'inv-exp', '1/x'}
        Method for spacing $\theta$ values on the spiral.

    Returns
    -------
    %(returns)s
    """
    x = np.linspace(0, 1, n)
    if spacing == "linear":
        x2 = x[1:]
    elif spacing == "log":
        x2 = np.log(x[1:] * 10 + 1) / 10
    elif spacing == "inv-exp":
        x2 = 1 - np.exp(-x[1:])
    elif spacing == "1/x":
        x2 = 6 - 1 / (x[1:] + 1 / 6)
    else:
        raise ValueError
    theta = np.append(0, np.deg2rad(x2 / x2.max() * 360 * revs))

    # Compute r from theta
    if kind == "Archimedean":
        r0 = theta / theta.max()
    elif kind == "Fermat's":
        r0 = np.sqrt(theta / theta.max())
    elif kind == "logarithmic":
        raise NotImplementedError
    else:
        raise ValueError
    r = rmin + r0 * (rmax - rmin)

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    return np.column_stack((x, y))


@_add_to_tracers(
    short="Create `Tracers` by sampling from normal distributions using `points_randn`."
)
@_add_snippets(snippets=dict(returns=_points_returns))
def points_randn(n, *, mu_x=0, mu_y=0, sig_x=1, sig_y=1):
    """Sample from normal distribution.

    Parameters
    ----------
    n : int
        Number of points.
    mu_x, mu_y : float
        Mean/center of the distribution in each direction.
    sig_x, sig_y : float
        Standard deviation of the distribution in each direction.

    Returns
    -------
    %(returns)s
    """
    x = np.random.normal(mu_x, sig_x, (n,))
    y = np.random.normal(mu_y, sig_y, (n,))
    return np.column_stack((x, y))


# TODO: sample from any scipy dist, optionally different for x and y


@_add_to_tracers(short="Create gridded arrangement of `Tracers` using `points_grid`.")
@_add_snippets(snippets=dict(returns=_points_returns))
def points_grid(nx, ny, *, xbounds=(-2, 2), ybounds=(-2, 2), dxy=None):
    """Points on a grid.

    Parameters
    ----------
    nx, ny : int
        Number of points in the grid in each direction.
    xbounds, ybounds : array_like
        Inclusive bounds in each direction (lower, upper).
    dxy : float, optional
        Overrides `xbounds` and `ybounds`, setting both to `(-dxy, dxy)`;
        more convenient if finer-grained control is not needed.

    Returns
    -------
    %(returns)s
    """
    if dxy:
        xbounds = ybounds = (-dxy, dxy)
    x = np.linspace(*xbounds, nx)
    y = np.linspace(*ybounds, ny)
    X, Y = np.meshgrid(x, y)
    return np.column_stack((X.ravel(), Y.ravel()))


@_add_to_tracers(short="Create concentric circle arrangement of `Tracers` using `points_circles`.")
@_add_snippets(snippets=dict(returns=_points_returns))
def points_circles(ns=(10, 20, 34, 50), rs=(0.5, 1, 1.5, 2)):
    """Concentric circles.

    Parameters
    ----------
    ns : array_like
        Number of points in each circle.
    rs : array_like
        Radii of each circle (one for each value of `ns`).

    Returns
    -------
    %(returns)s
    """
    x = []
    y = []
    for n, r in zip(ns, rs):
        dtheta = 360 / n
        thetas = np.deg2rad(np.linspace(0, 360 - dtheta, n))
        x = np.append(x, r * np.cos(thetas))
        y = np.append(y, r * np.sin(thetas))

    return np.column_stack((x, y))


def rotmat_2d(ang_deg):  # TODO: could lru_cache?
    """Return rotation matrix for rotation `ang_deg` in degrees.
    For left-multiplication of a column position vector.

    .. note::
       [`scipy.spatial.transform.Rotation`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html)
       can be used for 3-d rotations.
    """
    ang = np.deg2rad(ang_deg)
    c, s = np.cos(ang), np.sin(ang)
    R = np.array([[c, -s], [s, c]])
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


@_add_to_vortons(
    short="Create polygonal arrangement of `Vortons` using `vertices_regular_polygon`."
)
@_add_snippets(snippets=dict(returns=_points_returns))
def vertices_regular_polygon(n, *, r_c=1):
    """Regular polygon vertices.

    Parameters
    ----------
    n : int
        Polygon order (number of sides/vertices).
    r_c : float, int
        Radius $r_c$ of the inscribing circle.

    Returns
    -------
    %(returns)s
    """
    # initial vertex
    vert0 = np.r_[0, r_c]

    # rotation matrix -- left-multiplies a column position vector to give rotated position
    rotmat = rotmat_2d(360 / n)

    verts = np.full((n, 2), vert0, dtype=float)
    # successive rotations
    for i in range(1, n):
        verts[i, :] = rotate_2d(verts[i - 1, :], rotmat=rotmat)

    return verts


@_add_to_vortons(
    short="Create isosceles triangle arrangement of `Vortons` using `vertices_isos_triangle`."
)
@_add_snippets(snippets=dict(returns=_points_returns))
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

    Returns
    -------
    %(returns)s
    """
    if (theta_deg is not None and Lambda is not None) or (theta_deg is None and Lambda is None):
        raise Exception("Specify either `theta_deg` or `Lambda` (not both).")

    if Lambda:
        assert Lambda > 0 and Lambda <= 1
        theta_deg = 180 / (Lambda**2 + 2)

    theta = np.deg2rad(theta_deg)

    xb = 1.5 / np.tan(theta)  # one half of x base

    xi = [-xb, 0, xb]
    yi = [-0.5, 1, -0.5]

    Lambda = np.sqrt((180 - 2 * theta_deg) / float(theta_deg))  # Marcelo eqns 17--19

    return np.column_stack((xi, yi))


@_add_to_vortons(short="Create asterisk arrangement of `Vortons` using `points_asterisk`.")
@_add_snippets(snippets=dict(returns=_points_returns))
def points_asterisk(n_limbs=5, n_per_limb=3, *, rmax=1):
    """Asterisk with `n_limbs` number of limbs and `n_per_limb` points per limb.

    Parameters
    ----------
    n_limbs : int
        For example, `5` to get a 5-pointed asterisk
    n_per_limb : int
        Number of evenly-spaced points in the limb, not including the center!
    rmax : float
        Limb length (maximum radius for on-limb points)

    Returns
    -------
    %(returns)s
    """
    assert n_limbs >= 1 and n_per_limb >= 1, "both n's must be >= 1"
    # Center point
    x = [0]
    y = [0]

    # Limbs (could use regular_polygon for this really)
    theta = np.deg2rad(np.linspace(90, 360 + 90, n_limbs + 1))[:-1]
    rs = np.linspace(0, rmax, n_per_limb + 1)[1:]
    for r in rs:
        x.extend((r * np.cos(theta)).tolist())
        y.extend((r * np.sin(theta)).tolist())

    return np.column_stack((x, y))


def _maybe_fill_G(G, n):
    if G is None:  # this first part maybe shouldn't be here? or kwarg for default G val?
        G = 1.0
    G = np.asarray(G)
    if G.size == 1:  # single G provided, or using the default
        G = np.full((n,), G)  # TODO: could also the constructor to accept single G
    if G.size != n:
        raise ValueError(f"`G` must have size `n` or 1, but is {G.size!r}")

    return G
