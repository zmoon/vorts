import ctypes as ct
from pathlib import Path

import numpy as np

_here = Path(__file__).parent
_builddir = _here / "_build"
_libbasename = "libvortsf"


try:
    _lib = np.ctypeslib.load_library(_libbasename, _builddir)
except OSError:
    import platform
    import sys

    if platform.system() == "Windows" and sys.version_info >= (3, 8):
        _lib = ct.CDLL((_builddir / f"{_libbasename}.dll").as_posix(), winmode=0)
        # `winmode=0` is supposedly a default
        #   https://docs.python.org/3/library/ctypes.html#ctypes.CDLL
        # but it doesn't work without specifying it.
        # Related issue: https://bugs.python.org/issue42114
        # Looks like it really defaults to `nt._LOAD_LIBRARY_SEARCH_DEFAULT_DIRS`,
        # which is nonzero on my Windows machine.
        # From https://docs.python.org/3/library/ctypes.html#loading-shared-libraries
        # `winmode` was added in Python 3.8.

    else:
        raise


_run = _lib.asdf
_run.argtypes = [
    ct.c_int,
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
    ct.c_double,
    ct.c_int,
    ct.c_int,
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags="F, W"),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags="F, W"),
]

# Must match order in `int_scheme_from_index`!
_methods = ["FT", "RK4"]


def run(G, x, y, *, dt, nt, method="RK4"):

    nv = G.size
    assert nv == x.size == y.size

    imethod = _methods.index(method) + 1

    # Allocate outputs
    xout = np.full((nv, nt + 1), 0.0, order="F")
    yout = np.full((nv, nt + 1), 0.0, order="F")

    _run(
        nv,
        x,
        y,
        G,
        dt,
        nt,
        imethod,
        xout,
        yout,
    )

    return xout, yout


if __name__ == "__main__":

    x = np.array([-0.5, 0, 0.5])
    y = np.array([-0.3, 0.8, 0.3])
    G = np.ones_like(x)
    dt = 0.1
    nt = 5

    xout, yout = run(G, x, y, dt=dt, nt=nt)

    print("From Python: output")
    print(xout)
    print(yout)
