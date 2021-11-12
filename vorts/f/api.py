import ctypes as ct
import numpy as np

_lib = np.ctypeslib.load_library("libvortsf", "./_build")

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
    xout = np.full((nv, nt+1,), 0., order="F")
    yout = np.full((nv, nt+1,), 0., order="F")

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
