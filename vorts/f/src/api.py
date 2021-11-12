import ctypes as ct
import numpy as np

lib = ct.CDLL("./vorts.so")
f = lib.asdf

f.argtypes = [
    ct.c_int,
    ct.POINTER(ct.c_double),    
    ct.POINTER(ct.c_double),    
    ct.POINTER(ct.c_double),
    ct.c_double,
    ct.c_int,
    ct.c_int,
    ct.POINTER(ct.c_double),
    ct.POINTER(ct.c_double),
]

nv = 3
x = np.full((nv,), 1., order="F")
y = np.full((nv,), 2., order="F")
G = np.full((nv,), 3., order="F")

dt = 0.1
nt = 3
imethod = 0

nout = nt + 1
xout = np.full((nout,), 0., order="F")
yout = np.full((nout,), 0., order="F")


f(
    ct.c_int(nv),
    x.ctypes.data_as(ct.POINTER(ct.c_double)), 
    y.ctypes.data_as(ct.POINTER(ct.c_double)),
    G.ctypes.data_as(ct.POINTER(ct.c_double)),
    ct.c_double(dt),
    ct.c_int(nt),
    ct.c_int(imethod),
    xout.ctypes.data_as(ct.POINTER(ct.c_double)), 
    yout.ctypes.data_as(ct.POINTER(ct.c_double)), 
)
