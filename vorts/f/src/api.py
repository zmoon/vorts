import ctypes as ct
import numpy as np

lib = np.ctypeslib.load_library("libvorts", ".")
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
x = np.array([-0.5, 0, 0.5], order="F")
y = np.array([-0.3, 0.8, 0.3], order="F")
G = np.full((nv,), 1., order="F")

dt = 0.1
nt = 5
imethod = 2

nout = nt + 1
xout = np.full((nv, nout,), 0., order="F")
yout = np.full((nv, nout,), 0., order="F")


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

print("From Python: output")
print(xout)
print(yout)

