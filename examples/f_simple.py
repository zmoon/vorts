# -*- coding: utf-8 -*-
"""
A test run with tracers.
"""
import sys
sys.path.append("../")

import matplotlib.pyplot as plt
# import numpy as np

import vorts

plt.close("all")


# %% create

# vs = vorts.Vortons.isos_triangle(G=1, Lambda=1)  # <-> equilateral triangle (theta_deg=60)

vs = vorts.Vortons([1, 1], [0, 0], [-0.5, 1])  # <-> Lambda=0 (no longer a triangle)

vs.add_tracers(100)

vs.plot()

m = vorts.Model_f(
    vs,
    dt=0.005, nt=2e6,
    int_scheme_name='RK4',
    write_vortons=True,  # default `True`
    write_tracers=False,  # default `False`
    write_ps=False,  # default `False`
)


# %% run

m.run()


# %% plot results

m.plot()
m.plot("tracers")


# %% Poincare?

m.plot("poincare")
m.plot("poincare", iv_ref=1)  # this is the one with initial position (0, 1)
m.plot("poincare", iv_ref=1, xtol=1e-3)
m.plot("poincare", iv_ref=1, xtol=1e-6)
