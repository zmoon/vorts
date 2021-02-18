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

vs.plot()

ts = vorts.Tracers.randu(100)

ts.plot()

m = vorts.Model_f(
    vs, ts,
    dt=0.05, nt=2e5,  # 0.005, 2e6 to get better Poincare, but only use this with `write_ps` only to avoid writing large tracer files!
    int_scheme_name='RK4',
    write_vortons=True,  # default `True`
    write_tracers=True,  # default `False`
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
vorts.plot.frame_only()

m.plot("poincare", iv_ref=1, xtol=1e-3)
vorts.plot.remove_frame(keep_title=False)
