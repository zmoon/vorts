# -*- coding: utf-8 -*-
"""
A test run with tracers.
"""
import sys
sys.path.append("../")

import matplotlib.pyplot as plt
import numpy as np

import vorts

plt.close("all")


# %% create

G = np.ones(3) * 1

theta_deg = 60  # 72 gives a Lambda near the critical 1/sqrt(2) (I think); 60 -> equi tri
theta = np.deg2rad(theta_deg)  # angle between base and connections to the top point at (0,1)

xb = 1.5/np.tan(theta)  # one half of x base

xi = [-xb,  0,  xb]
yi = [-0.5, 1, -0.5]

Lambda = np.sqrt( (180-2*theta_deg) / float(theta_deg) )  # Marcelo eqns 17--19

vs = vorts.Vortons(G, xi, yi)
vs.add_tracers(100)

vs.plot()

m = vorts.Model_f(
    vs,
    dt=0.005, nt=2e5,
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
