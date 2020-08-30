# -*- coding: utf-8 -*-
"""
A test run with tracers.
"""
import sys
sys.path.append("../")

import matplotlib.pyplot as plt
import numpy as np

import vorts


# %% create

G = np.ones(3) * 1

theta_deg = 60  # 72 gives a Lambda near the critical 1/sqrt(2) (I think)
theta = np.deg2rad(theta_deg)  # angle between base and connections to the top point at (0,1)

xb = 1.5/np.tan(theta)  # one half of x base

xi = [-xb,  0,  xb]
yi = [-0.5, 1, -0.5]

Lambda = np.sqrt( (180-2*theta_deg) / float(theta_deg) )  # Marcelo eqns 17--19


n = 12  # sqrt of num tracers
d = 3.0  # displacement off center to select initial tracer coords from
xit = np.random.uniform(-d, d, (n, n))
yit = np.random.uniform(-d, d, (n, n))

m = vorts.Model_f(
    G, xi, yi,
    xit=xit, yit=yit,
    dt=0.005, nt=2e5,
    int_scheme_name='RK4',
    write_vortons=True,  # default `True`
    write_tracers=False,  # default `False`
    write_ps=True,  # default `False`
)


# %% run

m.run()


# %% plot

fig, ax = plt.subplots()

# plot tracer positions at Poincare section times
xt, yt = m.ps.T
ax.plot(xt, yt, ".", c="0.5", ms=4, alpha=0.5)

# plot vorton tracks
colors = plt.cm.Dark2(np.linspace(0, 1, 8)[2:2+len(G)+1])
for i, v in enumerate(m.vortons):
    x = v.xhist
    y = v.yhist
    ax.plot(x, y, color=colors[i], lw=0.5, alpha=0.5)
    ax.plot(x[0], y[0], 'o', color=colors[i])

# ax.axis("equal")
ax.set_aspect("equal", "box")

# plt.show()
