# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 15:18:29 2020

@author: zmoon
"""

import numpy as np

import vorts


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

m = vorts.model_f(
    G, xi, yi,
    xit=xit, yit=yit,
    dt=0.005, nt=1000,
    int_scheme_name='RK4',
    # write_ps=True,
)

m.run()

