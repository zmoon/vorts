#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 13:32:15 2018

@author: zmoon
"""

import matplotlib.pyplot as plt
import numpy as np

import vortspy as vpy

plt.close('all')



#%% create

#G = np.ones(2)
#G = [1, 5]
#xi = [0, 0]
#yi = [-1, 0.5]
#yi = [-1, 1]

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

m = vpy.model_f(G, xi, yi,
                xit=xit, yit=yit,
                dt=0.005, nt=10000000,
                int_scheme_name='RK4')



#%% run

m.run()


#%%

colors = plt.cm.Dark2(np.linspace(0, 1, 8)[2:2+len(G)+1])


#%% plot tracks

#fig1, f1a1 = plt.subplots(figsize=(4.5, 4), num='vorton orbits')
#
#
#for i, v in enumerate(m.vortons):
#
#    x = v.xhist
#    y = v.yhist
#
#    f1a1.plot(x, y, color=colors[i], lw=0.5, alpha=0.5)
#    f1a1.plot(x[0], y[0], 'o', color=colors[i])
#
#    f1a1.set_xlabel('x')
#    f1a1.set_ylabel('y')
#
#    f1a1.axis('equal')


#%% plot tracers
# need to find a way to cut this down if there are too many to see what is happening

#fig2, f2a1 = plt.subplots(figsize=(4.5, 4), num='tracer orbits')
#
#offset = 300
#
#for i, t in enumerate(m.tracers):
#
#    x = t.xhist
#    y = t.yhist
#
#    xdisp = x[offset:]
#    ydisp = y[offset:]
#
#    f2a1.plot(xdisp, ydisp, '.', color='0.35', ms=1, alpha=0.1)
#
#    f2a1.set_xlabel('x')
#    f2a1.set_ylabel('y')
#
#    f2a1.axis('equal')


#%% Poincare sections

#> find indices of return

t_crosses = []
psx = []
psy = []

for i, v in enumerate(m.vortons):

    x = v.xhist
    y = v.yhist

    #> to make the Poincare section
    #  we want to mark when the vorton returns to original position

    # using loop
#    t_cross = []
#    psxi = []
#    psyi = []
#    for j in range(1, x.size):
#        if x[j] < 0 and x[j-1] > 0 and y[j] > 0:  # really only works for the one at the top (0, 1)
#            psyi.append(y[j])
#            psxi.append(x[j])
#            t_cross.append(j)

    # using fancy indexing
    inds = np.arange(0, x.size)
    j   = inds[1:]
    jm1 = inds[:-1]
    if i == 1: # top one
        ret = (x[j] < xi[i]) & (x[j-1] > xi[i]) & (y[j] > 0)  # tells us where return is true, coming from top
    else:  # bottom ones
        ret = (x[j] > xi[i]) & (x[j-1] < xi[i]) & (y[j] < 0)
    ret = np.insert(ret, 0, True)

    psxi = x[ret]
    psyi = y[ret]
    t_cross = inds[ret]

    t_crosses.append(t_cross)


#> record tracer positions at indices of return to construct Poincare Section

psty = []
pstx = []

i_tcross = 1

for i, t in enumerate(m.tracers):

    x = t.xhist
    y = t.yhist

    # combine Poincare Section results for the individual tracers
    pstx += list(x[t_crosses[i_tcross]])
    psty += list(y[t_crosses[i_tcross]])


#> plot

fig3, f3a1 = plt.subplots(figsize=(8, 7), num='ps')

f3a1.plot(pstx, psty, '.', color='0.35', ms=1, alpha=0.1)

for i, v in enumerate(m.vortons):
    x0 = v.xhist[0]
    y0 = v.yhist[0]
    f3a1.plot(x0, y0, 'o', color=colors[i])

    f3a1.set_xlabel('x')
    f3a1.set_ylabel('y')




f3a1.axis('equal')

fig3.savefig('ps_theta{:d}deg.png'.format(theta_deg), dpi=600,
             transparent=True, bbox_inches='tight', pad_inches=0.05)
