# -*- coding: utf-8 -*-
"""
Kind of testbed for making sure the options work...
"""
import sys
sys.path.append("../")

import matplotlib.pyplot as plt
import numpy as np

import vorts



#%% create

#G = np.ones(2)
# G = [1, 5]
# xi = [0, 0]
# yi = [-1, 0.5]

# G = np.ones(3)
# xi = [-0.3, 0, 0.3]
# yi = [-0.5, 1, 0]

G = [1, 1, 1, 0, 0]
xi = [-0.3, 0, 0.3, 0.5, 0.1]
yi = [-0.5, 1, 0, 0, 0]


m = vorts.model_py(
    G, xi, yi,
    dt=0.1, nt=2000,
    # int_scheme_name="not-a-scheme",
    # int_scheme_name='FT_2',
    # int_scheme_name="RK4_2",
    # int_scheme_name="RK4_3",
    # int_scheme_name='scipy_DOP853',
    # adapt_tstep=False,
    # adapt_tstep=True,
)



#%% run

m.run()


#%% plot tracks

colors = plt.cm.Dark2(np.linspace(0, 1, 8)[2:2+len(G)+1])


fig1, f1a1 = plt.subplots(figsize=(4.5, 4))

for i, v in enumerate(m.vortons):

    x = v.xhist
    y = v.yhist

    f1a1.plot(x, y, color=colors[i], lw=0.5, alpha=0.5)
    f1a1.plot(x[0], y[0], 'o', color=colors[i])

    f1a1.set_xlabel('x')
    f1a1.set_ylabel('y')

    f1a1.axis('equal')
