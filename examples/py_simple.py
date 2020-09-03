# -*- coding: utf-8 -*-
"""
Kind of testbed for making sure the options work...
"""
import sys
sys.path.append("../")

import matplotlib.pyplot as plt
import numpy as np

import vorts

plt.close("all")


#%% create

G = [1, 5]
xi = [0, 0]
yi = [-1, 0.5]

# G = np.ones(3)
# xi = [-0.3, 0, 0.3]
# yi = [-0.5, 1, 0]


vs = vorts.Vortons(G, xi, yi)
vs.plot()  # plot initial state

ts = vorts.Tracers.spiral(100, c=vs.cm())
ts.plot()

m = vorts.Model_py(
    vs,
    ts,
    dt=0.1, nt=2000,
    # int_scheme_name="not-a-scheme",  # raises ValueError
    # int_scheme_name="FT",
    # int_scheme_name="FT_1b1",
    int_scheme_name="RK4",
    # int_scheme_name="RK4_1b1",
    # int_scheme_name="scipy_RK45",
    # int_scheme_name='scipy_DOP853',
    # adapt_tstep=False,
    # adapt_tstep=True,
)



#%% run

m.run()



# %% plot run results

m.plot()  # vortons by default

m.plot("tracers")
