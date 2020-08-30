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

#G = np.ones(2)
# G = [1, 5]
# xi = [0, 0]
# yi = [-1, 0.5]

G = np.ones(3)
xi = [-0.3, 0, 0.3]
yi = [-0.5, 1, 0]

# G = [1, 1, 1, 0, 0]
# xi = [-0.3, 0, 0.3, 0.5, 0.1]
# yi = [-0.5, 1, 0, 0, 0]


vs = vorts.Vortons(G, xi, yi)
vs.add_tracers(100)

vs.plot()  # plot initial state

m = vorts.Model_py(
    vs,
    dt=0.1, nt=200,
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



# %% plot run results

m.plot()  # vortons by default

m.plot("tracers")
