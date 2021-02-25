# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.10.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Simple cases using the Python integrator(s)

# %%
import sys
import warnings
from timeit import default_timer as timer

import matplotlib.pyplot as plt
import numpy as np
from ipywidgets import interact

sys.path.append("../")
warnings.filterwarnings("ignore")  # TODO: warned about passing un-needed kwargs to solve_ivp

import vorts

# %matplotlib widget


# %% [markdown]
# ## General method

# %% [markdown]
# ### Create case and examine

# %%
G = [1, 5]
xi = [0, 0]
yi = [-1, 0.5]

fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(7, 3.5), sharex=True, sharey=True)

vs = vorts.Vortons(G, xi, yi)
vs.plot(ax=ax1, adjustable=None)  # plot initial state

ts = vorts.Tracers.spiral(6, c=vs.cm(), revs=1)
ts.plot(ax=ax2)

m = vorts.Model_py(
    vs, ts,
    dt=0.1, nt=20000,
    use_tqdm="notebook",
)

# %% [markdown]
# ### Run and plot

# %%
m.run()

# %% plot run results
m.plot()  # vortons by default

# %%
m.plot("tracers")

# %% [markdown]
# ## Interactive, examining impact of options
#
# In this case (the default, equilateral triangle), the vorton trajectories in the solution should trace out a perfect circle. We can see which methods give us this. Additionally, we estimate the run time to see which are faster.
#
# Note that `adapt_tstep` and `use_tqdm` are not used by the SciPy integrator options. Also, for the SciPy integrators, `dt` is used as the maximum step size.

# %%
fig1 = plt.figure(figsize=(5, 5))

def run_plot(dt=0.1, int_scheme="RK4", adapt_tstep=False):
    fig = plt.figure(fig1.number); fig.clf(); ax = fig.add_subplot()
    m = vorts.Model_py(dt=dt, nt=int(100/dt), int_scheme_name=int_scheme, adapt_tstep=adapt_tstep, use_tqdm=False)
    start = timer()
    m.run()
    finish = timer()
    m.plot(ax=ax)
    ax.set_title(f"Model time: {m.hist.t.isel(t=-1).values}\nWall time: {finish-start:.4g} s\n", loc="left")
    fig.tight_layout()

interact(run_plot, dt=[0.01, 0.1, 1.0, 2.0, 5.0], int_scheme=[n for n in vorts.Model_py._allowed_int_scheme_names if "1b1" not in n])

# %% [markdown]
# 👆 We can see that FT with adaptive time stepping leads to an accurate solution, but it takes much longer than other methods.
