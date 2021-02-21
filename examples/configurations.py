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
# # Built-in configurations

# %%
import matplotlib.pyplot as plt
from ipywidgets import interact

import sys; sys.path.append("../"); import vorts

# %matplotlib widget

# %% [markdown]
# ## Tracer

# %% [markdown]
# ### Spiral

# %%
fig1 = plt.figure()

def plot_spiral(n=100, revs=2, connect=False):
    fig = plt.figure(fig1.number); fig.clf(); ax = fig.add_subplot()
    vorts.Tracers.spiral(n=n, revs=revs).plot(ax=ax, connect=connect)

interact(plot_spiral, n=(1, 500), revs=(1, 10))

# %% [markdown]
# ### Uniform random

# %%
vorts.Tracers.randu(200).plot()

# %% [markdown]
# ### Gaussian random

# %%
vorts.Tracers.randn(400, sig_y=0.5).plot()

# %% [markdown]
# ### Grid

# %%
vorts.Tracers.grid(10, 20).plot()  # x first, then y

# %% [markdown]
# ### Concentric circles

# %%
vorts.Tracers.circles().plot()

# %% [markdown]
# ## Vortons

# %% [markdown]
# ### Isosceles triangle

# %%
vorts.Vortons.isos_triangle(theta_deg=72).plot()

# %% [markdown]
# ### Regular polygon

# %%
fig2 = plt.figure()

def plot_spiral(n=6):
    fig = plt.figure(fig2.number); fig.clf(); ax = fig.add_subplot()
    vorts.Vortons.regular_polygon(n=n).plot(ax=ax)

interact(plot_spiral, n=(3, 20))
