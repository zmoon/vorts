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

import vorts

# %matplotlib widget

# %% [markdown]
# ## Tracer

# %% [markdown]
# ### Spiral

# %%
fig1 = plt.figure()


def plot_spiral(n=100, revs=3, rmin=0, kind="Archimedean", spacing="linear", connect=True):
    fig = plt.figure(fig1.number); fig.clf(); ax = fig.add_subplot()  # noqa: E702
    vorts.Tracers.spiral(n=n, rmin=rmin, revs=revs, kind=kind, spacing=spacing).plot(ax=ax, connect=connect)


interact(plot_spiral, n=(1, 500), revs=(1, 10), rmin=(0, 1, 0.1), kind=["Archimedean", "Fermat's"], spacing=["linear", "log", "inv-exp", "1/x"])

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


def plot_polygon(n=6):
    fig = plt.figure(fig2.number); fig.clf(); ax = fig.add_subplot()  # noqa: E702
    vorts.Vortons.regular_polygon(n=n).plot(ax=ax)


interact(plot_polygon, n=(3, 20))

# %% [markdown]
# ### *n*-pointed asterisk

# %%
fig3 = plt.figure()


def plot_asterisk(n_limbs=5, n_per_limb=3):
    fig = plt.figure(fig3.number); fig.clf(); ax = fig.add_subplot()  # noqa: E702
    vorts.Vortons.asterisk(n_limbs, n_per_limb).plot(ax=ax)


interact(plot_asterisk, n_limbs=(1, 12), n_per_limb=(0, 10))

# %% [markdown]
# ## Combining

# %%
(vorts.Tracers.circles() + vorts.Tracers.randn(50, sig_x=0.3, sig_y=0.3)).plot()

# %% [markdown]
# ## Transforms
#
# Note that these transforms also work for `Tracers`.

# %%
v0 = vorts.Vortons.regular_polygon(4)

# %%
(v0 + (1, 2)).plot()

# %%
(0.1*v0).plot()

# %%
v0.rotate(20).plot()
