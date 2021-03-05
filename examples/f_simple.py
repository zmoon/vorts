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
# # Fortran model

# %%
import vorts

# %matplotlib widget


# %% [markdown]
# ## Create case

# %% create
vs = vorts.Vortons.regular_polygon(3) + vorts.Vortons(-1, 0, 0)
vs.plot()

ts = vorts.Tracers.randu(100, dx=1.5, dy=1.5)
ts.plot()

m = vorts.Model_f(
    vs, ts,
    dt=0.05, nt=2e4,
    int_scheme_name='RK4',
    write_tracers=True,
)


# %% [markdown]
# ## Run and plot

# %% run
m.run()


# %% plot results
m.plot()  # vortons are stationary!
m.plot.tracers()


# %% [markdown]
# In this special completely stationary case (non-rotating), all points of the tracer trajectories can be used in the Poincare section.
#
# Note that in rotating cases, `dt=0.005, nt=2e6` would give a better Poincare section, but we would only want to use this with `write_ps` only to avoid writing large tracer files. And, on the usual Binder, a longer run than what we have set would use too much memory and fail.

# %%
m.plot.poincare(c=["teal", "tomato"])
