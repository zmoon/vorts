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
# Note that on Binder we need to install `gfortran`. In the future, we might use a repo `postBuild` file for that. For now, uncomment the below to install with `conda`:

# %%
# # !conda install -c conda-forge gfortran_linux-64 --yes
# # !mkdir -p ~/.local/bin
# # !ln -sf /srv/conda/envs/notebook/bin/x86_64-conda-linux-gnu-gfortran.bin ~/.local/bin/gfortran

# %% [markdown]
# ## Create case

# %% create
vs = vorts.Vortons([1, 1], [0, 0], [-0.5, 1])  # <-> Lambda=0 (no longer a triangle)
vs.plot()

ts = vorts.Tracers.randu(100)
ts.plot()

m = vorts.Model_f(
    vs, ts,
    dt=0.05, nt=2e5,  # 0.005, 2e6 to get better Poincare, but only use this with `write_ps` only to avoid writing large tracer files!
    int_scheme_name='RK4',
    write_vortons=True,  # default `True`
    write_tracers=True,  # default `False`
    write_ps=False,  # default `False`
)


# %% [markdown]
# ## Run and plot

# %% run
m.run()


# %% plot results
m.plot()
m.plot("tracers")


# %% [markdown]
# Since both have $x=0$, either one works to make a nice Poincare map.

# %% Poincare?
m.plot("poincare")

m.plot("poincare", iv_ref=1)  # this is the one with initial position (0, 1)
vorts.plot.remove_frame(keep_title=False)  # doesn't show up with `%matplotlib notebook`
