
## Setup

0. Install Julia

1. Follow the instructions [here](https://julialang.github.io/Pkg.jl/v1/environments/#Using-someone-else's-project)
   to create an environment based on the [`Project.toml`](./Project.toml).

2. With [pyjulia](https://github.com/JuliaPy/pyjulia) installed (included in `vorts`'s deps), you can run
   ```python
   import julia; julia.install()
   ```
   to check that your Julia/Python installations are ready for pyjulia.
   This will install [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) in  your base Julia environment if you don't already have it.
