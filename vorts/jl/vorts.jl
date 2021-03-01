"""
Solving the vorts problem using [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
"""
# Note: pre-compiling DifferentialEquations takes *long* time (for Plots too!)
using DifferentialEquations
using Plots


"Pre-factor used for both the ``x`` and ``y`` tend calculations"
const TEND_PRE = 1 / 2pi


"Parameters used in the position tendency calculation"
struct TendParams
  G :: Vector{Float64}
  # n_vortons :: Integer
end


# Note: `@benchmark` estimates only 2 allocs, but should try to improve speed / Julian-ify
"""
    tend!

Update the position time tendency `dr`, with `dr` and `r` passed as arrays of size `(2, n)`.
That is, each col is the ``xy`` position vector for a given vorton (or tracer).
"""
function tend!(dr, r, p::TendParams, t)
  G = p.G
  n = length(G)

  # Compute tendencies for each vorton
  for i ∈ 1:n

    # Sum all contributions to the tendencies for vorton `i`
    dr[1, i] = 0  # dxdt
    dr[2, i] = 0  # dydt
    for j ∈ 1:n
      Gj = G[j]
      (i == j || Gj == 0) && continue  # skip if same vorton or tracer

      dx = r[1, i] - r[1, j]
      dy = r[2, i] - r[2, j]
      lij_sqd = dx^2 + dy^2

      dr[1, i] += -TEND_PRE*Gj * dy/lij_sqd
      dr[2, i] += TEND_PRE*Gj * dx/lij_sqd
    end
  end
end

# TODO: complex number formulation + benchmark


"""
    integrate(r₀, G, dt, nt; int_scheme_name)

Integrate and return the solution object.

# Arguments
* `r₀`: initial positions, where each row is an ``xy`` coordinate pair
* `G`: vector of ``\\Gamma`` values for each position
* `dt`: time step that we want the output to have
* `nt`: number of time steps taken from ``t=0``
* `int_scheme_name`: string of a solver name, e.g., from the
  [ODE solvers list](https://diffeq.sciml.ai/stable/solvers/ode_solve/)
"""
function integrate(r₀, G, dt, nt; int_scheme_name="Tsit5")
  # Inputs
  r₀ = permutedims(r₀)  # for Julia, we want coords as cols
  @assert length(G) == size(r₀, 2)
  p = TendParams(G)
  tspan = (0.0, nt*dt)
  solver = getfield(Main, Symbol(int_scheme_name))

  # Use DifferentialEquations
  prob = ODEProblem(tend!, r₀, tspan, p)
  sol = solve(prob, solver(); saveat=0:dt:tspan[2])
  # Notes:
  # * with default settings, ends up using Tsit5 (`all(sol.alg_choice .== 1)` is true)
  # * might want to enable dense for Poincare purposes in the future (not compatible with `saveat`)
  # * could consider using `dtmax=dt`

  # Note that `sol.u` is the solution, as a vec of arrays like r₀ (each r(t) is its own array)
  return sol
end


"Plot the solution on Julia side for testing."
function plot_sol(sol)
  xymax = maximum(abs.(hcat(sol.u...)))
  lims = (-xymax, xymax)
  p_plot = (aspect_ratio=:equal, xlabel=raw"$x$", ylabel=raw"$y$", xlims=lims, ylims=lims)
  plot()
  for i ∈ 1:size(sol.u[1], 2)
    # Column-major, so (1, 3) is first row, (1, 2) is first col
    plot!(sol, vars=(2i-1, 2i), label="$i")
  end
  plot!(;p_plot...)
  gui()
end


# try it out
function main()
  sol = integrate([0 1 ; 0 0.01; 0 -1], ones(3), 0.1, 500)
  plot_sol(sol)
end


if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
