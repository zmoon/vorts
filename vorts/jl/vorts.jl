"""
Solving the vorts problem using [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
"""
# Note: pre-compiling DifferentialEquations takes *long* time (for Plots too!)
using DifferentialEquations: ODEProblem, solve, Vern7
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


# Set up problem inputs
r₀ = [
  0 3;
  0 -0.1;
  0 -3;
]'  # coords as rows -> coords as cols
G = [1, 5, 1]
@assert length(G) == size(r₀, 2)
p = TendParams(G)
tspan = (0.0, 1e4)

# Integrate
prob = ODEProblem(tend!, r₀, tspan, p)
sol = solve(prob)  # ends up using Tsit5 by default (`all(sol.alg_choice .== 1)` is true)
# sol = solve(prob, Vern7())

# Plot
lims = (-5, 5)
p_plot = (aspect_ratio=:equal, xlabel=raw"$x$", ylabel=raw"$y$", xlims=lims, ylims=lims)
plot()
for i ∈ 1:length(p.G)
  # Column-major, so (1, 3) is first row, (1, 2) is first col
  plot!(sol, vars=(2i-1, 2i), label="$i")
end
plot!(;p_plot...)
gui()
