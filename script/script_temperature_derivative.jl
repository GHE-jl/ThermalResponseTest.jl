# Compute the time derivative of a temperature signal in two ways: a centered finite difference,
# and the Bourdet et al. (1989) three-point formula (Beier 2020, Appendix C). Both are the basis of
# the CFOA-Ṫ interpretation methods (Pasquier 2018). Here they are illustrated on a synthetic
# infinite-line-source signal where the analytical derivative is known.

import Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate() # Once per computer

using ThermalResponseTest
using GroundHeatExchanger          # ils model
using CairoMakie
using Random

# Synthetic ILS temperature signal
t  = collect(60.0:60.0:365 * 24 * 3600)   # 1 min → 1 year, uniform 60 s
rb = 0.075
k  = 2.5
Cs = 2.0e6
q  = 50.0                                   # heat injection per length [W/m]
T  = q .* ils(t, rb, k, Cs)                 # borehole-wall temperature rise [°C]

# Numerical derivatives vs the analytical ILS derivative  dT/dt = q/(4πk) · 1/t
dT_centered = centered_finite_difference(t, T)
dT_bourdet  = bourdet_derivative(t, T)                # default δ = 0.3
dT_exact    = q ./ (4π * k) ./ t

fig = Figure()
ax = Axis(fig[1, 1], title = "Time derivative of an ILS temperature signal",
    xlabel = "Time (s)", ylabel = "dT/dt (°C/s)", xscale = log10, yscale = log10)
lines!(ax, t, dT_exact,    color = :black, linewidth = 3, label = "Analytical  q/(4πk·t)")
lines!(ax, t, dT_centered, color = :red,   linewidth = 2, linestyle = :dash,
    label = "Centered finite difference")
lines!(ax, t, dT_bourdet,  color = :blue,  linewidth = 2, linestyle = :dot,
    label = "Bourdet derivative (δ = 0.3)")
axislegend(ax, position = :rt)
display(fig)

# Effect of the smoothing parameter δ on the Bourdet derivative.
Random.seed!(1)
T_noisy = T .+ 0.02 .* randn(length(T))     # ±0.02 °C measurement noise
δ_values = [0.1, 0.3, 0.5, 1.0]

fig2 = Figure()

ax2 = Axis(fig2[1, 1], title = "Effect of δ on the Bourdet derivative (noisy signal)",
    xlabel = "Time (s)", ylabel = "dT/dt (°C/s)", xscale = log10, yscale = Makie.pseudolog10)
lines!(ax2, t, dT_exact, color = :black, linewidth = 3, label = "Analytical  q/(4πk·t)")
for δ in δ_values
    dT_δ = bourdet_derivative(t, T_noisy, δ)
    lines!(ax2, t, dT_δ, linewidth = 2, label = "δ = $δ")
end
axislegend(ax2, position = :rt)
display(fig2)
