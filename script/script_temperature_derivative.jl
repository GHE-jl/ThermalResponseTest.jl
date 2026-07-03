"""
Compute the time derivative of a temperature signal with a centered finite difference. The
derivative is the basis of the CFOA-Ṫ interpretation methods (Pasquier 2018); here it is
illustrated on a synthetic infinite-line-source signal where the analytical derivative is known.

Run from the package root:
    julia --project=script script/script_temperature_derivative.jl
"""

using ThermalResponseTest
using GroundHeatExchanger          # ils (re-exported from GroundResponse.jl)
using CairoMakie

# Synthetic ILS temperature signal
t  = collect(60.0:60.0:365 * 24 * 3600)   # 1 min → 1 year, uniform 60 s
rb = 0.075
k  = 2.5
Cs = 2.0e6
q  = 50.0                                   # heat injection per length [W/m]
T  = q .* ils(t, rb, k, Cs)                 # borehole-wall temperature rise [°C]

# Numerical derivative vs the analytical ILS derivative  dT/dt = q/(4πk) · 1/t
dT_num   = centered_finite_difference(t, T)
dT_exact = q ./ (4π * k) ./ t

fig = Figure(size = (900, 500))
ax = Axis(fig[1, 1], title = "Time derivative of an ILS temperature signal",
    xlabel = "Time (s)", ylabel = "dT/dt (°C/s)", xscale = log10, yscale = log10)
lines!(ax, t, dT_exact, color = :black, linewidth = 3, label = "Analytical  q/(4πk·t)")
lines!(ax, t, dT_num,   color = :red,   linewidth = 2, linestyle = :dash,
    label = "Centered finite difference")
axislegend(ax, position = :rt)

mkpath(joinpath(@__DIR__, "figures"))
outfile = joinpath(@__DIR__, "figures", "temperature_derivative.png")
save(outfile, fig)
println("Saved figure to $(outfile)")
# display(fig)   # uncomment for interactive use
