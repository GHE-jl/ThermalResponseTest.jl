"""
Script that tests methods to compute the time derivative of a temperature signal obtained from a
thermal response test.
"""

using CairoMakie

includet("../src/ThermalResponseTest.jl")
using .ThermalResponseTest

# Simple example with the infinite line source
t = 0:60:365 * 24 * 3600
rb = 0.075
k = 2.5
Cs = 2e6
T = ils(t, rb, k, Cs)

# Compute the time derivative using a centered finite difference method
dT = centered_finite_difference(t, T)

# Plot results
fig = Figure()
ax1 = Axis(fig[1, 1], title="Temperature vs Time", xlabel="Time (h)", ylabel="Temperature (°C)",
    xscale=log10, yscale=log10)
lines!(ax1, dT, color=:black, linewidth=2, label="Temperature")
axislegend(ax1, position=:rb)
display(fig)