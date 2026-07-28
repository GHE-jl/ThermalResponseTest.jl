# Segment a noisy variable-power signal into a small number of constant steps with `step_signal`
# (k-means clustering). This helps turn a measured, fluctuating heating power into the piecewise-
# constant states used by a non-stationary interpretation of a thermal response test.

import Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate() # Once per computer

using ThermalResponseTest
using CairoMakie

# Q_Var.csv holds a measured, fluctuating heating power [W] (no header: time, power).
datafile = joinpath(@__DIR__, "..", "data", "Q_Var.csv")
Q = [parse(Float64, split(line, ',')[2]) for line in readlines(datafile)]

nsteps = 100
Qstep = step_signal(Q, nsteps)
println("Segmented $(length(Q)) samples into $(nsteps) constant steps " *
        "($(length(unique(Qstep))) distinct levels)")

fig = Figure(size = (900, 500))
ax = Axis(fig[1, 1], xlabel = "Sample", ylabel = "Power (W)",
    title = "Step segmentation of a variable heating power")
lines!(ax, 1:length(Q), Q, color = (:blue, 0.6), label = "Measured")
lines!(ax, 1:length(Qstep), Qstep, color = :red, linewidth = 2, label = "$(nsteps)-step signal")
axislegend(ax, position = :rb)
display(fig)
