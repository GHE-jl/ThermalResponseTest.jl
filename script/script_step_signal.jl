"""
Segment a noisy variable-power signal into a small number of constant steps with `step_signal`
(k-means clustering). This helps turn a measured, fluctuating heating power into the piecewise-
constant states used by a non-stationary interpretation of a thermal response test.

Run from the package root:
    julia --project=script script/script_step_signal.jl
"""

using ThermalResponseTest
using CairoMakie

# QVar.csv holds a measured, fluctuating heating power [W] (no header: time, power).
datafile = joinpath(@__DIR__, "..", "trt_data", "QVar.csv")
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

mkpath(joinpath(@__DIR__, "figures"))
outfile = joinpath(@__DIR__, "figures", "step_signal.png")
save(outfile, fig)
println("Saved figure to $(outfile)")
# display(fig)   # uncomment for interactive use
