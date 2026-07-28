# Compare the mean-fluid-temperature averaging methods (p-linear mean of Marcotte & Pasquier 2008)
# on a thermal response test data set.

import Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate() # Once per computer

using ThermalResponseTest
using CairoMakie

datafile = joinpath(@__DIR__, "..", "data", "TRT_CL_2Phases.csv")
df = load_trt_data(datafile)

# Named methods
T_arithmetic  = mean_fluid_temperature(df.T_in, df.T_out, :arithmetic)
T_logarithmic = mean_fluid_temperature(df.T_in, df.T_out, :logarithmic)
T_geometric   = mean_fluid_temperature(df.T_in, df.T_out, :geometric)
T_harmonic    = mean_fluid_temperature(df.T_in, df.T_out, :harmonic)

# Continuum of p-linear exponents
p = range(-3, stop = 2, length = 20)
T_pLinear = hcat([mean_fluid_temperature(df.T_in, df.T_out, pᵢ) for pᵢ in p]...)

fig = Figure(size = (1200, 500))
ax1 = Axis(fig[1, 1], title = "Inlet / outlet / mean", xlabel = "Time (s)",
    ylabel = "Temperature (°C)")
lines!(ax1, df.elapsed_time, df.T_in,   color = :blue,  linewidth = 2, label = "T_in")
lines!(ax1, df.elapsed_time, df.T_out,  color = :red,   linewidth = 2, label = "T_out")
lines!(ax1, df.elapsed_time, df.T_mean, color = :green, linewidth = 2, label = "T_mean (p-linear)")
axislegend(ax1, position = :rb)

ax2 = Axis(fig[1, 2], title = "Averaging methods", xlabel = "Time (s)",
    ylabel = "Temperature (°C)")
lines!(ax2, df.elapsed_time, T_arithmetic,  linewidth = 2.5, label = "arithmetic")
lines!(ax2, df.elapsed_time, T_logarithmic, linewidth = 2,   label = "logarithmic")
lines!(ax2, df.elapsed_time, T_geometric,   linewidth = 1.5, label = "geometric")
lines!(ax2, df.elapsed_time, T_harmonic,    linewidth = 1,   label = "harmonic")
for i in eachindex(p)
    lines!(ax2, df.elapsed_time, T_pLinear[:, i], color = (:gray, 0.4), linewidth = 0.5)
end
axislegend(ax2, position = :rb)
display(fig)
