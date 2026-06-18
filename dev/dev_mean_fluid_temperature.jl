"""
Script to compute the mean fluid temperature in the thermal response test.
"""

using CairoMakie, DataFrames

includet("../src/ThermalResponseTest.jl")
using .ThermalResponseTest

# Load data
file_path = "trt_data/DataCL_TRT.csv"
# file_path = "trt_data/DataCL.csv"
df = load_trt_data(file_path)

# Decompose into heating and cooling phases
dataset = decompose_trt(df)

# Compute mean fluid temperature using different methods
T_arithmetic = mean_fluid_temperature(df.T_in, df.T_out, :arithmetic)
T_logarithmic = mean_fluid_temperature(df.T_in, df.T_out, :logarithmic)
T_geometric = mean_fluid_temperature(df.T_in, df.T_out, :geometric)
T_harmonic = mean_fluid_temperature(df.T_in, df.T_out, :harmonic)
p = range(-3, stop=2, length=20)
T_pLinear = zeros(length(df.T_in), 20)
for i in 1:20
    T_pLinear[:, i] = mean_fluid_temperature(df.T_in, df.T_out, p[i])
end

# Plot Heating Phase
fig = Figure()
ax1 = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Temperature (°C)")
lines!(ax1, df.t, df.T_in, color=:blue, linewidth=2, label="T_in")
lines!(ax1, df.t, df.T_out, color=:red, linewidth=2, label="T_out")
lines!(ax1, df.t, df.T_mean, color=:green, linewidth=2, label="Mean")
axislegend(ax1)

ax2 = Axis(fig[1, 2], xlabel="Time (s)", ylabel="Temperature (°C)")
lines!(ax2, df.t, T_arithmetic, color=:cyan, linewidth=2.5, label="arithmetic")
lines!(ax2, df.t, T_logarithmic, color=:orange, linewidth=2, label="logarithmic")
lines!(ax2, df.t, T_geometric, color=:magenta, linewidth=1.5, label="geometric")
lines!(ax2, df.t, T_harmonic, color=:brown, linewidth=1, label="harmonic")
axislegend(ax2)

ax3 = Axis(fig[2, 1:2], xlabel="Time (s)", ylabel="Temperature (°C)")
for i in 1:20
    lines!(ax3, df.t, T_pLinear[:, i], color=:gray, linewidth=0.5, label="p=$(round(p[i],
        digits=2))")
end
axislegend(ax3)

display(fig)