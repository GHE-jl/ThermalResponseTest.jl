"""
Load a thermal response test data file, decompose it into heating and recovery phases, and plot the
inlet / outlet / mean fluid temperatures and the heating power.

Run from the package root:
    julia --project=script script/script_load_trt.jl
"""

using ThermalResponseTest
using DataFrames
using CairoMakie

# Plot helper (plotting lives in the scripts, not the library).
function plot_trt(df::DataFrame; title = "")
    fig = Figure(size = (900, 600))
    ax1 = Axis(fig[1, 1], title = title, xlabel = "Time (s)", ylabel = "Temperature (°C)")
    lines!(ax1, df.elapsed_time, df.T_in,   color = :red,   linewidth = 2, label = "T_in")
    lines!(ax1, df.elapsed_time, df.T_out,  color = :blue,  linewidth = 2, label = "T_out")
    lines!(ax1, df.elapsed_time, df.T_mean, color = :green, linewidth = 2, label = "T_mean")
    axislegend(ax1, position = :rb)
    ax2 = Axis(fig[2, 1], xlabel = "Time (s)", ylabel = "Power (W)")
    lines!(ax2, df.elapsed_time, df.power, color = :orange, linewidth = 2, label = "Power")
    axislegend(ax2, position = :rt)
    return fig
end

datafile = joinpath(@__DIR__, "..", "trt_data", "DataCL_TRT.csv")
trt = load_trt_data(datafile)
dataset = decompose_trt(trt)

println("Loaded $(nrow(trt)) samples: $(nrow(dataset.heating)) heating, " *
        "$(nrow(dataset.cooling)) recovery")

fig = plot_trt(trt; title = "TRT — full data")

mkpath(joinpath(@__DIR__, "figures"))
outfile = joinpath(@__DIR__, "figures", "load_trt.png")
save(outfile, fig)
println("Saved figure to $(outfile)")
# display(fig)   # uncomment for interactive use
