# Load thermal response test data files and decompose them into heating and recovery phases.
# This script also validates two `load_trt_data` behaviors:
#   1. A TRT logged with a pre-heating recirculation phase (`TRT_CL_3Phases.csv`) loads to the same
#      "t = 0 at the start of heating" DataFrame as one logged without it (`TRT_CL_2Phases.csv`),
#      since `load_trt_data` trims the recirculation phase and rebases time by default.
#   2. A TRT with non-uniform time steps (`TRT_CL_Num.csv`, a Comsol export with an adaptive-time-
#      stepping ramp at the start) is automatically interpolated onto a uniform grid.

import Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate() # Once per computer

using ThermalResponseTest
using CSV, DataFrames
using CairoMakie

# 1. Recirculation trim/rebase: with vs without a logged recirculation phase
file_2phases = joinpath(@__DIR__, "..", "data", "TRT_CL_2Phases.csv")   # heating + recovery
file_3phases = joinpath(@__DIR__, "..", "data", "TRT_CL_3Phases.csv")   # recirc + heat + recovery

trt_3_raw = load_trt_data(file_3phases; trim_recirculation = false)     # untrimmed, for comparison
trt_2 = load_trt_data(file_2phases)                                     # default: trim + rebase
trt_3 = load_trt_data(file_3phases)                                     # default: trim + rebase
d_2 = decompose_trt(trt_2)
d_3 = decompose_trt(trt_3)

println("3-phase file, untrimmed : $(nrow(trt_3_raw)) samples, ",
    "heating starts at t = $(trt_3_raw.elapsed_time[findfirst(>(100), trt_3_raw.power)]) s")
println("2-phase file, trimmed   : $(nrow(trt_2)) samples, heating starts at t = ",
    "$(trt_2.elapsed_time[findfirst(>(100), trt_2.power)]) s")
println("3-phase file, trimmed   : $(nrow(trt_3)) samples, heating starts at t = ",
    "$(trt_3.elapsed_time[findfirst(>(100), trt_3.power)]) s")
println("Heating duration        : 2-phase = $(d_2.heating.t_rel[end]) s, ",
    "3-phase = $(d_3.heating.t_rel[end]) s")

fig1 = Figure()

ax0 = Axis(fig1[1, 1], xlabel="Time (s)", ylabel="Temperature (°C)",
    title="3-phase file, untrimmed (trim_recirculation = false)")
lines!(ax0, trt_3_raw.elapsed_time, trt_3_raw.T_mean, color=:purple, linewidth=2,
    label="T_mean")
vlines!(ax0, [trt_3_raw.elapsed_time[findfirst(>(100), trt_3_raw.power)]], color=:gray,
    linestyle=:dash, label="heating onset")
axislegend(ax0, position=:rt)

ax1 = Axis(fig1[2, 1], xlabel="Time since start of heating (s)", ylabel="Temperature (°C)",
    title="2-phase vs 3-phase, after load_trt_data's default trim + rebase")
lines!(ax1, trt_2.elapsed_time, trt_2.T_mean, color=:green, linewidth=3,
    label="2-phase (no recirculation logged)")
lines!(ax1, trt_3.elapsed_time, trt_3.T_mean, color=:black, linewidth=1, linestyle=:dash,
    label="3-phase (recirculation trimmed)")
axislegend(ax1, position=:rt)

ax2 = Axis(fig1[3, 1], xlabel="Time since start of heating (s)", ylabel="Power (W)")
lines!(ax2, trt_2.elapsed_time, trt_2.power, color=:green, linewidth=3, label="2-phase")
lines!(ax2, trt_3.elapsed_time, trt_3.power, color=:black, linewidth=1, linestyle=:dash,
    label="3-phase")
axislegend(ax2, position=:rt)
display(fig1)

# 2. Non-uniform time steps: a Comsol numerical export

file_num = joinpath(@__DIR__, "..", "data", "Data_CL_Num.csv")

trt_num_raw = CSV.read(file_num, DataFrame; delim=',', header=1)
rename!(trt_num_raw, [1 => :t, 2 => :power, 3 => :T_in, 4 => :T_out])

trt_num = load_trt_data(file_num)   # prints a warning: adaptive-time-stepping ramp at the start
d_num = decompose_trt(trt_num)

fig2 = Figure()
ax3 = Axis(fig2[1, 1], xlabel="Time (s)", ylabel="Temperature (°C)",
    title="Numerical (Comsol) TRT, after uniform-grid interpolation")
lines!(ax3, trt_num.elapsed_time, trt_num.T_in, color=:red, linewidth=2, label="T_in")
lines!(ax3, trt_num.elapsed_time, trt_num.T_out, color=:blue, linewidth=2, label="T_out")
lines!(ax3, trt_num.elapsed_time, trt_num.T_mean, color=:green, linewidth=2, label="T_mean")
lines!(ax3, d_num.heating.elapsed_time, d_num.heating.T_mean, color=:black, linewidth=1,
    linestyle=:dot)
lines!(ax3, d_num.recovery.elapsed_time, d_num.recovery.T_mean, color=:black, linewidth=1,
    linestyle=:dot)
axislegend(ax3, position=:rt)
display(fig2)
