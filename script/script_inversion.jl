# Model inversion of a thermal response test using `fit_ground_response` directly, called once per
# ground model so it's clear what goes into each one. Every model is fitted by superimposing the
# measured load with its GroundResponse.jl g-function (Eq. 1 of Pasquier 2018).

import Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate() # Once per computer

using ThermalResponseTest
using GroundHeatExchanger
using DataFrames
using CairoMakie

# Load and decompose a thermal response test
datafile = joinpath(@__DIR__, "..", "data", "TRT_CL_Num.csv")
trt = load_trt_data(datafile)
dataset = decompose_trt(trt)

# Borehole / ground parameters
H  = 150.0          # borehole depth [m]
D  = 0.0            # buried depth [m]
rb = 0.075          # borehole radius [m]
Cs = 2.0e6          # ground volumetric heat capacity [J/m³K]
Cf = 4.2e6          # groundwater volumetric heat capacity [J/m³K]
T0 = trt.T_mean[1]  # undisturbed ground temperature [°C]

t = dataset.data.elapsed_time
T = dataset.data.T_mean
q = dataset.data.power / H

# Fit each ground model.
@time ils  = fit_ground_response(t, T, q, rb, T0, p -> ILSModel(p[1], Cs), [2.5], [0.2], [7.0])
@time ics  = fit_ground_response(t, T, q, rb, T0, p -> ICSModel(rb, p[1], Cs), [2.5], [0.2], [7.0])
@time fls  = fit_ground_response(t, T, q, rb, T0, p -> FLSModel(H, D, p[1], Cs), [2.5], [0.2], [7.0])
@time mils = fit_ground_response(t, T, q, rb, T0, p -> MILSModel(rb, p[1], Cs, Cf, p[2]),
    [2.5, 1e-7], [0.2, 1e-9], [7.0, 1e-5])
@time mfls = fit_ground_response(t, T, q, rb, T0, p -> MFLSModel(H, rb, D, p[1], Cs, Cf, p[2]),
    [2.5, 1e-7], [0.2, 1e-9], [7.0, 1e-5])

# Reconstruct each fitted temperature profile
T_ils = fluid_temperature(t, q, ground_response(t, rb, [0.0 0.0], ils.model), T0, ils.Rbₑ)
T_ics = fluid_temperature(t, q, ground_response(t, rb, [0.0 0.0], ics.model), T0, ics.Rbₑ)
T_fls = fluid_temperature(t, q, ground_response(t, rb, [0.0 0.0], fls.model), T0, fls.Rbₑ)
T_mils = fluid_temperature(t, q, ground_response(t, rb, [0.0 0.0], mils.model), T0, mils.Rbₑ)
T_mfls = fluid_temperature(t, q, ground_response(t, rb, [0.0 0.0], mfls.model), T0, mfls.Rbₑ)

# Summary table — one row per model.
rmse(T_fit) = sqrt(sum(abs2, T_fit .- T) / length(T))

summary_df = DataFrame(
    Model = ["ILS", "ICS", "FLS", "MILS", "MFLS"],
    k = round.([ils.params[1], ics.params[1], fls.params[1], mils.params[1], mfls.params[1]],
        digits = 4),
    Rbₑ = round.([ils.Rbₑ, ics.Rbₑ, fls.Rbₑ, mils.Rbₑ, mfls.Rbₑ], digits = 5),
    vD = [missing, missing, missing, round(mils.params[2], digits = 3), round(mfls.params[2], 
        digits = 3)],
    RMSE = round.([rmse(T_ils), rmse(T_ics), rmse(T_fls), rmse(T_mils), rmse(T_mfls)], digits = 4),
)
println(summary_df)

# Figure — measured vs every fitted model, overlaid for direct comparison.
fig = Figure()
ax = Axis(fig[1, 1], title = "TRT model inversion — measured vs fitted mean fluid temperature",
    xlabel = "Time since start of heating (s)", ylabel = "Temperature (°C)")
scatter!(ax, t, T, color = (:black, 0.4), markersize = 10, label = "Measured")
lines!(ax, t, T_ils, color = :dodgerblue, linewidth = 2,
    label = "ILS (k=$(round(ils.params[1], digits=2)))")
lines!(ax, t, T_ics, color = :seagreen, linewidth = 2,
    label = "ICS (k=$(round(ics.params[1], digits=2)))")
lines!(ax, t, T_fls, color = :orange, linewidth = 2,
    label = "FLS (k=$(round(fls.params[1], digits=2)))")
lines!(ax, t, T_mils, color = :purple, linestyle = :dash, linewidth = 2,
    label = "MILS (k=$(round(mils.params[1], digits=2)))")
lines!(ax, t, T_mfls, color = :firebrick, linestyle = :dot, linewidth = 2,
    label = "MFLS (k=$(round(mfls.params[1], digits=2)))")
axislegend(ax, position = :rt)
display(fig)
