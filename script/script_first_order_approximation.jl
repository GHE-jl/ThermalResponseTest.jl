# Application of the first-order approximation (FOA) interpretation methods of a thermal response
# test, following Pasquier (2018). Four methods are demonstrated on a real TRT data set:

#   1. UFOA-T-H  — temperature, heating phase             (`fit_ils_foa_T`,           Eq. 3)
#   2. UFOA-T-R  — temperature, recovery phase            (`fit_ils_foa_T_recovery`,  Eq. 13)
#   3. CFOA-Ṫ-H  — temperature derivative, heating phase  (`fit_ils_foa_dT`,          Eqs. 8–10)
#   4. CFOA-Ṫ-R  — temperature derivative, recovery phase (`fit_ils_foa_dT_recovery`, Eqs. 14–18)

import Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate() # Once per computer

using ThermalResponseTest
using CairoMakie

# Load and decompose the TRT data.
datafile = joinpath(@__DIR__, "..", "data", "TRT_CL_Num.csv")
trt = load_trt_data(datafile)
dataset = decompose_trt(trt)

# Borehole / ground parameters
H  = 150.0          # borehole depth [m]
rb = 0.075          # borehole radius [m]
ri = 0.0137         # pipe inner radius [m]
Cs = 2.4e6          # ground volumetric heat capacity [J/m³K]
V  = 30.0 / 6e4     # volumetric flow rate [m³/s] (30 L/min)

# 1. Heating phase — temperature  (UFOA-T-H, Eq. 3)
k_H, Rbₑ_H, reg_H, ind_H = fit_ils_foa_T(dataset, H, rb, Cs)
tc = critical_time(rb, k_H, Cs)        # FOA only valid for t > tc

println("UFOA-T-H  (heating, temperature)   : k = $(round(k_H, digits=3)) W/mK, " *
        "Rbₑ = $(round(Rbₑ_H, digits=4)) mK/W,  tc = $(round(tc/3600, digits=1)) h")

# 2. Recovery phase — temperature  (UFOA-T-R, Eq. 13)
k_R, reg_R, ind_R = fit_ils_foa_T_recovery(dataset, H, rb, Cs)

println("UFOA-T-R  (recovery, temperature)  : k = $(round(k_R, digits=3)) W/mK")

# 3. Heating phase — temperature derivative  (CFOA-Ṫ-H, Eqs. 8–10)
dT = bourdet_derivative(dataset.heating.t_rel, dataset.heating.T_mean, 0.2)
k_dH, reg_dH, ind_dH = fit_ils_foa_dT(dataset, V, H, ri)

tr = residence_time(V, H, ri)          # derivative method valid in [64 tr, 512 tr] (capped to data)
println("CFOA-Ṫ-H  (heating, derivative)    : k = $(round(k_dH, digits=3)) W/mK, " *
        "tr = $(round(tr/3600, digits=2)) h")

# 4. Recovery phase — temperature derivative  (CFOA-Ṫ-R, Eqs. 14–18)
t̄ = dataset.heating.t_rel[end]  # heating duration, start of the recovery validity window
dT_R = bourdet_derivative(dataset.recovery.t_rel .- t̄, dataset.recovery.T_mean, 0.2) # t_rel .- t̄
k_dR, reg_dR, ind_dR = fit_ils_foa_dT_recovery(dataset, V, H, ri)

println("CFOA-Ṫ-R  (recovery, derivative)   : k = $(round(k_dR, digits=3)) W/mK")

# Figure — the four methods side by side
fig = Figure()

# Panel 1 — heating temperature
ax1 = Axis(fig[1, 1], title = "UFOA-T-H (heating)", xlabel = "Time since start of heating (s)",
    ylabel = "Temperature (°C)", xscale = log10)
lines!(ax1, dataset.heating.t_rel, dataset.heating.T_mean, color = :black,
    linewidth = 2, label = "Measured")
lines!(ax1, reg_H[:, 1], reg_H[:, 2], color = :red, linewidth = 2, label = "FOA fit")
vlines!(ax1, [tc], color = :gray, linestyle = :dash, label = "tc")
axislegend(ax1, position = :rb)
text!(ax1, 0.05, 0.95, space = :relative, align = (:left, :top),
    text = "k = $(round(k_H, digits=2)) W/mK\nRbₑ = $(round(Rbₑ_H, digits=3)) mK/W")

# Panel 2 — recovery temperature
ax2 = Axis(fig[1, 2], title = "UFOA-T-R (recovery)", xlabel = "Time since start of heating (s)",
    ylabel = "Temperature (°C)", xscale = log10)
lines!(ax2, dataset.recovery.t_rel, dataset.recovery.T_mean, color = :black,
    linewidth = 2, label = "Measured")
lines!(ax2, reg_R[:, 1], reg_R[:, 2], color = :blue, linewidth = 2, label = "FOA fit")
axislegend(ax2, position = :rt)
text!(ax2, 0.05, 0.95, space = :relative, align = (:left, :top),
    text = "k = $(round(k_R, digits=2)) W/mK")

# Panel 3 — heating derivative
ax3 = Axis(fig[2, 1], title = "CFOA-Ṫ-H (heating derivative)",
    xlabel = "Time since start of heating (s)", ylabel = "dT/dt (°C/s)", xscale = log10)
lines!(ax3, dataset.heating.t_rel, dT, color = :black, linewidth = 2, label = "Measured")
lines!(ax3, reg_dH[:, 1], reg_dH[:, 2], color = :green, linewidth = 3, label = "FOA fit")
vlines!(ax3, [64tr, min(512tr, dataset.heating.t_rel[end])], color = :gray, linestyle = :dash,
    label = "[64 tr, 512 tr]")
axislegend(ax3, position = :rt)
text!(ax3, 0.05, 0.05, space = :relative, align = (:left, :bottom),
    text = "k = $(round(k_dH, digits=2)) W/mK")

# Panel 4 — recovery derivative
ax4 = Axis(fig[2, 2], title = "CFOA-Ṫ-R (recovery derivative)",
    xlabel = "Time since start of heating (s)", ylabel = "-dT/dt (°C/s)", xscale = log10,
    yscale = log10)
lines!(ax4, dataset.recovery.t_rel, -dT_R, color = :black, linewidth = 2, label = "Measured")
lines!(ax4, reg_dR[:, 1], -reg_dR[:, 2], color = :purple, linewidth = 3, label = "FOA fit")
vlines!(ax4, [t̄ + 64tr, min(t̄ + 512tr, dataset.recovery.t_rel[end])], color = :gray,
    linestyle = :dash, label = "t̄ + [64 tr, 512 tr]")
axislegend(ax4, position = :rt)
text!(ax4, 0.05, 0.05, space = :relative, align = (:left, :bottom),
    text = "k = $(round(k_dR, digits=2)) W/mK")

display(fig)
