"""
Application of the first-order approximation (FOA) interpretation methods of a thermal response
test, following Pasquier (2018). Three methods are demonstrated on a real TRT data set:

  1. UFOA-T-H  — temperature, heating phase            (`fit_ils_foa_T`,          Eq. 3)
  2. UFOA-T-R  — temperature, recovery phase           (`fit_ils_foa_T_recovery`, Eq. 13)
  3. CFOA-Ṫ-H  — temperature derivative, heating phase (`fit_ils_foa_dT`,         Eqs. 8–10)

Run from the package root:
    julia --project=script script/script_first_order_approximation.jl
"""

using ThermalResponseTest
using CairoMakie

# --- Load and decompose the TRT data --------------------------------------------------
datafile = joinpath(@__DIR__, "..", "trt_data", "DataCL_TRT.csv")
trt = load_trt_data(datafile)
dataset = decompose_trt(trt)

# --- Borehole / ground parameters -----------------------------------------------------
H  = 138.0          # borehole depth [m]
rb = 0.075          # borehole radius [m]
ri = 0.02           # pipe inner radius [m]
Cs = 2.0e6          # ground volumetric heat capacity [J/m³K]
V  = 30 / 60000     # volumetric flow rate [m³/s] (30 L/min)

# ======================================================================================
# 1. Heating phase — temperature  (UFOA-T-H, Eq. 3)
# ======================================================================================
k_H, Rb_H, reg_H, ind_H = fit_ils_foa_T(dataset.heating, H, rb, Cs)
tc = critical_time(rb, k_H, Cs)        # FOA only valid for t > tc
println("UFOA-T-H  (heating, temperature)   : k = $(round(k_H, digits=3)) W/mK, " *
        "Rb = $(round(Rb_H, digits=4)) mK/W,  tc = $(round(tc/3600, digits=1)) h")

# ======================================================================================
# 2. Recovery phase — temperature  (UFOA-T-R, Eq. 13)
# ======================================================================================
k_R, T0_R, reg_R, ind_R = fit_ils_foa_T_recovery(dataset, H)
println("UFOA-T-R  (recovery, temperature)  : k = $(round(k_R, digits=3)) W/mK, " *
        "T0 = $(round(T0_R, digits=2)) °C")

# ======================================================================================
# 3. Heating phase — temperature derivative  (CFOA-Ṫ-H, Eqs. 8–10)
# ======================================================================================
dT = centered_finite_difference(dataset.heating.elapsed_time, dataset.heating.T_mean)
k_dH, reg_dH, ind_dH = fit_ils_foa_dT(dataset.heating, V, H, ri)
tr = residence_time(V, H, ri)          # derivative method valid in [4 tr, 16 tr]
println("CFOA-Ṫ-H  (heating, derivative)    : k = $(round(k_dH, digits=3)) W/mK, " *
        "tr = $(round(tr/3600, digits=2)) h")

# ======================================================================================
# Figure — the three methods side by side
# ======================================================================================
fig = Figure(size = (1500, 450))

# Panel 1 — heating temperature
ax1 = Axis(fig[1, 1], title = "UFOA-T-H (heating)", xlabel = "Time (s)",
    ylabel = "Temperature (°C)", xscale = log10)
lines!(ax1, dataset.heating.elapsed_time, dataset.heating.T_mean, color = :black,
    linewidth = 2, label = "Measured")
lines!(ax1, reg_H[:, 1], reg_H[:, 2], color = :red, linewidth = 2, label = "FOA fit")
vlines!(ax1, [tc], color = :gray, linestyle = :dash, label = "tc")
axislegend(ax1, position = :rb)
text!(ax1, 0.05, 0.95, space = :relative, align = (:left, :top),
    text = "k = $(round(k_H, digits=2)) W/mK\nRb = $(round(Rb_H, digits=3)) mK/W")

# Panel 2 — recovery temperature
ax2 = Axis(fig[1, 2], title = "UFOA-T-R (recovery)", xlabel = "Time (s)",
    ylabel = "Temperature (°C)", xscale = log10)
lines!(ax2, dataset.cooling.elapsed_time, dataset.cooling.T_mean, color = :black,
    linewidth = 2, label = "Measured")
lines!(ax2, reg_R[:, 1], reg_R[:, 2], color = :blue, linewidth = 2, label = "FOA fit")
axislegend(ax2, position = :rt)
text!(ax2, 0.05, 0.95, space = :relative, align = (:left, :top),
    text = "k = $(round(k_R, digits=2)) W/mK\nT0 = $(round(T0_R, digits=2)) °C")

# Panel 3 — heating derivative
ax3 = Axis(fig[1, 3], title = "CFOA-Ṫ-H (heating derivative)", xlabel = "Time (s)",
    ylabel = "dT/dt (°C/s)", xscale = log10)
lines!(ax3, dataset.heating.elapsed_time, dT, color = :black, linewidth = 2, label = "Measured")
lines!(ax3, reg_dH[:, 1], reg_dH[:, 2], color = :green, linewidth = 3, label = "FOA fit")
vlines!(ax3, [4tr, 16tr], color = :gray, linestyle = :dash, label = "[4 tr, 16 tr]")
axislegend(ax3, position = :rt)
text!(ax3, 0.05, 0.05, space = :relative, align = (:left, :bottom),
    text = "k = $(round(k_dH, digits=2)) W/mK")

mkpath(joinpath(@__DIR__, "figures"))
outfile = joinpath(@__DIR__, "figures", "first_order_approximation.png")
save(outfile, fig)
println("Saved figure to $(outfile)")
# display(fig)   # uncomment for interactive use
