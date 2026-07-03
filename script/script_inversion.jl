"""
Model inversion of a thermal response test with Optimization.jl (Optim.jl backend). The ground
thermal conductivity `k` is the main unknown for every model; the moving models (MILS, MFLS) also
recover the groundwater Darcy velocity `vD`. Each model is fitted by superimposing the measured load
with the corresponding GroundResponse.jl g-function (Eq. 1 of Pasquier 2018).

Run from the package root:
    julia --project=script script/script_inversion.jl
"""

using ThermalResponseTest
using GroundHeatExchanger          # ils, ics, fls, mils, mfls and convolution (for the fitted curves)
using CairoMakie

# --- Load and decompose ----------------------------------------------------------------
datafile = joinpath(@__DIR__, "..", "trt_data", "DataCL_TRT.csv")
trt = load_trt_data(datafile)
dataset = decompose_trt(trt)

# --- Parameters ------------------------------------------------------------------------
H  = 138.0          # borehole depth [m]
D  = 4.0            # buried depth [m]
rb = 0.075          # borehole radius [m]
Cs = 2.0e6          # ground volumetric heat capacity [J/m³K]
Cf = 4.2e6          # groundwater volumetric heat capacity [J/m³K]
T0 = trt.T_mean[1]  # undisturbed ground temperature [°C]

# Decimate the heating phase uniformly (keeps the equal time spacing the superposition needs)
# so the quadrature-based models (FLS, MFLS) stay fast in the optimization loop.
heating = dataset.heating
step = max(1, fld(nrow(heating), 400))
h = heating[1:step:nrow(heating), :]
t = collect(Float64, h.elapsed_time)
T = collect(Float64, h.T_mean)
q = collect(Float64, h.power ./ H)

# Helper to rebuild the fitted mean fluid temperature from a g-function vector.
predict(g, Rb) = T0 .+ Rb .* q .+ convolution(q, g)

# --- Run every model inversion ---------------------------------------------------------
println("Inverting the TRT ($(length(t)) samples) with each ground model ...\n")

r_ils  = fit_ils(t, T, q, rb, T0, Cs)
g_ils  = predict(ils(t, rb, r_ils.k, Cs), r_ils.Rb)
println("ILS  : k = $(round(r_ils.k,  digits=3)) W/mK,  Rb = $(round(r_ils.Rb,  digits=4)) mK/W")

r_ics  = fit_ics(t, T, q, rb, T0, Cs)
g_ics  = predict(ics(t, rb, rb, r_ics.k, Cs), r_ics.Rb)
println("ICS  : k = $(round(r_ics.k,  digits=3)) W/mK,  Rb = $(round(r_ics.Rb,  digits=4)) mK/W")

r_fls  = fit_fls(t, T, q, rb, H, D, T0, Cs)
g_fls  = predict(fls(t, rb, H, D, r_fls.k, Cs), r_fls.Rb)
println("FLS  : k = $(round(r_fls.k,  digits=3)) W/mK,  Rb = $(round(r_fls.Rb,  digits=4)) mK/W")

r_mils = fit_mils(t, T, q, rb, T0, Cs, Cf)
g_mils = predict(mils(t, [0.0, 0.0], rb, r_mils.k, Cs, Cf, r_mils.vD), r_mils.Rb)
println("MILS : k = $(round(r_mils.k, digits=3)) W/mK,  Rb = $(round(r_mils.Rb, digits=4)) mK/W,  " *
        "vD = $(round(r_mils.vD, sigdigits=3)) m/s")

r_mfls = fit_mfls(t, T, q, rb, H, D, T0, Cs, Cf)
g_mfls = predict(mfls(t, [0.0, 0.0], H, rb, D, r_mfls.k, Cs, Cf, r_mfls.vD), r_mfls.Rb)
println("MFLS : k = $(round(r_mfls.k, digits=3)) W/mK,  Rb = $(round(r_mfls.Rb, digits=4)) mK/W,  " *
        "vD = $(round(r_mfls.vD, sigdigits=3)) m/s")

# --- Figure — measured vs fitted -------------------------------------------------------
fig = Figure(size = (900, 600))
ax = Axis(fig[1, 1], title = "TRT model inversion — measured vs fitted mean fluid temperature",
    xlabel = "Time (s)", ylabel = "Temperature (°C)", xscale = log10)
scatter!(ax, t, T, color = (:black, 0.4), markersize = 5, label = "Measured")
lines!(ax, t, g_ils,  linewidth = 2, label = "ILS  (k=$(round(r_ils.k, digits=2)))")
lines!(ax, t, g_ics,  linewidth = 2, label = "ICS  (k=$(round(r_ics.k, digits=2)))")
lines!(ax, t, g_fls,  linewidth = 2, label = "FLS  (k=$(round(r_fls.k, digits=2)))")
lines!(ax, t, g_mils, linewidth = 2, linestyle = :dash, label = "MILS (k=$(round(r_mils.k, digits=2)))")
lines!(ax, t, g_mfls, linewidth = 2, linestyle = :dot,  label = "MFLS (k=$(round(r_mfls.k, digits=2)))")
axislegend(ax, position = :rb)

mkpath(joinpath(@__DIR__, "figures"))
outfile = joinpath(@__DIR__, "figures", "inversion.png")
save(outfile, fig)
println("\nSaved figure to $(outfile)")
# display(fig)   # uncomment for interactive use
