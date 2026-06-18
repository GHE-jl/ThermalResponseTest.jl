"""
Script that tests the first order approximation (FOA) method for estimating thermal conductivity and
borehole thermal resistance from thermal response test data.
"""

using CairoMakie

includet("../src/load_trt_data.jl")
includet("../src/first_order_approximation.jl")

# Load TRT data
file1 = "trt_data/DataCL_TRT.csv"
trt = load_trt_data(file1)

# Decompose into heating and cooling phases
trt_dataset = decompose_trt(trt)

# Apply FOA to heating phase
H = 138.0               # Example borehole depth in meters
rb = 0.075              # Example borehole radius in meters
Cs = 2e6                # Volumetric heat capacity of the ground in J/m³K
k1, Rb1, reg1, ind1 = fit_ils_foa_T(
    trt_dataset.heating.elapsed_time,
    trt_dataset.heating.T_mean,
    trt_dataset.heating.power ./ H,  # Convert total power to power per unit length
    rb,
    Cs
)
k, Rb, reg, ind = fit_ils_foa_T(trt_dataset.heating, H, rb, Cs)

# Fit
Texp = trt_dataset.heating.T_mean[ind]
rmse = sqrt(sum((reg[:, 2] .- Texp) .^ 2)/length(Texp))
tc = 5 * rb^2 / (k / Cs)                # Critical time

@assert k1 ≈ k "FOA results from direct data and DataFrame input should match"
@assert Rb1 ≈ Rb "FOA results from direct data and DataFrame input should match"
@assert reg1 ≈ reg "FOA regression data from direct data and DataFrame input should match"

# Plot results
fig = Figure()
ax1 = Axis(fig[1, 1], title="FOA Fit - Heating Phase", xlabel="Time (s)", ylabel="Temperature (°C)",
    xscale=log10)
lines!(ax1, trt_dataset.heating.elapsed_time, trt_dataset.heating.T_mean, color=:black, 
    linewidth=2, label="Measured")
lines!(ax1, reg[:, 1], reg[:, 2], color=:red, linewidth=2, label="FOA Fit")
vlines!(ax1, [tc], color=:black, linestyle=:dash, label="Critical Time")
axislegend(ax1, position=:rb)
text!(ax1, 0.75, 0.5, text=
"k: $(round(k, digits=2)) W/mK
Rb: $(round(Rb, digits=3)) mK/W
RMSE: $(round(rmse, digits=2)) °C",
    space=:relative, align=(:left, :top))
text!(ax1, tc + 3600, minimum(trt_dataset.heating.T_mean) + 2, text=
    "tc = $(round(tc / 3600, digits=1)) h", align=(:left, :top))
display(fig)