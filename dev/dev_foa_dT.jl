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

# Inputs
V = 30 / 60000          # Example volumetric flow rate [m³/s] (30 L/min)
H = 138.0               # Example borehole depth in meters
rb = 0.075              # Example borehole radius in meters
ri = 0.02               # Example inner radius of the borehole in meters
Cs = 2e6                # Volumetric heat capacity of the ground in J/m³K

# Decompose into heating and cooling phases
trt_dataset = decompose_trt(trt)
dT = centered_finite_difference(trt_dataset.heating.elapsed_time, trt_dataset.heating.T_mean)
dT_mask = copy(dT)
dT_mask[dT .< 0] .= 0.0  # Set negative derivatives to zero for FOA fitting

# Identifying residence time and indices for FOA fitting
tr = residence_time(V, H, ri)
tr1 = 4 * tr
tr2 = 16 * tr
indices = findall(x -> x > tr1 && x < tr2, trt_dataset.heating.elapsed_time)

# Apply FOA to heating phase
k, reg, indices = fit_ils_foa_dT(
    trt_dataset.heating.elapsed_time,
    dT,
    trt_dataset.heating.power ./ H,  # Convert total power to power per unit length
    indices)
k1, reg1, indices1 = fit_ils_foa_dT(
    trt_dataset.heating.elapsed_time,
    dT,
    trt_dataset.heating.power ./ H,  # Convert total power to po wer per unit length
    V, H, ri)
k2, reg2, indices2 = fit_ils_foa_dT(trt_dataset.heating, V, H, ri)

@assert k ≈ k1 ≈ k2 "FOA thermal conductivity estimates should be consistent across methods."
@assert reg ≈ reg1 ≈ reg2 "FOA regression outputs should be consistent across methods."
@assert indices ≈ indices1 ≈ indices2 "FOA regression indices should be consistent across methods."

# Fit
rmse = sqrt(sum((reg[:, 2] .- dT[indices]) .^ 2) / length(indices))

# Plot results
fig = Figure()
ax1 = Axis(fig[1, 1], title="FOA Fit - Heating Phase", xlabel="Time (s)",
    ylabel="Temperature Derivative (°C/s)", xscale=log10)#, yscale=log10)
lines!(ax1, trt_dataset.heating.elapsed_time, dT, color=:black, 
    linewidth=2, label="Measured")
lines!(ax1, reg[:, 1], reg[:, 2], color=:red, linewidth=2, label="FOA Fit")
vlines!(ax1, [tr1, tr2], color=:black, linestyle=:dash, label="Critical Time")
axislegend(ax1, position=:rc)
text!(ax1, 0.1, 0.1, text="k: $(round(k, digits=2)) W/mK\nRMSE: $(round(rmse, digits=4)) °C/s",
    space=:relative, align=(:left, :bottom))
text!(ax1, tr1 .+ 60, dT[1], text="$(round(tr1 / 3600, digits=1)) h", align=(:left, :top))
text!(ax1, tr2 .+ 3600, dT[1], text="$(round(tr2 / 3600, digits=1)) h", align=(:left, :top))
display(fig)