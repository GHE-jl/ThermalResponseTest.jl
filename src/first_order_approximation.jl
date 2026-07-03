using LinearAlgebra
using DataFrames: AbstractDataFrame

"""
    fit_ils_foa_T(t, T, q, r_b, Cs)
    fit_ils_foa_T(t, T, Q, H, r_b, Cs=2e6)
    fit_ils_foa_T(trt, H, r_b, Cs)

Fit the first order approximation (FOA) of the infinite line source (ILS) model to experimental
temperature from a thermal response test. The function iteratively estimates the effective thermal
conductivity (k) and borehole thermal resistance (Rb). This corresponds to the unconstrained FOA
of temperature during the heating phase (UFOA-T-H) as described in Pasquier (2018) (Eq. 3).
# Arguments
    - `t`: Time vector [s]
    - `Temp`: Temperature vector [°C]
    - `Q`: Total heat injection rate [W]
    - `q`: Heat injection rate per borehole length (Q/H) [W/m]
    - `H`: Borehole depth [m]
    - `rb`: Borehole radius [m]
    - `Cs`: Volumetric heat capacity of the ground [J/m³K]
    - `trt`: DataFrame containing the TRT data with columns `:elapsed_time`, `:T_mean`, and `:power`
# Output
    - `k`: Estimated effective thermal conductivity [W/mK]
    - `Rb`: Estimated borehole thermal resistance [mK/W]
    - `reg`: Matrix [N, 2] with the time and fitted temperature values from the FOA regression
    - `indices`: Indices of the data points used in the regression
# Reference
    - Austin, W. A. (1998). DEVELOPMENT OF AN IN SITU SYSTEM FOR MEASURING GROUND THERMAL PROPERTIES
        [Master of Science, Oklahoma State University].
        https://hvac.okstate.edu/sites/default/files/Austin_thesis.pdf
    - Pasquier, P. (2018). Interpretation of the first hours of a thermal response test using the
        time derivative of the temperature. Applied Energy, 213, 56–75.
        https://doi.org/10.1016/j.apenergy.2018.01.022
"""
function fit_ils_foa_T(t::AbstractVector{<:Real}, T::AbstractVector{<:Real},
    q::AbstractVector{<:Real}, rb::Real, Cs::Real)
    # Initial guesses
    tᵢ = Inf
    kᵢ = 2.5
    tc = 5 * rb^2 / (kᵢ / Cs)                # Critical time
    iₘ = 25
    i = 0
    
    # Pre-allocate for regression results
    slope, intercept, indices = 0.0, 0.0, Int[]

    # Effective thermal conductivity estimation loop
    while abs(tc - tᵢ) > 1 && i < iₘ
        # Find indices where t > tc
        indices = findall(>(tc), t)
        # println("Iteration $i: k = $(round(kᵢ, digits=2)) W/mK, tc = $(round(tc, digits=2)) s")
        
        # Divergence handling
        if length(indices) ≤ 2  # Require at least 5 points for regression
            n = length(t)
            start_idx = max(1, n - round(Int, 0.3 * n)) # Use last 30% of data
            indices = collect(start_idx:n)
            @warn "tc diverges at iteration $i"
        end

        # Perform Linear Regression: T = slope * log(t) + intercept
        ln_t = log.(view(t, indices))
        X = [ln_t ones(length(ln_t))]
        y = view(T, indices)
        
        coef = X \ y
        # coef = pinv(X) * y  # Use pseudo-inverse for stability if needed
        slope, intercept = coef

        # Update k and tc
        kᵢ = sum(view(q, indices)) / (4 * pi * slope * length(indices)) # k = q / (4π * slope)
        tᵢ = tc
        tc = 5 * rb^2 / (kᵢ / Cs)           # Update critical time
        i += 1
    end

    # Effective borehole thermal resistance estimation
    γ = Base.MathConstants.eulergamma
    mean_q = sum(q) / length(q)
    Rb = (intercept - T[1] - (slope * (log((4 * kᵢ / Cs) / (rb^2)) - γ))) / mean_q

    # Handle complex results
    k = real(kᵢ)
    if !isreal(kᵢ)
        @warn "Complex result in FOA; taking real part."
    end

    # Prepare regression output data
    reg = hcat(t[indices], slope .* log.(t[indices]) .+ intercept)

    return k, Rb, reg, indices
end
function fit_ils_foa_T(trt::AbstractDataFrame, H::Real, rb::Real, Cs::Real)
    q = trt.power ./ H                      # Convert total power to power per unit length
    return fit_ils_foa_T(trt.elapsed_time, trt.T_mean, q, rb, Cs)
end

"""
    fit_ils_foa_T_recovery(t, T, q, t̄)
    fit_ils_foa_T_recovery(dataset::TRTDataset, H)

Fit the unconstrained first order approximation (FOA) of the infinite line source (ILS) to the
temperature measured during the **recovery phase** of a thermal response test, after the heating
power has been switched off while circulation continues. This is the UFOA-T-R method described in
Pasquier (2018) (Eq. 13).

Setting a null heating power for `t > t̄` (with `t̄` the heating phase duration) and applying the
ILS superposition principle, the recovery-phase fluid temperature simplifies to
`Tf(t) = T₀ + q/(4π·k)·ln(t / (t - t̄))`. A linear regression of `T` against
`x = ln(t / (t - t̄))` then has slope `m = q/(4π·k)` and intercept `T₀`, so the ground thermal
conductivity is `k = q / (4π·m)`. Unlike the heating-phase method, the recovery phase carries no
information on the borehole resistance.
# Arguments
    - `t`: Time vector measured from the **start of heating** [s] (recovery portion, `t > t̄`)
    - `T`: Mean fluid temperature during recovery [°C]
    - `q`: Mean heat injection rate per borehole length during heating (Q/H) [W/m] (scalar)
    - `t̄`: Heating phase duration [s]
    - `dataset`: `TRTDataset` from `decompose_trt` (heating + cooling phases)
    - `H`: Borehole depth [m]
# Output
    - `k`: Estimated effective thermal conductivity [W/mK]
    - `T0`: Estimated undisturbed ground temperature (regression intercept) [°C]
    - `reg`: Matrix [N, 2] with the time and fitted temperature from the FOA regression
    - `indices`: Indices of the recovery data points used in the regression
# Reference
    - Pasquier, P. (2018). Interpretation of the first hours of a thermal response test using the
        time derivative of the temperature. Applied Energy, 213, 56–75.
        https://doi.org/10.1016/j.apenergy.2018.01.022
"""
function fit_ils_foa_T_recovery(t::AbstractVector{<:Real}, T::AbstractVector{<:Real},
    q::Real, t̄::Real)
    # Regression variable x = ln(t / (t - t̄)); valid only strictly after heating stops.
    indices = findall(tᵢ -> tᵢ > t̄, t)
    if length(indices) < 5
        throw(ArgumentError("Not enough recovery data points (t > t̄) for UFOA-T-R."))
    end
    x = log.(view(t, indices) ./ (view(t, indices) .- t̄))
    X = [x ones(length(x))]
    y = view(T, indices)
    coef = X \ y
    slope, intercept = coef                 # slope = q/(4π·k), intercept = T₀

    k = q / (4 * pi * slope)
    if !isfinite(k) || k <= 0
        throw(ArgumentError("UFOA-T-R produced a non-physical thermal conductivity estimate."))
    end
    T0 = intercept

    reg = hcat(t[indices], slope .* x .+ intercept)
    return k, T0, reg, indices
end
function fit_ils_foa_T_recovery(dataset::TRTDataset, H::Real)
    nrow(dataset.cooling) > 0 || throw(ArgumentError("Dataset has no recovery (cooling) phase."))
    t̄ = dataset.heating.elapsed_time[end]                       # heating duration
    q = (sum(dataset.heating.power) / nrow(dataset.heating)) / H  # mean heating power per length
    return fit_ils_foa_T_recovery(dataset.cooling.elapsed_time, dataset.cooling.T_mean, q, t̄)
end

"""
    fit_ils_foa_dT(t, dT, q, indices)
    fit_ils_foa_dT(t, dT, q, V, H, ri)
    fit_ils_foa_dT(trt, V, H, ri)

Fit a linear regression of the infinite line source derivative (ILSd) model to the time
derivative of the temperature during the heating phase of a thermal response test. This corresponds
to the constrained FOA method CFOA-Ṫ-H as described in Pasquier (2018) (Eqs. 4 and 8 to 12). This 
method uses indices as inputs, typically corresponding to 4 and 16 times the residence time (tr) to
constrain the regression to a specific valid time range due to measurement errors in the time
derivative of the temperature.
# Arguments
    - `t`: Time vector [s]
    - `dT`: Time derivative of the temperature vector [°C/s]
    - `q`: Heat injection rate per borehole length (Q/H) [W/m]
    - `indices`: Indices of the data points to use in the regression
    - `ind`: Pair of indices corresponding to the valid time range for the regression
        - e.g., [t1, t2]
    - `V`: Volumetric flow rate (m³/s)
    - `H`: Borehole depth (m)
    - `ri`: Inner radius of the borehole (m)
    - `trt`: DataFrame containing the TRT data with columns `:elapsed_time`, `:T_mean`, and `:power`
# Output
    - `k`: Estimated effective thermal conductivity [W/mK]
    - `reg`: Matrix [N, 2] with the time and fitted temperature derivative from the FOA regression
    - `indices`: Indices of the data points used in the regression
"""
function fit_ils_foa_dT(t::AbstractVector{<:Real}, dT::AbstractVector{<:Real},
    q::AbstractVector{<:Real}, indices::AbstractVector{<:Integer})
    # Evaluate b̃ (Eq. 9 in Pasquier 2018) and update k (Eq. 10 in Pasquier 2018)
    b̃ = sum(log.(view(t, indices)) .+ log.(view(dT, indices))) / length(indices)
    avg_q = sum(view(q, indices)) / length(indices)
    k = avg_q / (4 * pi * exp(b̃))
    if !isfinite(k) || k <= 0
        throw(ArgumentError("CFOA-Ṫ-H produced a non-physical thermal conductivity estimate."))
    end
    # Prepare regression output data
    Tdot_fit = (avg_q / (4 * pi * k)) ./ t[indices] # Regression (Eq. 6 in Pasquier 2018)
    reg = hcat(t[indices], Tdot_fit)
    return k, reg, indices
end
function fit_ils_foa_dT(t::AbstractVector{<:Real}, dT::AbstractVector{<:Real},
    q::AbstractVector{<:Real}, V::Real, H::Real, ri::Real)
    v = V / (π * ri^2)                     # Calculate fluid velocity from volumetric flow rate and inner radius
    tr = residence_time(v, H)               # Calculate residence time
    tc1 = 4 * tr                            # First critical time for FOA validity
    tc2 = 16 * tr                           # Second critical time for FOA validity
    indices = findall(x -> x > tc1 && x < tc2, t) # Find indices corresponding to valid time range
    if length(indices) < 5
        throw(ArgumentError("Not enough data points in the valid time range for CFOA-Ṫ-H."))
    end
    return fit_ils_foa_dT(t, dT, q, indices)
end
function fit_ils_foa_dT(trt::AbstractDataFrame, V::Real, H::Real, ri::Real)
    # Compute dT/dt from measured mean fluid temperature before CFOA optimization.
    dT = centered_finite_difference(trt.elapsed_time, trt.T_mean)
    q = trt.power ./ H
    return fit_ils_foa_dT(trt.elapsed_time, dT, q, V, H, ri)
end