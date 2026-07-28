using LinearAlgebra

"""
    fit_ils_foa_T(t, T, q, rb, T0, Cs=2e6)
    fit_ils_foa_T(dataset::TRTDataset, H, rb, Cs=2e6)

Fit the first order approximation (FOA) of the infinite line source (ILS) model to experimental
temperature from a thermal response test. The function iteratively estimates the effective thermal
conductivity (k) and effective borehole thermal resistance (Rbₑ). This corresponds to the
unconstrained FOA of temperature during the heating phase (UFOA-T-H) as described in Pasquier (2018)
(Eq. 3).

**Important**: `t` must be measured from the start of the heating phase (`t = 0` the instant the
heater turns on). `T0` is the undisturbed ground temperature, known ahead of the regression (e.g.
the sample right before heating starts, `dataset.data.T_mean[1]`) — it is not something this method
can estimate from the heating phase alone.
# Arguments
    - `t`: Time vector, from the start of heating [s]
    - `T`: Temperature vector [°C]
    - `q`: Heat injection rate per borehole length (Q/H) [W/m]
    - `rb`: Borehole radius [m]
    - `T0`: Undisturbed ground temperature [°C]
    - `Cs`: (default 2e6) Volumetric heat capacity of the ground [J/m³K]
    - `H`: Borehole depth [m]
    - `dataset`: `TRTDataset` from `decompose_trt` (uses the heating phase; `T0` is taken from
        `dataset.data.T_mean[1]`)
# Output
    - `k`: Estimated effective thermal conductivity [W/mK]
    - `Rbₑ`: Estimated effective borehole thermal resistance [mK/W]
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
    q::AbstractVector{<:Real}, rb::Real, T0::Real, Cs::Real=2e6)
    tᵢ = Inf
    kᵢ = 2.5
    tc = critical_time(rb, kᵢ, Cs)
    iₘ = 25
    i = 0

    # Pre-allocate for regression results
    slope, intercept, indices = 0.0, 0.0, Int[]

    # Effective thermal conductivity estimation loop
    while abs(tc - tᵢ) > 1 && i < iₘ
        # Find indices where t > tc
        indices = findall(>(tc), t)

        # Divergence handling
        if length(indices) ≤ 2  # Handle diversion in fitting
            n = length(t)
            start_idx = max(1, n - round(Int, 0.3 * n)) # Use last 30% of data
            indices = collect(start_idx:n)
            @warn "tc diverges at iteration $i"
        end

        # Perform Linear Regression: T = slope * log(t) + intercept
        ln_t = log.(view(t, indices))
        X = [ln_t ones(length(ln_t))]
        y = view(T, indices)
        slope, intercept = X \ y

        # Update k and tc
        kᵢ = sum(view(q, indices)) / (4 * pi * slope * length(indices)) # k = q / (4π * slope)
        tᵢ = tc
        tc = critical_time(rb, kᵢ, Cs)
        i += 1
    end

    if !isfinite(kᵢ) || kᵢ <= 0
        throw(ArgumentError("UFOA-T-H produced a non-physical thermal conductivity estimate."))
    end
    k = kᵢ

    # Effective borehole thermal resistance estimation (Eq. 3 in Pasquier 2018)
    γ = Base.MathConstants.eulergamma
    mean_q = sum(q) / length(q)
    Rbₑ = (intercept - T0 - (slope * (log((4 * k / Cs) / (rb^2)) - γ))) / mean_q

    # Prepare regression output data
    reg = hcat(t[indices], slope .* log.(t[indices]) .+ intercept)

    return k, Rbₑ, reg, indices
end
function fit_ils_foa_T(dataset::TRTDataset, H::Real, rb::Real, Cs::Real=2e6)
    heating = dataset.heating
    q = heating.power ./ H                  # Convert total power to power per unit length
    T0 = dataset.data.T_mean[1]             # Undisturbed ground temperature, known from loading
    return fit_ils_foa_T(heating.t_rel, heating.T_mean, q, rb, T0, Cs)
end

"""
    fit_ils_foa_T_recovery(t, T, q, rb, t̄, Cs=2.4e6)
    fit_ils_foa_T_recovery(dataset::TRTDataset, H, rb, Cs=2.4e6)

Fit the unconstrained first order approximation (FOA) of the infinite line source (ILS) to the
temperature measured during the **recovery phase** of a thermal response test, after the heating
power has been switched off while circulation continues. This is the UFOA-T-R method described in
Pasquier (2018) (Eq. 13). Like [`fit_ils_foa_T`](@ref), the function iteratively estimates the
critical time `tc` (offset by the heating duration `t̄`) and restricts the regression to recovery
data past it, since the first-order approximation is only valid there.

**Important**: `t` and `t̄` must be measured from the start of the **heating** phase, not from the
start of the recovery phase or the raw log (see [`fit_ils_foa_T`](@ref)). Unlike the heating phase,
this regression does not need a known `T0`: its intercept directly estimates it.
# Arguments
    - `t`: Time vector measured from the **start of heating** [s] (recovery portion, `t > t̄`)
    - `T`: Mean fluid temperature during recovery [°C]
    - `q`: Mean heat injection rate per borehole length during heating (Q/H) [W/m] (scalar)
    - `rb`: Borehole radius (m)
    - `t̄`: Heating phase duration [s]
    - `Cs`: (default 2.4e6) Volumetric heat capacity of the ground [J/m³K]
    - `dataset`: `TRTDataset` from `decompose_trt` (heating + recovery phases)
    - `H`: Borehole depth [m]
# Output
    - `k`: Estimated effective thermal conductivity [W/mK]
    - `reg`: Matrix [N, 2] with the time and fitted temperature from the FOA regression
    - `indices`: Indices of the recovery data points used in the regression
# Reference
    - Pasquier, P. (2018). Interpretation of the first hours of a thermal response test using the
        time derivative of the temperature. Applied Energy, 213, 56–75.
        https://doi.org/10.1016/j.apenergy.2018.01.022
"""
function fit_ils_foa_T_recovery(t::AbstractVector{<:Real}, T::AbstractVector{<:Real},
    q::Real, rb::Real, t̄::Real, Cs::Real=2.4e6)
    if count(tᵢ -> tᵢ > t̄, t) < 5
        throw(ArgumentError("Not enough recovery data points (t > t̄) for UFOA-T-R."))
    end

    tᵢ = Inf
    kᵢ = 2.5
    tc = critical_time(rb, kᵢ, Cs) + t̄
    iₘ = 25
    i = 0

    # Pre-allocate for regression results
    slope, intercept, indices = 0.0, 0.0, Int[]

    # Effective thermal conductivity estimation loop (mirrors fit_ils_foa_T, offset by t̄)
    while abs(tc - tᵢ) > 1 && i < iₘ
        # Find indices where t > tc (i.e. past the recovery-phase critical time)
        indices = findall(tⱼ -> tⱼ > tc, t)

        # Divergence handling
        if length(indices) ≤ 2  # Handle divergence in fitting
            n = length(t)
            start_idx = max(1, n - round(Int, 0.3 * n)) # Use last 30% of data
            indices = collect(start_idx:n)
            @warn "tc diverges at iteration $i"
        end

        # Perform Linear Regression: T = slope * ln(t / (t - t̄)) + intercept
        x = log.(view(t, indices) ./ (view(t, indices) .- t̄))
        X = [x ones(length(x))]
        y = view(T, indices)
        slope, intercept = X \ y

        # Update k and tc
        kᵢ = q / (4 * pi * slope)
        tᵢ = tc
        tc = critical_time(rb, kᵢ, Cs) + t̄
        i += 1
    end

    if !isfinite(kᵢ) || kᵢ <= 0
        throw(ArgumentError("UFOA-T-R produced a non-physical thermal conductivity estimate."))
    end
    k = kᵢ

    # Prepare regression output data
    x = log.(t[indices] ./ (t[indices] .- t̄))
    reg = hcat(t[indices], slope .* x .+ intercept)
    return k, reg, indices
end
function fit_ils_foa_T_recovery(dataset::TRTDataset, H::Real, rb::Real, Cs::Real=2.4e6)
    nrow(dataset.recovery) > 0 || throw(ArgumentError("Dataset has no recovery phase."))
    t̄ = dataset.heating.t_rel[end]                                # heating duration
    q = (sum(dataset.heating.power) / nrow(dataset.heating)) / H  # mean heating power per length
    return fit_ils_foa_T_recovery(dataset.recovery.t_rel, dataset.recovery.T_mean, q, rb, t̄, Cs)
end

"""
    fit_ils_foa_dT(t, dT, q, indices)
    fit_ils_foa_dT(t, dT, q, tr)
    fit_ils_foa_dT(t, dT, q, V, H, ri)
    fit_ils_foa_dT(dataset::TRTDataset, V, H, ri; derivative_method=:bourdet, δ=0.2)

Fit a linear regression of the infinite line source derivative (ILSd) model to the time
derivative of the temperature during the heating phase of a thermal response test. This corresponds
to the constrained FOA method CFOA-Ṫ-H as described in Pasquier (2018) (Eqs. 4 and 8 to 12). This
method uses indices as inputs, defaulting to between 64 and 512 times the residence time (tr, capped
to the end of the available data): Pasquier's original 4–16 tr window is only "purely arbitrary,"
chosen for being computable from readily available parameters, and on real/numerical TRT signals the
derivative's log-log trend often has not yet settled onto its asymptotic unit slope by 16 tr — the
wider window lets the regression run past that transient. A power user can override this by calling
with an explicit `indices` or `tr`.

**Important**: `t` must be measured from the start of the heating phase.
# Arguments
    - `t`: Time vector, from the start of heating [s]
    - `dT`: Time derivative of the temperature vector [°C/s]. Either
        [`centered_finite_difference`](@ref) or [`bourdet_derivative`](@ref) can be used to obtain
        it — both return the same physical quantity (`dT/dt`) and are interchangeable here. The
        Bourdet derivative trades some local resolution for robustness to measurement noise, which
        is preferable on real (as opposed to synthetic) TRT signals, and is what lets the wider
        default window below be fit at all without the raw finite difference's sign-flipping noise
        breaking the regression.
    - `q`: Heat injection rate per borehole length (Q/H) [W/m]
    - `indices`: Indices of the data points to use in the regression corresponding to between 64 and
        512 fluid residence times in the ground heat exchanger.
    - `tr`: Fluid residence time inside the ground heat exchanger [s]
    - `V`: Volumetric flow rate (m³/s)
    - `H`: Borehole depth (m)
    - `ri`: Inner radius of the borehole (m)
    - `dataset`: `TRTDataset` from `decompose_trt` (uses the heating phase)
    - `derivative_method`: (default `:bourdet`) Method used to compute `dT` from
        `dataset.heating.T_mean` — `:bourdet` for [`bourdet_derivative`](@ref) or `:centered` for
        [`centered_finite_difference`](@ref).
    - `δ`: (default 0.2) Smoothing parameter forwarded to [`bourdet_derivative`](@ref); unused when
        `derivative_method = :centered`.
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
        throw(ArgumentError("CFOA-Ṫ-H produced a non-physical thermal conductivity estimate."))
    end
    # Prepare regression output data
    Tdot_fit = (avg_q / (4 * pi * k)) ./ t[indices] # Regression (Eq. 6 in Pasquier 2018)
    reg = hcat(t[indices], Tdot_fit)
    return k, reg, indices
end
function fit_ils_foa_dT(t::AbstractVector{<:Real}, dT::AbstractVector{<:Real},
    q::AbstractVector{<:Real}, tr::Real)
    tc1 = 64 * tr                            # First critical time for FOA validity
    tc2 = min(512 * tr, t[end])              # Second critical time, capped to the available data
    indices = findall(x -> x > tc1 && x <= tc2, t)
    if length(indices) < 5
        throw(ArgumentError("Not enough data points in the valid time range for CFOA-Ṫ-H."))
    end
    return fit_ils_foa_dT(t, dT, q, indices)
end
function fit_ils_foa_dT(t::AbstractVector{<:Real}, dT::AbstractVector{<:Real},
    q::AbstractVector{<:Real}, V::Real, H::Real, ri::Real)
    v = V / (π * ri^2)                      # Fluid velocity from the volumetric flow rate
    tr = residence_time(v, H)                # Fluid residence time
    return fit_ils_foa_dT(t, dT, q, tr)
end
function fit_ils_foa_dT(dataset::TRTDataset, V::Real, H::Real, ri::Real;
    derivative_method::Symbol=:bourdet, δ::Real=0.2)
    heating = dataset.heating
    # Compute dT/dt from measured mean fluid temperature before CFOA optimization.
    dT = if derivative_method == :centered
        centered_finite_difference(heating.t_rel, heating.T_mean)
    elseif derivative_method == :bourdet
        bourdet_derivative(heating.t_rel, heating.T_mean, δ)
    else
        error("Unsupported derivative_method: $derivative_method. Use :centered or :bourdet.")
    end
    q = heating.power ./ H
    return fit_ils_foa_dT(heating.t_rel, dT, q, V, H, ri)
end

"""
    fit_ils_foa_dT_recovery(t, dT, q, t̄, indices)
    fit_ils_foa_dT_recovery(t, dT, q, t̄, tr)
    fit_ils_foa_dT_recovery(t, dT, q, t̄, V, H, ri)
    fit_ils_foa_dT_recovery(dataset::TRTDataset, V, H, ri; derivative_method=:bourdet, δ=0.2)

Fit a linear regression of the infinite line source derivative (ILSd) model to the time derivative
of the temperature during the **recovery phase** of a thermal response test, under the constraint
of a negative unit slope. This is the CFOA-Ṫ-R method described in Pasquier (2018) (Eqs. 14-18),
the recovery-phase counterpart of [`fit_ils_foa_dT`](@ref) (CFOA-Ṫ-H) — see that docstring for why
the default valid window is 64–512 tr (capped to the available data) rather than Pasquier's original
4–16 tr, and why `:bourdet` is the default derivative method. When called on a `dataset`,
`:bourdet` is evaluated on time-since-heater-off rather than on `t_rel` directly, since its
neighbor bracket is multiplicative (unlike a centered difference) and would otherwise be dominated
by the heating-phase offset `t̄`.

**Important**: `t` (and `t̄`) must be measured from the start of the **heating** phase.
# Arguments
    - `t`: Time vector measured from the **start of heating** [s] (recovery portion, `t > t̄`)
    - `dT`: Time derivative of the recovery-phase temperature vector [°C/s] (negative). Either
        [`centered_finite_difference`](@ref) or [`bourdet_derivative`](@ref) can be used to obtain
        it — both return the same physical quantity (`dT/dt`) and are interchangeable here. The
        Bourdet derivative trades some local resolution for robustness to measurement noise, which
        is preferable on real (as opposed to synthetic) TRT signals.
    - `q`: Mean heat injection rate per borehole length during heating (Q/H) [W/m] (scalar)
    - `t̄`: Heating phase duration [s]
    - `indices`: Indices of the data points to use in the regression
    - `tr`: Fluid residence time inside the ground heat exchanger [s]
    - `V`: Volumetric flow rate (m³/s)
    - `H`: Borehole depth (m)
    - `ri`: Inner radius of the borehole (m)
    - `dataset`: `TRTDataset` from `decompose_trt` (uses the heating and recovery phases)
    - `derivative_method`: (default `:bourdet`) Method used to compute `dT` from
        `dataset.recovery.T_mean` — `:bourdet` for [`bourdet_derivative`](@ref) or `:centered` for
        [`centered_finite_difference`](@ref).
    - `δ`: (default 0.2) Smoothing parameter forwarded to [`bourdet_derivative`](@ref); unused when
        `derivative_method = :centered`.
# Output
    - `k`: Estimated effective thermal conductivity [W/mK]
    - `reg`: Matrix [N, 2] with the time and fitted temperature derivative from the FOA regression
    - `indices`: Indices of the recovery data points used in the regression
# Reference
    - Pasquier, P. (2018). Interpretation of the first hours of a thermal response test using the
        time derivative of the temperature. Applied Energy, 213, 56–75.
        https://doi.org/10.1016/j.apenergy.2018.01.022
"""
function fit_ils_foa_dT_recovery(t::AbstractVector{<:Real}, dT::AbstractVector{<:Real},
    q::Real, t̄::Real, indices::AbstractVector{<:Integer})
    # Evaluate b̃ (Eq. 17 in Pasquier 2018) and update k (Eq. 18 in Pasquier 2018)
    tw = view(t, indices)
    b̃ = sum(log.(-view(dT, indices)) .+ log.(tw .* (tw .- t̄))) / length(indices)
    k = q * t̄ / (4 * pi * exp(b̃))
    if !isfinite(k) || k <= 0
        throw(ArgumentError("CFOA-Ṫ-R produced a non-physical thermal conductivity estimate."))
    end
    # Prepare regression output data (Eq. 15 in Pasquier 2018)
    Tdot_fit = -(q * t̄ / (4 * pi * k)) ./ (tw .* (tw .- t̄))
    reg = hcat(t[indices], Tdot_fit)
    return k, reg, indices
end
function fit_ils_foa_dT_recovery(t::AbstractVector{<:Real}, dT::AbstractVector{<:Real},
    q::Real, t̄::Real, tr::Real)
    tc1 = t̄ + 64 * tr                            # First critical time for FOA validity
    tc2 = min(t̄ + 512 * tr, t[end])              # Second critical time, capped to the available data
    indices = findall(x -> x > tc1 && x <= tc2, t)
    if length(indices) < 5
        throw(ArgumentError("Not enough data points in the valid time range for CFOA-Ṫ-R."))
    end
    return fit_ils_foa_dT_recovery(t, dT, q, t̄, indices)
end
function fit_ils_foa_dT_recovery(t::AbstractVector{<:Real}, dT::AbstractVector{<:Real},
    q::Real, t̄::Real, V::Real, H::Real, ri::Real)
    v = V / (π * ri^2)                          # Fluid velocity from the volumetric flow rate
    tr = residence_time(v, H)                    # Fluid residence time
    return fit_ils_foa_dT_recovery(t, dT, q, t̄, tr)
end
function fit_ils_foa_dT_recovery(dataset::TRTDataset, V::Real, H::Real, ri::Real;
    derivative_method::Symbol=:bourdet, δ::Real=0.2)
    nrow(dataset.recovery) > 0 || throw(ArgumentError("Dataset has no recovery phase."))
    t̄ = dataset.heating.t_rel[end]
    q = (sum(dataset.heating.power) / nrow(dataset.heating)) / H  # mean heating power per length
    # Compute dT/dt from measured mean fluid temperature before CFOA optimization. Bourdet's
    # neighbor bracket t·exp(∓δ) is multiplicative, not shift-invariant like a centered difference
    # is, so it must be evaluated on time-since-heater-off (t_rel .- t̄), not on t_rel directly —
    # otherwise the huge t̄ offset dominates the bracket and biases/clamps the derivative.
    dT = if derivative_method == :centered
        centered_finite_difference(dataset.recovery.t_rel, dataset.recovery.T_mean)
    elseif derivative_method == :bourdet
        bourdet_derivative(dataset.recovery.t_rel .- t̄, dataset.recovery.T_mean, δ)
    else
        error("Unsupported derivative_method: $derivative_method. Use :centered or :bourdet.")
    end
    return fit_ils_foa_dT_recovery(dataset.recovery.t_rel, dT, q, t̄, V, H, ri)
end
