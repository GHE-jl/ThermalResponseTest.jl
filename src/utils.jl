"""
    mean_fluid_temperature(T_in, T_out, p)
    mean_fluid_temperature(T_in, T_out, method=:arithmetic)

Calculate the mean temperature of the fluid based on inlet and outlet temperatures. The function
is based on the p-linear mean of Marcotte & Pasquier (2008), which generalizes several common means
(arithmetic, geometric, harmonic, logarithmic, and pLinear).
For surface temperature averaging, all methods are suitable, as they give similar results. For depth
temperature averaging, the p-linear mean with p=-1.0000000001 is recommended as it provides a better
estimate of the effective temperature driving the heat transfer in the borehole.
# Arguments
    - `T_in`: Inlet temperature (°C)
    - `T_out`: Outlet temperature (°C)
    - `p`: Exponent of the p-linear mean for the numeric method (required for that method). The
        `:pLinear` symbol corresponds to `p = -1.0000000001`.
    - `method`: Method to calculate mean temperature. Options are:
        - `:arithmetic` -> p = 1: (T_in + T_out) / 2 (default)
        - `:logarithmic` -> p = 0: (T_out - T_in) / log(T_out / T_in)
        - `:geometric` -> p = -0.5: sqrt(T_in * T_out)
        - `:harmonic` -> p = -2: 2 / (1/T_in + 1/T_out)
        - `:pLinear` -> p = -1.0000000001: See Marcotte & Pasquier (2008) for details.
# Output
    - Mean temperature (°C)
# Reference
    - Marcotte, D., & Pasquier, P. (2008). On the estimation of thermal resistance in borehole
        thermal conductivity test. Renewable Energy, 33(11), 2407–2415.
        https://doi.org/10.1016/j.renene.2008.01.021
"""
function mean_fluid_temperature(T_in::AbstractVector{<:Real}, T_out::AbstractVector{<:Real},
    p::Real)
    function _plinear(Tin, Tout)
        Δ = abs(Tin - Tout)
        scale = max(abs(Tin), abs(Tout))
        # At T_in ≈ T_out it gives 0/0 (L'Hôpital's rule) and the mean is (T_in + T_out) / 2
        Δ ≤ sqrt(eps(typeof(float(Tin)))) * scale && return (Tin + Tout) / 2
        return (p * (abs(Tin)^(p + 1) - abs(Tout)^(p + 1))) /
            ((1 + p) * (abs(Tin)^p - abs(Tout)^p))
    end
    return _plinear.(T_in, T_out)
end
function mean_fluid_temperature(T_in::AbstractVector{<:Real}, T_out::AbstractVector{<:Real},
    method::Symbol=:arithmetic)

    if method == :arithmetic
        p = 1
    elseif method == :logarithmic
        p = 0
    elseif method == :geometric
        p = -0.5
    elseif method == :harmonic
        p = -2
    elseif method == :pLinear
        p = -1.0000000001     
    else
        error("Unsupported method: $method")
    end
    return mean_fluid_temperature(T_in, T_out, p)
end

"""
    centered_finite_difference(t, x)

Calculate the time derivative of a vector using a centered finite difference method. This is used
to compute the time derivative of the temperature in the first order approximation method.
# Arguments
    - `t`: Time vector (s)
    - `x`: Variable vector (e.g., temperature) corresponding to time vector
# Output
    - Time derivative of `x` with respect to `t`
"""
function centered_finite_difference(t::AbstractVector{<:Real}, x::AbstractVector{<:Real})
    n = length(t)
    dxdt = similar(x)
    # Forward difference for the first point
    dxdt[1] = (x[2] - x[1]) / (t[2] - t[1])
    # Backward difference for the last point
    dxdt[n] = (x[n] - x[n-1]) / (t[n] - t[n-1])
    # Centered difference for the interior points
    for i in 2:n-1
        dxdt[i] = (x[i+1] - x[i-1]) / (t[i+1] - t[i-1])
    end
    return dxdt
end

"""
    bourdet_derivative(t, T, δ=0.3)

Calculate the time derivative of a vector using the Bourdet et al. (1989) three-point formula,
adapted from well-test analysis. Instead of the immediate neighbors used by
[`centered_finite_difference`](@ref), the two points bracketing `t[i]` are taken at a fixed
natural-log spacing `δ` away (`t[i]·exp(∓δ)`), with the temperature there obtained by linear
interpolation. This trades a small amount of resolution for robustness to measurement noise and
uneven sampling, which is often preferable for real (as opposed to synthetic) TRT signals. The
output is `dT/dt`, the same physical quantity as `centered_finite_difference`, so the two are
interchangeable wherever a temperature derivative is required. The "logarithmic derivative" 
`dT/d(ln t) = t·dT/dt` used in Beier (2020) can identify the borehole-dominated/transition/
steady-flux periods of a TRT.
# Arguments
    - `t`: Time vector (s)
    - `T`: Variable vector (e.g., temperature) corresponding to time vector
    - `δ`: Natural-log spacing used to pick the bracketing points (default 0.3, following Beier
        2020, who uses 0.3–0.5). Larger `δ` gives a smoother but less locally-resolved derivative.
# Output
    - Time derivative of `T` with respect to `t`
# Reference
    - Bourdet, D., Ayoub, J.A., Pirard, Y.M., 1989. Use of pressure derivative in well-test
        interpretation. SPE Form. Eval. 4 (2), 293–302. https://doi.org/10.2118/12777-PA
    - Beier, R.A., 2020. Deconvolution and convolution methods for thermal response tests on
        borehole heat exchangers. Geothermics 86, 101786.
        https://doi.org/10.1016/j.geothermics.2019.101786 (Appendix C)
"""
function bourdet_derivative(t::AbstractVector{<:Real}, T::AbstractVector{<:Real}, δ::Real=0.3)
    n = length(t)
    t1 = clamp.(t .* exp(-δ), t[1], t[end])
    t2 = clamp.(t .* exp(δ), t[1], t[end])

    function _linear_interp(t::AbstractVector{<:Real}, x::AbstractVector{<:Real}, tq::Real)
        j = searchsortedlast(t, tq)
        j == 0 && return x[1]
        j == length(t) && return x[end]
        t[j] == tq && return x[j]
        w = (tq - t[j]) / (t[j+1] - t[j])
        return x[j] + w * (x[j+1] - x[j])
    end

    T1 = _linear_interp.(Ref(t), Ref(T), t1)
    T2 = _linear_interp.(Ref(t), Ref(T), t2)

    dX1, dX2 = t .- t1, t2 .- t
    dT1, dT2 = T .- T1, T2 .- T

    dTdt = similar(T, float(eltype(T)))
    dTdt[1] = dT2[1] / dX2[1]
    dTdt[n] = dT1[n] / dX1[n]
    for i in 2:n-1
        dTdt[i] = (dT1[i] * dX2[i] / dX1[i] + dT2[i] * dX1[i] / dX2[i]) / (dX1[i] + dX2[i])
    end
    return dTdt
end

"""
    residence_time(V̇, H)
    residence_time(V, H, ri)

Calculate the residence time of the fluid in the borehole, which is the time it takes for the fluid
to travel down and back up the borehole.
# Arguments
    - `V̇`: Fluid velocity in the borehole (m/s)
    - `V`: Volumetric flow rate (m³/s)
    - `H`: Borehole depth (m)
    - `ri`: Inner radius of the borehole (m)
# Output
    - Residence time (s)
"""
function residence_time(V̇::Real, H::Real)
    return 2 * H / V̇
end
function residence_time(V::Real, H::Real, ri::Real)
    A = π * ri^2 # Cross-sectional area of the borehole
    return residence_time(V / A, H)
end

"""
    critical_time(rb, k, Cs)

Calculate the ILS critical time `tc = 5 rb² / (k/Cs)`, the Fourier-number threshold (Fo ≳ 5) beyond
which the infinite line source (and its first-order approximation methods) become valid.
# Arguments
    - `rb`: Borehole radius (m)
    - `k`: Ground thermal conductivity (W/mK)
    - `Cs`: Volumetric heat capacity of the ground (J/m³K)
# Output
    - Critical time (s)
# Reference
    - Pasquier, P. (2018). Interpretation of the first hours of a thermal response test using the
        time derivative of the temperature. Applied Energy, 213, 56–75.
        https://doi.org/10.1016/j.apenergy.2018.01.022
"""
function critical_time(rb::Real, k::Real, Cs::Real)
    return 5 * rb^2 / (k / Cs)
end