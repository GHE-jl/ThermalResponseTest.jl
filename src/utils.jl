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
    critical_time(rb, k, Cs)

Calculate the critical time (tc) for a borehole thermal response test, which is the time after which
the infinite line source (ILS) model becomes valid.
# Arguments
    - `rb`: Borehole radius (m)
    - `k`: Thermal conductivity of the ground (W/mK)
    - `Cs`: Volumetric heat capacity of the ground (J/m³K)
# Output
    - Critical time (s)
"""
function critical_time(rb::Real, k::Real, Cs::Real)
    return 5 * rb^2 / (k / Cs)
end

"""
    residence_time(v, H)
    residence_time(V, H, ri)

Calculate the residence time of the fluid in the borehole, which is the time it takes for the fluid to travel down and back up the borehole.
# Arguments
    - `V`: Volumetric flow rate (m³/s)
    - `v`: Fluid velocity in the borehole (m/s)
    - `H`: Borehole depth (m)
    - `ri`: Inner radius of the borehole (m)
# Output
    - Residence time (s)
"""
function residence_time(v::Real, H::Real)
    return 2 * H / v # Time for fluid to travel down and back up the borehole
end
function residence_time(V::Real, H::Real, ri::Real)
    A = π * ri^2 # Cross-sectional area of the borehole
    v = V / A    # Fluid velocity in the borehole
    return residence_time(v, H)
end

"""
    residence_time_indice(t, tr)
    residence_time_indice(t, v, H)
    residence_time_indice(t, V, H, ri)

Find the index in the time vector `t` that corresponds to the residence time of the fluid in
the borehole.
# Arguments
    - `t`: Time vector (s)
    - `tr`: Residence time of the fluid in the borehole (s)
    - `V`: Volumetric flow rate (m³/s)
    - `v`: Fluid velocity in the borehole (m/s)
    - `H`: Borehole depth (m)
    - `ri`: Inner radius of the borehole (m)
# Output
    - Index in the time vector `t` that corresponds to the residence time `tr`
"""
function residence_time_indice(t::AbstractVector{<:Real}, tr::Real)
    return findfirst(x -> x >= tr, t)
end
function residence_time_indice(t::AbstractVector{<:Real}, v::Real, H::Real)
    tr = residence_time(v, H)
    return residence_time_indice(t, tr)
end
function residence_time_indice(t::AbstractVector{<:Real}, V::Real, H::Real, ri::Real)
    tr = residence_time(V, H, ri)
    return residence_time_indice(t, tr)
end