using Optimization, OptimizationOptimJL, FiniteDiff
using DataFrames: AbstractDataFrame

# =====================================================================================
# Model inversion of a thermal response test
#
# Each `fit_*` function estimates ground/borehole parameters by least-squares fitting one of the
# GroundResponse.jl ground models (re-exported through GroundHeatExchanger.jl) to the measured mean
# fluid temperature. The ground thermal conductivity `k` is always the main unknown; the moving
# models (MILS, MFLS) additionally recover the Darcy velocity `vD`.
#
# Forward model (single borehole, Eq. 1 of Pasquier 2018):
#
#     Tf(t) = T₀ + Rb·q(t) + Σᵢ (qᵢ - qᵢ₋₁) · g(t - tᵢ₋₁)
#           = T₀ .+ Rb .* q .+ convolution(q, g)
#
# where `g` is the model g-function [°Cm/W] and `convolution` (from GroundHeatExchanger.jl) performs
# the temporal superposition by FFT. The convolution assumes **uniformly spaced** time steps.
# =====================================================================================

"""
    _trt_forward(g, q, Rb, T0)

Mean fluid temperature predicted by the superposition forward model from a g-function vector `g`
[°Cm/W], the load vector `q` [W/m], the borehole resistance `Rb` [mK/W] and the undisturbed ground
temperature `T0` [°C]. Internal helper shared by all `fit_*` functions.
"""
_trt_forward(g, q, Rb, T0) = T0 .+ Rb .* q .+ convolution(q, g)

"""
    _solve_trt(loss, p0, lb, ub, optimizer, adtype)

Build and solve the bounded least-squares `OptimizationProblem` shared by all model inversions.
Returns the `SciMLBase` solution object.
"""
function _solve_trt(loss, p0, lb, ub, optimizer, adtype)
    f = OptimizationFunction(loss, adtype)
    prob = OptimizationProblem(f, p0; lb = lb, ub = ub)
    return solve(prob, optimizer)
end

# Default bounded gradient-based solver. AutoFiniteDiff keeps the special-function- and
# quadrature-heavy g-functions (Bessel, exponential integral, numerical integration) differentiable
# without requiring dual-number support, and Fminbox enforces the physical box constraints.
const _DEFAULT_OPT = Fminbox(LBFGS())
const _DEFAULT_AD = Optimization.AutoFiniteDiff()

# -------------------------------------------------------------------------------------
# Infinite line source (ILS) — unknowns [k, Rb]
# -------------------------------------------------------------------------------------
"""
    fit_ils(t, T, q, rb, T0, Cs; k0, Rb0, lb, ub, optimizer, adtype)
    fit_ils(trt, H, rb, T0, Cs; kwargs...)

Invert a thermal response test with the infinite line source (ILS) model, estimating the ground
thermal conductivity `k` and the borehole thermal resistance `Rb`. The mean fluid temperature is
reconstructed by temporal superposition of the measured load with the ILS g-function and fitted to
the measurements with `Optimization.jl` (Optim.jl `Fminbox(LBFGS())` backend by default).
# Arguments
    - `t`: Time vector, uniformly spaced [s]
    - `T`: Measured mean fluid temperature [°C]
    - `q`: Heat injection rate per borehole length (Q/H) [W/m]
    - `rb`: Borehole radius [m]
    - `T0`: Undisturbed ground temperature [°C]
    - `Cs`: Volumetric heat capacity of the ground [J/m³K]
    - `trt`: DataFrame with columns `:elapsed_time`, `:T_mean`, `:power`
    - `H`: Borehole depth [m] (used to convert total power to power per unit length)
# Keywords
    - `k0`, `Rb0`: Initial guesses for `k` [W/mK] and `Rb` [mK/W]
    - `lb`, `ub`: Lower/upper bounds for `[k, Rb]`
    - `optimizer`: Optim.jl optimizer (default `Fminbox(LBFGS())`)
    - `adtype`: Optimization.jl AD backend (default `AutoFiniteDiff()`)
# Output
    - `(k, Rb, sol)`: estimated conductivity, borehole resistance, and the optimization solution
"""
function fit_ils(t::AbstractVector{<:Real}, T::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    rb::Real, T0::Real, Cs::Real;
    k0::Real = 2.5, Rb0::Real = 0.1, lb = [0.2, 0.0], ub = [7.0, 0.5],
    optimizer = _DEFAULT_OPT, adtype = _DEFAULT_AD)
    tF = collect(Float64, t); qF = collect(Float64, q); TF = collect(Float64, T)
    function loss(p, _)
        k, Rb = p
        g = ils(tF, rb, k, Cs)
        return sum(abs2, _trt_forward(g, qF, Rb, T0) .- TF)
    end
    sol = _solve_trt(loss, [float(k0), float(Rb0)], lb, ub, optimizer, adtype)
    return (k = sol.u[1], Rb = sol.u[2], sol = sol)
end
function fit_ils(trt::AbstractDataFrame, H::Real, rb::Real, T0::Real, Cs::Real; kwargs...)
    return fit_ils(trt.elapsed_time, trt.T_mean, trt.power ./ H, rb, T0, Cs; kwargs...)
end

# -------------------------------------------------------------------------------------
# Infinite cylindrical source (ICS) — unknowns [k, Rb]
# -------------------------------------------------------------------------------------
"""
    fit_ics(t, T, q, rb, T0, Cs; kwargs...)
    fit_ics(trt, H, rb, T0, Cs; kwargs...)

Invert a thermal response test with the infinite cylindrical source (ICS) model (cylinder radius
taken as the borehole radius `rb`), estimating `k` and `Rb`. See [`fit_ils`](@ref) for the argument
and keyword conventions.
# Output
    - `(k, Rb, sol)`
"""
function fit_ics(t::AbstractVector{<:Real}, T::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    rb::Real, T0::Real, Cs::Real;
    k0::Real = 2.5, Rb0::Real = 0.1, lb = [0.2, 0.0], ub = [7.0, 0.5],
    optimizer = _DEFAULT_OPT, adtype = _DEFAULT_AD)
    tF = collect(Float64, t); qF = collect(Float64, q); TF = collect(Float64, T)
    function loss(p, _)
        k, Rb = p
        g = ics(tF, rb, rb, k, Cs)
        return sum(abs2, _trt_forward(g, qF, Rb, T0) .- TF)
    end
    sol = _solve_trt(loss, [float(k0), float(Rb0)], lb, ub, optimizer, adtype)
    return (k = sol.u[1], Rb = sol.u[2], sol = sol)
end
function fit_ics(trt::AbstractDataFrame, H::Real, rb::Real, T0::Real, Cs::Real; kwargs...)
    return fit_ics(trt.elapsed_time, trt.T_mean, trt.power ./ H, rb, T0, Cs; kwargs...)
end

# -------------------------------------------------------------------------------------
# Finite line source (FLS) — unknowns [k, Rb]
# -------------------------------------------------------------------------------------
"""
    fit_fls(t, T, q, rb, H, D, T0, Cs; kwargs...)
    fit_fls(trt, H, rb, D, T0, Cs; kwargs...)

Invert a thermal response test with the finite line source (FLS) model, estimating `k` and `Rb`.
The FLS additionally needs the borehole depth `H` and buried depth `D`. See [`fit_ils`](@ref) for
the remaining argument and keyword conventions.
# Output
    - `(k, Rb, sol)`
"""
function fit_fls(t::AbstractVector{<:Real}, T::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    rb::Real, H::Real, D::Real, T0::Real, Cs::Real;
    k0::Real = 2.5, Rb0::Real = 0.1, lb = [0.2, 0.0], ub = [7.0, 0.5],
    optimizer = _DEFAULT_OPT, adtype = _DEFAULT_AD)
    tF = collect(Float64, t); qF = collect(Float64, q); TF = collect(Float64, T)
    function loss(p, _)
        k, Rb = p
        g = fls(tF, rb, H, D, k, Cs)
        return sum(abs2, _trt_forward(g, qF, Rb, T0) .- TF)
    end
    sol = _solve_trt(loss, [float(k0), float(Rb0)], lb, ub, optimizer, adtype)
    return (k = sol.u[1], Rb = sol.u[2], sol = sol)
end
function fit_fls(trt::AbstractDataFrame, H::Real, rb::Real, D::Real, T0::Real, Cs::Real; kwargs...)
    return fit_fls(trt.elapsed_time, trt.T_mean, trt.power ./ H, rb, H, D, T0, Cs; kwargs...)
end

# -------------------------------------------------------------------------------------
# Moving infinite line source (MILS) — unknowns [k, Rb, vD]
# -------------------------------------------------------------------------------------
"""
    fit_mils(t, T, q, rb, T0, Cs, Cf; kwargs...)
    fit_mils(trt, H, rb, T0, Cs, Cf; kwargs...)

Invert a thermal response test with the moving infinite line source (MILS) model, estimating the
ground thermal conductivity `k`, the borehole resistance `Rb` **and the Darcy velocity `vD`** of the
groundwater flow. The model is evaluated at the borehole wall (circumferential-average branch, like
`GroundResponse.ground_response`). Requires the groundwater volumetric heat capacity `Cf`.
# Keywords
    - `k0`, `Rb0`, `vD0`: initial guesses (`vD0` defaults to `1e-7` m/s)
    - `lb`, `ub`: bounds for `[k, Rb, vD]`
    - `optimizer`, `adtype`: as in [`fit_ils`](@ref)
# Output
    - `(k, Rb, vD, sol)`
"""
function fit_mils(t::AbstractVector{<:Real}, T::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    rb::Real, T0::Real, Cs::Real, Cf::Real;
    k0::Real = 2.5, Rb0::Real = 0.1, vD0::Real = 1e-7,
    lb = [0.2, 0.0, 1e-9], ub = [7.0, 0.5, 1e-5],
    optimizer = _DEFAULT_OPT, adtype = _DEFAULT_AD)
    tF = collect(Float64, t); qF = collect(Float64, q); TF = collect(Float64, T)
    xy = [0.0, 0.0]
    function loss(p, _)
        k, Rb, vD = p
        g = mils(tF, xy, rb, k, Cs, Cf, vD)
        return sum(abs2, _trt_forward(g, qF, Rb, T0) .- TF)
    end
    sol = _solve_trt(loss, [float(k0), float(Rb0), float(vD0)], lb, ub, optimizer, adtype)
    return (k = sol.u[1], Rb = sol.u[2], vD = sol.u[3], sol = sol)
end
function fit_mils(trt::AbstractDataFrame, H::Real, rb::Real, T0::Real, Cs::Real, Cf::Real; kwargs...)
    return fit_mils(trt.elapsed_time, trt.T_mean, trt.power ./ H, rb, T0, Cs, Cf; kwargs...)
end

# -------------------------------------------------------------------------------------
# Moving finite line source (MFLS) — unknowns [k, Rb, vD]
# -------------------------------------------------------------------------------------
"""
    fit_mfls(t, T, q, rb, H, D, T0, Cs, Cf; kwargs...)
    fit_mfls(trt, H, rb, D, T0, Cs, Cf; kwargs...)

Invert a thermal response test with the moving finite line source (MFLS) model, estimating `k`,
`Rb` **and the Darcy velocity `vD`**. Combines the finite-depth geometry (`H`, `D`) with groundwater
advection (`Cf`, `vD`). See [`fit_mils`](@ref) for the keyword conventions.
# Output
    - `(k, Rb, vD, sol)`
"""
function fit_mfls(t::AbstractVector{<:Real}, T::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    rb::Real, H::Real, D::Real, T0::Real, Cs::Real, Cf::Real;
    k0::Real = 2.5, Rb0::Real = 0.1, vD0::Real = 1e-7,
    lb = [0.2, 0.0, 1e-9], ub = [7.0, 0.5, 1e-5],
    optimizer = _DEFAULT_OPT, adtype = _DEFAULT_AD)
    tF = collect(Float64, t); qF = collect(Float64, q); TF = collect(Float64, T)
    xy = [0.0, 0.0]
    function loss(p, _)
        k, Rb, vD = p
        g = mfls(tF, xy, H, rb, D, k, Cs, Cf, vD)
        return sum(abs2, _trt_forward(g, qF, Rb, T0) .- TF)
    end
    sol = _solve_trt(loss, [float(k0), float(Rb0), float(vD0)], lb, ub, optimizer, adtype)
    return (k = sol.u[1], Rb = sol.u[2], vD = sol.u[3], sol = sol)
end
function fit_mfls(trt::AbstractDataFrame, H::Real, rb::Real, D::Real, T0::Real, Cs::Real, Cf::Real;
    kwargs...)
    return fit_mfls(trt.elapsed_time, trt.T_mean, trt.power ./ H, rb, H, D, T0, Cs, Cf; kwargs...)
end
