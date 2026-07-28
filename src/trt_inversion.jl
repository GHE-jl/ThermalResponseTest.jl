using Optimization, OptimizationOptimJL, FiniteDiff

# Default bounded gradient-based solver. AutoFiniteDiff keeps the special-function- and
# quadrature-heavy g-functions (Bessel, exponential integral, numerical integration) differentiable
# without requiring dual-number support, and Fminbox enforces the physical box constraints.
const _DEFAULT_OPT = Fminbox(LBFGS())
const _DEFAULT_AD = Optimization.AutoFiniteDiff()

# Wall-clock cap for `_solve_trt` (seconds), see its docstring: generous enough to absorb Julia's
# one-time JIT compilation of the optimizer/AD/g-function call graph on the first fit of a session
# and still let the fit run to its (fast, sub-second) practical convergence.
const _DEFAULT_MAXTIME = 30.0

"""
    _solve_trt(loss, p0, lb, ub, optimizer, adtype, maxtime, abstol, reltol)

Build and solve the bounded least-squares `OptimizationProblem` shared by all model inversions.
`maxtime` remains a hard wall-clock backstop; `abstol`/`reltol` (each `nothing` by default, meaning
"use the optimizer's own default") are the actual convergence criteria and are what should be tuned
to trade fit accuracy for speed. Returns the `SciMLBase` solution object.
"""
function _solve_trt(loss, p0, lb, ub, optimizer, adtype, maxtime, abstol, reltol)
    f = OptimizationFunction(loss, adtype)
    prob = OptimizationProblem(f, p0; lb = lb, ub = ub)
    return solve(prob, optimizer; maxtime, abstol, reltol)
end

"""
    fit_ground_response(t, T, q, rb, T0, build_model, p0, lb, ub;
        xy, bc, solver, interp, optimizer, adtype, maxtime, abstol, reltol)

Generic model inversion: fits any `AbstractGroundModel` (from GroundResponse.jl, re-exported by
GroundHeatExchanger.jl) plus the effective borehole resistance `Rbₑ` to experimental mean fluid
temperature, by superimposing the measured load with the model's `ground_response` g-function
(Eq. 1 of Pasquier 2018) and least-squares fitting with `Optimization.jl`. This is the shared core
behind [`fit_ils`](@ref), [`fit_ics`](@ref), [`fit_fls`](@ref), [`fit_mils`](@ref) and
[`fit_mfls`](@ref).
# Arguments
    - `t`: Time vector, from the start of heating [s], uniformly spaced
    - `T`: Measured mean fluid temperature [°C]
    - `q`: Heat injection rate per borehole length (Q/H) [W/m]
    - `rb`: Borehole radius [m]
    - `T0`: Undisturbed ground temperature [°C]
    - `build_model`: Function mapping the ground-model unknowns `p` to an `AbstractGroundModel`
        (e.g. `p -> ILSModel(p[1], Cs)`); receives `p0`/`lb`/`ub`-sized vectors, i.e. **without**
        `Rbₑ`, which `fit_ground_response` always appends as the last decision variable, with its
        initial guess fixed to 0.1 mK/W and bounds fixed to `[0.0, 0.5]`
    - `p0`, `lb`, `ub`: Initial guess / bounds for the ground-model unknowns only (not `Rbₑ`)
# Keywords
    - `xy`: Borehole coordinates (nb x 2) [m] (default: a single borehole at the origin [0.0 0.0])
    - `bc`, `solver`: Forwarded to `ground_response` (boundary condition and solver backend; ignored
        for a single borehole)
    - `interp`: Forwarded to `ground_response` (default `true`: sub-sample onto an internal ~100-node
        grid and PCHIP-interpolate, see `ground_response`; pass `false` to evaluate exactly at `t`)
    - `optimizer`: Optim.jl optimizer (default `Fminbox(LBFGS())`)
    - `adtype`: Optimization.jl AD backend (default `AutoFiniteDiff()`)
    - `maxtime`: Wall-clock cap in seconds (default 30) for the underlying `solve` — a safety
        backstop, not the convergence criterion
    - `abstol`, `reltol`: Convergence tolerances forwarded to `solve` (default `nothing`, i.e. the
        optimizer's own default); these are what should be tuned to trade fit accuracy for speed
# Output
    - `(params, Rbₑ, model, sol)`: fitted ground-model unknowns, effective borehole resistance, the
        resulting `AbstractGroundModel`, and the optimization solution.
"""
function fit_ground_response(t::AbstractVector{<:Real}, T::AbstractVector{<:Real},
    q::AbstractVector{<:Real}, rb::Real, T0::Real, model::Function,
    p0::AbstractVector{<:Real}, lb::AbstractVector{<:Real}, ub::AbstractVector{<:Real};
    xy::AbstractMatrix{<:Real} = [0.0 0.0], bc::Symbol = :II, solver::Symbol = :successive,
    interp::Bool = true, optimizer = _DEFAULT_OPT, adtype = _DEFAULT_AD,
    maxtime::Real = _DEFAULT_MAXTIME, abstol::Union{Real, Nothing} = nothing,
    reltol::Union{Real, Nothing} = nothing)

    # Basic parameters
    tF = collect(Float64, t)
    qF = collect(Float64, q)
    TF = collect(Float64, T)
    n = length(p0)
    Rbₑ0 = 0.1

    # Loss function
    function loss(p, _)
        m = model(p[1:n])
        Rbₑ = p[n + 1]
        g = ground_response(tF, rb, xy, m; bc, solver, interp)
        return sum(abs2, fluid_temperature(tF, qF, g, T0, Rbₑ) .- TF)
    end

    # Solution
    sol = _solve_trt(loss, [collect(Float64, p0); Rbₑ0], [collect(Float64, lb); 0.0],
        [collect(Float64, ub); 0.5], optimizer, adtype, maxtime, abstol, reltol)
    params = sol.u[1:n]
    return (params = params, Rbₑ = sol.u[n + 1], model = model(params), sol = sol)
end

"""
    fit_ils(t, T, q, rb, T0, Cs; k0, lb, ub, kwargs...)
    fit_ils(dataset::TRTDataset, H, rb, T0, Cs; kwargs...)

Invert a thermal response test with the infinite line source (ILS) model, estimating the ground
thermal conductivity `k` and the effective borehole thermal resistance `Rbₑ`. Thin wrapper around
[`fit_ground_response`](@ref) with an `ILSModel`; see its docstring for the shared arguments/keywords
(`xy`, `bc`, `solver`, `interp`, `optimizer`, `adtype`, `maxtime`, `abstol`, `reltol`).

**Important**: `t` must be measured from the start of the heating phase (`t = 0` the instant the
heater turns on) and be strictly positive and uniformly spaced. If the raw log includes a
recirculation phase before heating, either trim it first or use the `TRTDataset`-based method below,
which rebases time automatically via [`decompose_trt`](@ref).
# Arguments
    - `t`, `T`, `q`, `rb`, `T0`: see [`fit_ground_response`](@ref)
    - `Cs`: Volumetric heat capacity of the ground [J/m³K]
    - `dataset`: `TRTDataset` from `decompose_trt` — `t`/`T`/`q` are taken from `dataset.data`, the
        full test (heating *and* recovery, already rebased/trimmed by `load_trt_data`); since the
        superposition forward model handles any time-varying load, this uses the whole test instead
        of restricting the fit to the heating phase alone, unlike the first-order-approximation
        methods which need heating/recovery kept separate
    - `H`: Borehole depth [m] (used to convert total power to power per unit length)
# Keywords
    - `k0`, `lb`, `ub`: Initial guess / bounds for `k` [W/mK] (default `2.5`, `0.2`, `7.0`)
# Output
    - `(k, Rbₑ, sol)`: estimated conductivity, effective borehole resistance, and the optimization
        solution
"""
function fit_ils(t::AbstractVector{<:Real}, T::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    rb::Real, T0::Real, Cs::Real; k0::Real = 2.5, lb::Real = 0.2, ub::Real = 7.0, kwargs...)
    r = fit_ground_response(t, T, q, rb, T0, p -> ILSModel(p[1], Cs), [k0], [lb], [ub]; kwargs...)
    return (k = r.params[1], Rbₑ = r.Rbₑ, sol = r.sol)
end
function fit_ils(dataset::TRTDataset, H::Real, rb::Real, T0::Real, Cs::Real; kwargs...)
    data = dataset.data
    return fit_ils(data.elapsed_time, data.T_mean, data.power ./ H, rb, T0, Cs; kwargs...)
end

"""
    fit_ics(t, T, q, rb, T0, Cs; k0, lb, ub, kwargs...)
    fit_ics(dataset::TRTDataset, H, rb, T0, Cs; kwargs...)

Invert a thermal response test with the infinite cylindrical source (ICS) model (cylinder radius
taken as the borehole radius `rb`), estimating `k` and `Rbₑ`. Thin wrapper around
[`fit_ground_response`](@ref) with an `ICSModel`. See [`fit_ils`](@ref) for the argument and keyword
conventions.
# Output
    - `(k, Rbₑ, sol)`
"""
function fit_ics(t::AbstractVector{<:Real}, T::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    rb::Real, T0::Real, Cs::Real; k0::Real = 2.5, lb::Real = 0.2, ub::Real = 7.0, kwargs...)
    r = fit_ground_response(t, T, q, rb, T0, p -> ICSModel(rb, p[1], Cs), [k0], [lb], [ub]; kwargs...)
    return (k = r.params[1], Rbₑ = r.Rbₑ, sol = r.sol)
end
function fit_ics(dataset::TRTDataset, H::Real, rb::Real, T0::Real, Cs::Real; kwargs...)
    data = dataset.data
    return fit_ics(data.elapsed_time, data.T_mean, data.power ./ H, rb, T0, Cs; kwargs...)
end

"""
    fit_fls(t, T, q, rb, H, D, T0, Cs; k0, lb, ub, kwargs...)
    fit_fls(dataset::TRTDataset, H, rb, D, T0, Cs; kwargs...)

Invert a thermal response test with the finite line source (FLS) model, estimating `k` and `Rbₑ`.
The FLS additionally needs the borehole depth `H` and buried depth `D`. Thin wrapper around
[`fit_ground_response`](@ref) with an `FLSModel`. See [`fit_ils`](@ref) for the remaining argument
and keyword conventions.
# Output
    - `(k, Rbₑ, sol)`
"""
function fit_fls(t::AbstractVector{<:Real}, T::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    rb::Real, H::Real, D::Real, T0::Real, Cs::Real;
    k0::Real = 2.5, lb::Real = 0.2, ub::Real = 7.0, kwargs...)
    r = fit_ground_response(t, T, q, rb, T0, p -> FLSModel(H, D, p[1], Cs), [k0], [lb], [ub]; kwargs...)
    return (k = r.params[1], Rbₑ = r.Rbₑ, sol = r.sol)
end
function fit_fls(dataset::TRTDataset, H::Real, rb::Real, D::Real, T0::Real, Cs::Real; kwargs...)
    data = dataset.data
    return fit_fls(data.elapsed_time, data.T_mean, data.power ./ H, rb, H, D, T0, Cs; kwargs...)
end

"""
    fit_mils(t, T, q, rb, T0, Cs, Cf; k0, vD0, lb, ub, kwargs...)
    fit_mils(dataset::TRTDataset, H, rb, T0, Cs, Cf; kwargs...)

Invert a thermal response test with the moving infinite line source (MILS) model, estimating the
ground thermal conductivity `k`, the effective borehole resistance `Rbₑ` **and the Darcy velocity
`vD`** of the groundwater flow. The model is evaluated at the borehole wall (circumferential-average
branch, like `ground_response`). Requires the groundwater volumetric heat capacity `Cf`. Thin wrapper
around [`fit_ground_response`](@ref) with a `MILSModel`. See [`fit_ils`](@ref) for the shared
argument/keyword conventions.
# Keywords
    - `k0`, `vD0`: Initial guesses for `[k, vD]` (default `2.5`, `1e-7` m/s)
    - `lb`, `ub`: bounds for `[k, vD]` (default `[0.2, 1e-9]`, `[7.0, 1e-5]`)
# Output
    - `(k, Rbₑ, vD, sol)`
"""
function fit_mils(t::AbstractVector{<:Real}, T::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    rb::Real, T0::Real, Cs::Real, Cf::Real;
    k0::Real = 2.5, vD0::Real = 1e-7, lb = [0.2, 1e-9], ub = [7.0, 1e-5], kwargs...)
    build = p -> MILSModel(rb, p[1], Cs, Cf, p[2])
    r = fit_ground_response(t, T, q, rb, T0, build, [k0, vD0], lb, ub; kwargs...)
    return (k = r.params[1], Rbₑ = r.Rbₑ, vD = r.params[2], sol = r.sol)
end
function fit_mils(dataset::TRTDataset, H::Real, rb::Real, T0::Real, Cs::Real, Cf::Real; kwargs...)
    data = dataset.data
    return fit_mils(data.elapsed_time, data.T_mean, data.power ./ H, rb, T0, Cs, Cf; kwargs...)
end

"""
    fit_mfls(t, T, q, rb, H, D, T0, Cs, Cf; k0, vD0, lb, ub, kwargs...)
    fit_mfls(dataset::TRTDataset, H, rb, D, T0, Cs, Cf; kwargs...)

Invert a thermal response test with the moving finite line source (MFLS) model, estimating `k`,
`Rbₑ` **and the Darcy velocity `vD`**. Combines the finite-depth geometry (`H`, `D`) with groundwater
advection (`Cf`, `vD`). Thin wrapper around [`fit_ground_response`](@ref) with a `MFLSModel`. See
[`fit_mils`](@ref) for the keyword conventions.
# Output
    - `(k, Rbₑ, vD, sol)`
"""
function fit_mfls(t::AbstractVector{<:Real}, T::AbstractVector{<:Real}, q::AbstractVector{<:Real},
    rb::Real, H::Real, D::Real, T0::Real, Cs::Real, Cf::Real;
    k0::Real = 2.5, vD0::Real = 1e-7, lb = [0.2, 1e-9], ub = [7.0, 1e-5], kwargs...)
    build = p -> MFLSModel(H, rb, D, p[1], Cs, Cf, p[2])
    r = fit_ground_response(t, T, q, rb, T0, build, [k0, vD0], lb, ub; kwargs...)
    return (k = r.params[1], Rbₑ = r.Rbₑ, vD = r.params[2], sol = r.sol)
end
function fit_mfls(dataset::TRTDataset, H::Real, rb::Real, D::Real, T0::Real, Cs::Real, Cf::Real;
    kwargs...)
    data = dataset.data
    return fit_mfls(data.elapsed_time, data.T_mean, data.power ./ H, rb, H, D, T0, Cs, Cf; kwargs...)
end
