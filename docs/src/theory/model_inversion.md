# Model inversion

Where the first-order approximation linearizes the ILS ground model, model inversion fits multiple
 models exactly, at the cost of a nonlinear optimization instead of a linear regression.

## The generic core

[`fit_ground_response`](@ref) is the shared engine behind every `fit_*` wrapper on this page. Given
a `model` function mapping the ground-model unknowns `p` to an `AbstractGroundModel`, it
minimizes the sum of squared residuals between the measured mean fluid temperature and the
superposition forward model from [Overview](@ref):

```math
\min_{p,\, R_b^\ast} \; \sum_i \Big[\big(T_0 +  q(t_i)R_b^\ast + \left(q \ast g_p\right)(t_i)\big) - T_{\text{meas}}(t_i)\Big]^2
```

`Rbₑ` is always appended as the last decision variable. The model-specific unknowns `p` (e.g. `k`
for the ILS, `[k, vD]` for the moving models) come from `model`. Guessed parameters have initial values `p0`, and bounds `lb` and `ub`
**without** `Rbₑ`, whose own initial guess (0.1 m·K/W) and bounds (`[0.0, 0.5]`) are fixed
internally. #TODO retype this paragraph. It is unclear.

The convolution `q ∗ gₚ` is evaluated through `GroundHeatExchanger.jl`'s fluid_response() (#TODO confirm, the convolution is in GroundHeatExchanger.jl, not in GroundResponse.jl). So, a borefield (`xy` with more than one row), a
custom `AbstractGroundModel`, or any of `ground_response`'s `bc`/`solver`/`interp` options are all
reachable through `fit_ground_response`'s keywords without new fitting code.

## Optimizer and AD backend

The bounded problem is solved with [Optimization.jl](https://github.com/SciML/Optimization.jl)
using the `Optim.jl` `Fminbox(LBFGS())` backend by default. `Fminbox` enforces the physical box
constraints (conductivity, Darcy velocity and `Rbₑ` all bounded away from non-physical values), and
`LBFGS` is a solid default for a smooth, low-dimensional (2–3 unknowns) least-squares problem.
Gradients use `AutoFiniteDiff()` rather than automatic differentiation, since the ground g-functions
involve special functions and numerical quadrature (Bessel functions, the exponential integral, or
adaptive integration for the finite line source) that are not guaranteed differentiable with a
forward-mode dual number.

`maxtime` (default 30 s) is a **wall-clock backstop**, not the convergence criterion. It is generous
enough to absorb Julia's JIT compilation of the optimizer/AD/g-function call graph on the first fit
of a session. `abstol`/`reltol` (default `nothing`, the optimizer's own default) are the actual
tolerances to tune if a fit needs to trade accuracy for speed.

>**Note** A `maxtime`-triggered `Failure` is expected, not a bad sign"
>`Fminbox`'s outer barrier-method loop essentially never reports `Success` for this problem once
>the finite-difference gradient becomes noise-dominated near the optimum (even though the
>estimate itself has already stabilized, typically in well under a second. `sol.retcode ==
>Failure` from hitting `maxtime` does not by itself indicate the fit is wrong. Check the fitted
>curve against the data instead. Raise `maxtime` only for a harder problem (e.g. a multi-borehole
>field) where the estimate genuinely hasn't stabilized yet.

## The five model wrappers

Each wrapper below is a thin, few-line adapter around [`fit_ground_response`](@ref): it builds the
right `AbstractGroundModel` from `p` and supplies that model's own required geometry.

| Function | Ground model | Unknowns | Extra geometry needed |
|---|---|---|---|
| [`fit_ils`](@ref) | Infinite line source | `k`, `Rbₑ` | — |
| [`fit_ics`](@ref) | Infinite cylindrical source (radius `rb`) | `k`, `Rbₑ` | — |
| [`fit_fls`](@ref) | Finite line source | `k`, `Rbₑ` | `H`, `D` (buried depth) |
| [`fit_mils`](@ref) | Moving infinite line source | `k`, `Rbₑ`, `vD` | `Cf` (groundwater heat capacity) |
| [`fit_mfls`](@ref) | Moving finite line source | `k`, `Rbₑ`, `vD` | `H`, `D`, `Cf` |

The moving variants add the groundwater Darcy velocity `vD` as a second unknown. Since the search
space grows and the two-parameter interaction is less well-conditioned, expect these fits to need
more of the `maxtime` budget than the single-unknown models.

`TRTDataset`-based overloads of every wrapper fit the **whole test** (`dataset.data`: heating and
recovery together), unlike the first-order-approximation methods, which need the phases kept
separate.

## Functions on this page

```@docs
fit_ground_response
fit_ils
fit_ics
fit_fls
fit_mils
fit_mfls
```
