# Tutorial

This tutorial walks through interpreting a thermal response test end to end: loading the raw log,
running every first-order-approximation (FOA) method, and inverting a full ground model. Every
function used here is documented in the [API reference](@ref).

## 1. Loading and decomposing a TRT

A TRT log is a CSV with four columns: time (in seconds or as a DateTime format), heating power (in watts), inlet and outlet fluid temperature (in degree Celcius).
[`load_trt_data`](@ref) reads it, computes the mean fluid temperature, and (by default)
interpolates non-uniform time steps and trims any recirculation phase logged before heating:

```julia
using ThermalResponseTest

trt = load_trt_data("data/TRT_CL_2Phases.csv")
```

`trt` is a `DataFrame` with `:elapsed_time`, `:power`, `:T_in`, `:T_out`, `:T_mean`, starting at
`t = 0` (the undisturbed ground temperature `T0`, before heating) with `power = 0`.

[`decompose_trt`](@ref) then splits it into heating and recovery phases:

```julia
dataset = decompose_trt(trt)   # TRTDataset(data, heating, recovery)
```

`dataset.heating` and `dataset.recovery` carry an extra `:t_rel` column, time rebased so that
`t_rel = 0` at the instant heating starts. Every fitting function below expects time measured from
that instant, so a recirculation phase logged beforehand needs no manual trimming when working
through a `TRTDataset`.

>**Note** Unconventional tests:
>[`TRTDataset`](@ref) assumes a single heating phase followed by a single recovery phase. A test
>with several on/off cycles doesn't fit that shape — interpret it with the vector-based methods
>below directly on `load_trt_data`'s output instead of through `decompose_trt`.

## 2. Mean fluid temperature

[`load_trt_data`](@ref) uses the arithmetic mean by default. [`mean_fluid_temperature`](@ref)
exposes the full p-linear family of Marcotte & Pasquier (2008) if a different convention is needed:

```julia
T_mean = mean_fluid_temperature(trt.T_in, trt.T_out, :pLinear)   # or :arithmetic, :logarithmic, ...
```

## 3. First-order approximation — heating phase (UFOA-T-H)

[`fit_ils_foa_T`](@ref) regresses the heating-phase temperature against `log(t)` past the ILS
critical time `tc`, iterating `tc` and the fitted conductivity together until they converge:

```julia
H, rb, Cs = 138.0, 0.075, 2.0e6          # depth [m], radius [m], ground heat capacity [J/m³K]

k_H, Rbₑ_H, reg, indices = fit_ils_foa_T(dataset, H, rb, Cs)
```

`T0`, the undisturbed ground temperature needed by the regression's intercept, is taken
automatically from `dataset.data.T_mean[1]`. `reg` is the fitted `[t, T]` curve over `indices`, handy
for overlaying on the measured signal in visualisation.

## 4. First-order approximation — recovery phase (UFOA-T-R)

[`fit_ils_foa_T_recovery`](@ref) is the recovery-phase counterpart: it regresses against
`log(t / (t - t̄))`, where `t̄` is the heating duration. It needs no `T0` (the intercept estimates it
directly) but does need at least 5 recovery samples:

```julia
k_R, reg_R, indices_R = fit_ils_foa_T_recovery(dataset, H, rb)
```

## 5. First-order approximation — temperature derivative (CFOA-Ṫ)

The derivative methods regress `log(t) + log(dT/dt)` for heating or an analogous recovery form,
restricted to a window expressed in fluid residence times `tr` (by default `[64 tr, 512 tr]`,
wider than Pasquier's original `[4 tr, 16 tr]` because on real signals the
derivative's log-log trend often hasn't settled onto its asymptotic slope that early). `dT/dt` itself
can come from two interchangeable estimators:

- [`centered_finite_difference`](@ref) — a plain three-point finite difference.
- [`bourdet_derivative`](@ref) — Bourdet et al.'s (1989) log-spaced, noise-robust estimator; the
  default inside the `TRTDataset` overloads because real TRT logs are noisier than synthetic ones.

```julia
V, ri = 30 / 60000, 0.02   # flow rate [m³/s], borehole inner radius [m]

k_dH, reg_dH, indices_dH = fit_ils_foa_dT(dataset, V, H, ri)                    # heating
k_dR, reg_dR, indices_dR = fit_ils_foa_dT_recovery(dataset, V, H, ri)           # recovery

# Force the plain finite difference instead of the default Bourdet derivative:
k_dH2, _, _ = fit_ils_foa_dT(dataset, V, H, ri; derivative_method = :centered)
```

See [First-order approximation](@ref) for the equations and validity windows behind all four
methods.

## 6. Model inversion

Model inversion fits a full ground-response model to the **entire** measured signal (heating and
recovery together) by bounded least squares, instead of restricting itself to one phase and one
closed-form regime. [`fit_ground_response`](@ref) is the generic core: it fits any
`AbstractGroundModel` (from GroundResponse.jl, re-exported by GroundHeatExchanger.jl) plus `Rbₑ` to
the measured mean fluid temperature, by superimposing the measured load with the model's g-function
through `ground_response` (Eq. 1 of Pasquier 2018). Driving the fit through `ground_response` (rather
than calling a low-level g-function directly) means a borefield (`xy` with more than one
row), a custom `AbstractGroundModel`, or `ground_response`'s `bc`/`solver`/`interp` options are all
reachable through the same keyword arguments, with no new fitting code required:

```julia
using GroundHeatExchanger: ILSModel

res = fit_ground_response(dataset.data.elapsed_time, dataset.data.T_mean, dataset.data.power ./ H,
    rb, trt.T_mean[1], p -> ILSModel(p[1], Cs), [2.5], [0.2], [7.0])
```

[`fit_ils`](@ref), [`fit_ics`](@ref), [`fit_fls`](@ref), [`fit_mils`](@ref) and [`fit_mfls`](@ref)
are thin wrappers around `fit_ground_response`, each building the right model from the unknowns
being fit. [`fit_ils`](@ref) is the simplest case (the infinite line source, same model as the FOA
methods) but fit over the whole test:

```julia
res = fit_ils(dataset, H, rb, trt.T_mean[1], Cs)
res.k, res.Rbₑ    # fitted conductivity and effective borehole resistance
```

Four more wrappers cover the other `GroundResponse.jl` models: [`fit_ics`](@ref) (infinite
cylindrical source), [`fit_fls`](@ref) (finite line source, needs `H` and buried depth `D`), and the
moving (groundwater-advection) variants [`fit_mils`](@ref) / [`fit_mfls`](@ref), which additionally
estimate the Darcy velocity `vD`:

```julia
D, Cf = 4.0, 4.2e6                         # buried depth [m], groundwater heat capacity [J/m³K]

res_fls  = fit_fls(dataset, H, rb, D, trt.T_mean[1], Cs)
res_mfls = fit_mfls(dataset, H, rb, D, trt.T_mean[1], Cs, Cf)
res_mfls.k, res_mfls.Rbₑ, res_mfls.vD
```

See [Model inversion](@ref) for the superposition equation, the optimizer/AD defaults, and why a
`maxtime`-triggered `Failure` return code does not by itself indicate a bad fit.

## 7. Variable-power signals

Some TRT rigs don't hold power perfectly constant. [`step_signal`](@ref) segments a noisy
variable-power signal into a chosen number of constant steps by k-means clustering, useful as a
pre-processing step before feeding `q` to any of the methods above:

```julia
q_steps = step_signal(trt.power, 5)   # 5 constant power steps
```

## Validation scripts

The `script/` directory contains runnable, plotted examples that exercise every method on real TRT
data. Run them from the package root:

```
julia --project=script -e 'using Pkg; Pkg.instantiate()'
julia --project=script script/script_first_order_approximation.jl
```

| Script | What it shows |
|---|---|
| `script_first_order_approximation.jl` | The four FOA methods (heating, recovery, derivative on each) on a real TRT |
| `script_inversion.jl` | Model inversion with all five ground models; measured vs fitted |
| `script_load_trt.jl` | Loading, decomposition, and plotting of a TRT data file, including recirculation trimming and non-uniform-time-step interpolation |
| `script_mean_fluid_temperature.jl` | Comparison of the mean-fluid-temperature averaging methods |
| `script_temperature_derivative.jl` | Centered finite-difference vs Bourdet derivative, against the analytical ILS derivative |
| `script_step_signal.jl` | k-means step segmentation of a variable heating power |
