# ThermalResponseTest.jl

A Julia package to interpret **thermal response tests (TRT)** performed on ground heat exchangers.
It estimates the **ground thermal conductivity** (and the effective borehole thermal resistance, and when using *moving* such as the moving finite line source, the groundwater Darcy velocity) from a measured TRT, using two complementary
families of methods:

1. **First-order approximation (FOA)** — fast, regression-based interpretation with the infinite
   line source, following Pasquier (2018): on the heating phase, the recovery phase, and the time
   derivative of the temperature.
2. **Model inversion** — least-squares fitting of any
   [GroundResponse.jl](https://github.com/GHE-jl/GroundResponse.jl) ground model (ILS, ICS,
   FLS, MILS, MFLS, or a custom `AbstractGroundModel`) with
   [Optimization.jl](https://github.com/SciML/Optimization.jl) (Optim.jl backend), with the ground
   conductivity as the main unknown.

`ThermalResponseTest.jl` is the interpretation layer of the GHE-jl ecosystem: it depends on
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl), which re-exports
the GroundResponse.jl ground models and the temporal-superposition `convolution` that the inversions
are built on.

## Quick start

```julia
using ThermalResponseTest

# Load a TRT file [Time, Power (W), T_in (°C), T_out (°C)] and split heating / recovery
trt     = load_trt_data("data/TRT_CL_2Phases.csv")
dataset = decompose_trt(trt)

H, rb, Cs = 138.0, 0.075, 2.0e6          # depth [m], radius [m], ground heat capacity [J/m³K]

# First-order approximation (ILS, Pasquier 2018)
k_H,  Rbₑ, _, _ = fit_ils_foa_T(dataset, H, rb, Cs)                 # heating  (UFOA-T-H)
k_R,  _, _      = fit_ils_foa_T_recovery(dataset, H, rb)            # recovery (UFOA-T-R)
k_dH, _,  _     = fit_ils_foa_dT(dataset, 30/6e4, H, 0.02)          # derivative (CFOA-Ṫ-H)
k_dR, _, _      = fit_ils_foa_dT_recovery(dataset, 30/6e4, H, 0.02) # derivative (CFOA-Ṫ-R)

# Model inversion (Optimization.jl + Optim.jl)
res = fit_fls(dataset, H, rb, 4.0, trt.T_mean[1], Cs)   # FLS: res.k, res.Rbₑ
```

## Data loading and pre-processing

| Function | Purpose |
|---|---|
| `load_trt_data(file; date_format, delim, header, dt_tol, trim_recirculation, threshold)` | Read a TRT CSV structured as `[Time, Power, T_in, T_out]` into a DataFrame. Adds `:elapsed_time` [s] and `:T_mean` (arithmetic mean). Interpolates non-uniform time steps and trims any pre-heating recirculation phase. |
| `decompose_trt(trt, threshold=100.0)` | Split into a `TRTDataset(data, heating, recovery)` custom structure based on a power threshold. `heating`/`recovery` also carry `:t_rel`, time rebased to `t = 0` at the start of heating, a recirculation phase logged before heating needs no manual pre-treatment. |
| `mean_fluid_temperature(T_in, T_out, method)` | *p-linear* mean of Marcotte & Pasquier (2008): `:arithmetic`, `:logarithmic`, `:geometric`, `:harmonic`, `:pLinear`, or a numeric exponent. |
| `centered_finite_difference(t, x)` | Time derivative of a signal via simple centered finite differences. |
| `bourdet_derivative(t, T, δ=0.3)` | Time derivative of a signal, robust to measurement noise (Bourdet et al. 1989). Default for the derivative FOA methods. |
| `step_signal(x, n)` | Segment a noisy variable-power signal into `n` constant steps (k-means). Usefull for shifts in heating power or flow rate signals. |
| `critical_time`, `residence_time` | Validity-window helpers for the FOA methods. |

## First-order approximation (FOA) methods

All four FOA methods use the infinite line source and are derived in Pasquier (2018). They return the
estimated thermal conductivity together with the regression curve and the indices used.

| Function | Method | Phase / quantity | Reference |
|---|---|---|---|
| `fit_ils_foa_T` | UFOA-T-H | temperature, heating | Eq. 3 |
| `fit_ils_foa_T_recovery` | UFOA-T-R | temperature, recovery | Eq. 13 |
| `fit_ils_foa_dT` | CFOA-Ṫ-H | temperature derivative, heating | Eqs. 8–10 |
| `fit_ils_foa_dT_recovery` | CFOA-Ṫ-R | temperature derivative, recovery | Eqs. 14–18 |

```julia
k, Rbₑ, reg, ind = fit_ils_foa_T(t, T, q, rb, T0, Cs)               # or fit_ils_foa_T(dataset, H, rb, Cs)
k, reg, ind       = fit_ils_foa_T_recovery(t, T, q, rb, t̄)          # or fit_ils_foa_T_recovery(dataset, H, rb)
k, reg, ind       = fit_ils_foa_dT(t, dT, q, V, H, ri)              # or fit_ils_foa_dT(dataset, V, H, ri)
k, reg, ind       = fit_ils_foa_dT_recovery(t, dT, q, t̄, V, H, ri)  # or fit_ils_foa_dT_recovery(dataset, V, H, ri)
```

For the `TRTDataset`-based structure overload, `T_0` it is taken automatically from `dataset.data.T_mean[1]` (see `load_trt_data`).

The heating-phase and recovery-phase temperature methods (UFOA-T-H, UFOA-T-R) are valid only after
the critical time `t_c = 5 rb²/α` (the recovery check offset by the heating duration `t̄`). The derivative
methods default to the window `[64 tr, 512 tr]` after the start of the respective phase (heating
or recovery), where `tr` is the fluid residence time in the ground heat exchanger (wider than Pasquier's original `[4 tr, 16 tr]`interval),
which is often too tight for the derivative's log-log trend to have settled on real (non-synthetic)
signals. Pass an explicit `tr` or `indices` to override.

## Model inversion

`fit_ground_response` is the generic core inversion: it fits any `AbstractGroundModel` (from GroundResponse.jl,
re-exported by GroundHeatExchanger.jl) plus `Rbₑ` to the measured mean fluid temperature, by
superimposing the measured load with the model's g-function through `ground_response` (Eq. 1 of
Pasquier 2018):

```
Tf(t) = T0 + q(t)·Rbₑ + Σᵢ (qᵢ − qᵢ₋₁)·g(t − tᵢ₋₁)
```

Driving the fit through `ground_response` (rather than calling the low-level
`ils`/`ics`/`fls`/`mils`/`mfls` g-functions directly) means a borefield (`xy` with more than one row),
a custom `AbstractGroundModel`, or `ground_response`'s `bc`/`solver`/`interp` options are all reachable
through the same keyword arguments, with no new fitting code required.

The optimization itself uses `Optimization.jl` with the `Optim.jl` `Fminbox(LBFGS())` backend and finite differences by default. Bounds, initial guesses, the optimizer and the AD backend are all keyword arguments.

`fit_ils`, `fit_ics`, `fit_fls`, `fit_mils` and `fit_mfls` are thin wrappers around
`fit_ground_response`, each building the right model from the unknowns being fit:

| Function | Model | Unknowns | Extra inputs |
|---|---|---|---|
| `fit_ground_response` | any `AbstractGroundModel` | user-defined + `Rbₑ` | `rb, T0, model, p0, lb, ub` |
| `fit_ils` | infinite line source | `k`, `Rbₑ` | `rb, T0, Cs` |
| `fit_ics` | infinite cylindrical source | `k`, `Rbₑ` | `rb, T0, Cs` |
| `fit_fls` | finite line source | `k`, `Rbₑ` | `rb, H, D, T0, Cs` |
| `fit_mils` | moving infinite line source | `k`, `Rbₑ`, `vD` | `rb, T0, Cs, Cf` |
| `fit_mfls` | moving finite line source | `k`, `Rbₑ`, `vD` | `rb, H, D, T0, Cs, Cf` |

```julia
using GroundHeatExchanger: ILSModel

res = fit_ground_response(t, T, q, rb, T0, p -> ILSModel(p[1], Cs), [2.5], [0.2], [7.0])
res = fit_ils(t, T, q, rb, T0, Cs)
res = fit_mfls(t, T, q, rb, H, D, T0, Cs, Cf)
res = fit_fls(dataset, H, rb, 4.0, T0, Cs)        # TRTDataset overloads
```

> **Note** The superposition forward model assumes **uniformly spaced** time steps (a standard TRT
> logging rate) and `t` measured **from the start of heating** (`t = 0` the instant the heater turns
> on, strictly positive). Passing raw vectors directly (rather than through `decompose_trt` /
> `TRTDataset`) requires trimming or rebasing any recirculation phase first.

> **Note** The moving models require `vD > 0`.

> **Note** Every fit is capped by a `maxtime` keyword (default 30 s).

## Scripts

Runnable, plotted examples live in `script/`. Run them from the package root:

```
julia --project=script -e 'using Pkg; Pkg.instantiate()'
julia --project=script script/script_first_order_approximation.jl
```

| Script | What it shows |
|---|---|
| `script_load_trt.jl` | Loading, decomposition, and plotting of a TRT data file |
| `script_first_order_approximation.jl` | The four FOA methods (heating, recovery, derivative on each) on a real TRT |
| `script_inversion.jl` | Model inversion with all five ground models. Measured vs fitted |
| `script_mean_fluid_temperature.jl` | Comparison of the mean-fluid-temperature averaging methods |
| `script_temperature_derivative.jl` | Centered finite-difference derivative vs the analytical ILS derivative |
| `script_step_signal.jl` | k-means step segmentation of a variable heating power |

## Installation

The GHE-jl packages are not yet registered. Clone the ecosystem packages side by side and
develop them locally:

```julia
using Pkg
Pkg.develop(path = "../GroundHeatExchanger.jl")
Pkg.develop(path = ".")
Pkg.instantiate()
```

## Dependencies

### Library

| Package | Used for |
|---|---|
| [GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl) | Ground models (re-exported from GroundResponse.jl) and temporal-superposition `convolution` |
| [Optimization.jl](https://github.com/SciML/Optimization.jl) / [OptimizationOptimJL.jl](https://github.com/SciML/Optimization.jl) | Bounded least-squares model inversion |
| [FiniteDiff.jl](https://github.com/JuliaDiff/FiniteDiff.jl) | Finite-difference gradients for the inversion |
| [CSV.jl](https://github.com/JuliaData/CSV.jl) / [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) | TRT data loading and handling |
| [Clustering.jl](https://github.com/JuliaStats/Clustering.jl) | k-means step segmentation (`step_signal`) |

### Scripts only

| Package | Used in |
|---|---|
| [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) | All visualisation scripts |

## References

- Pasquier, P. (2018). Interpretation of the first hours of a thermal response test using the time
  derivative of the temperature. *Applied Energy*, 213, 56–75.
  https://doi.org/10.1016/j.apenergy.2018.01.022
- Austin, W. A. (1998). *Development of an in situ system for measuring ground thermal properties.*
  M.Sc. thesis, Oklahoma State University.
- Marcotte, D., & Pasquier, P. (2008). On the estimation of thermal resistance in borehole thermal
  conductivity test. *Renewable Energy*, 33(11), 2407–2415.
  https://doi.org/10.1016/j.renene.2008.01.021
- Bourdet, D., Ayoub, J. A., & Pirard, Y. M. (1989). Use of pressure derivative in well-test
  interpretation. *SPE Formation Evaluation*, 4(2), 293–302. https://doi.org/10.2118/12777-PA
- Beier, R. A. (2020). Deconvolution and convolution methods for thermal response tests on
  borehole heat exchangers. *Geothermics*, 86, 101786.
  https://doi.org/10.1016/j.geothermics.2019.101786
