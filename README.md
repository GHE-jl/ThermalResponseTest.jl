# ThermalResponseTest.jl

A Julia package to interpret **thermal response tests (TRT)** performed on ground heat exchangers.
It estimates the **ground thermal conductivity** (and the borehole thermal resistance, and — with
the moving models — the groundwater Darcy velocity) from a measured TRT, using two complementary
families of methods:

1. **First-order approximation (FOA)** — fast, regression-based interpretation with the infinite
   line source, following Pasquier (2018): on the heating phase, the recovery phase, and the time
   derivative of the heating phase.
2. **Model inversion** — least-squares fitting of any
   [GroundResponse.jl](https://github.com/GeothermalJL/GroundResponse.jl) ground model (ILS, ICS,
   FLS, MILS, MFLS) with [Optimization.jl](https://github.com/SciML/Optimization.jl) (Optim.jl
   backend), with the ground conductivity as the main unknown.

`ThermalResponseTest.jl` is the interpretation layer of the GeothermalJL ecosystem: it depends on
[GroundHeatExchanger.jl](https://github.com/GeothermalJL/GroundHeatExchanger.jl), which re-exports
the GroundResponse.jl ground models and the temporal-superposition `convolution` that the inversions
are built on.

## Quick start

```julia
using ThermalResponseTest

# Load a TRT file [Time, Power (W), T_in (°C), T_out (°C)] and split heating / recovery
trt     = load_trt_data("trt_data/DataCL_TRT.csv")
dataset = decompose_trt(trt)

H, rb, Cs = 138.0, 0.075, 2.0e6          # depth [m], radius [m], ground heat capacity [J/m³K]

# --- First-order approximation (ILS, Pasquier 2018) ---
k_H,  Rb, _, _ = fit_ils_foa_T(dataset.heating, H, rb, Cs)   # heating  (UFOA-T-H)
k_R,  T0, _, _ = fit_ils_foa_T_recovery(dataset, H)          # recovery (UFOA-T-R)
k_dH, _,  _    = fit_ils_foa_dT(dataset.heating, 30/60000, H, 0.02)  # derivative (CFOA-Ṫ-H)

# --- Model inversion (Optimization.jl + Optim.jl) ---
res = fit_fls(dataset.heating, H, rb, 4.0, trt.T_mean[1], Cs)   # FLS: res.k, res.Rb
```

## Data loading and pre-processing

| Function | Purpose |
|---|---|
| `load_trt_data(file; date_format, delim, header)` | Read a TRT CSV `[Time, Power, T_in, T_out]` into a DataFrame; adds `:elapsed_time` [s] and `:T_mean` (p-linear mean) |
| `decompose_trt(trt, threshold=100.0)` | Split into a `TRTDataset(full_data, heating, cooling)` on a power threshold |
| `mean_fluid_temperature(T_in, T_out, method)` | p-linear mean of Marcotte & Pasquier (2008): `:arithmetic`, `:logarithmic`, `:geometric`, `:harmonic`, `:pLinear`, or a numeric exponent |
| `centered_finite_difference(t, x)` | Time derivative of a signal (for the derivative methods) |
| `step_signal(x, n)` | Segment a noisy variable-power signal into `n` constant steps (k-means) |
| `critical_time`, `residence_time`, `residence_time_indice` | Validity-window helpers for the FOA methods |

## First-order approximation methods

All three use the infinite line source and are derived in Pasquier (2018). They return the
estimated thermal conductivity together with the regression curve and the indices used.

| Function | Method | Phase / quantity | Reference |
|---|---|---|---|
| `fit_ils_foa_T` | UFOA-T-H | temperature, heating | Eq. 3 |
| `fit_ils_foa_T_recovery` | UFOA-T-R | temperature, recovery | Eq. 13 |
| `fit_ils_foa_dT` | CFOA-Ṫ-H | temperature derivative, heating | Eqs. 8–10 |

```julia
k, Rb, reg, ind    = fit_ils_foa_T(t, T, q, rb, Cs)          # or fit_ils_foa_T(heating_df, H, rb, Cs)
k, T0, reg, ind    = fit_ils_foa_T_recovery(t, T, q, t̄)      # or fit_ils_foa_T_recovery(dataset, H)
k, reg, ind        = fit_ils_foa_dT(t, dT, q, V, H, ri)      # or fit_ils_foa_dT(heating_df, V, H, ri)
```

The heating-phase methods are valid only after the critical time `tc = 5 rb²/α`; the derivative
method is constrained to roughly `[4 tr, 16 tr]` where `tr` is the fluid residence time.

## Model inversion

Each `fit_*` function fits a ground model to the measured mean fluid temperature by superimposing
the measured load with the model g-function (Eq. 1 of Pasquier 2018):

```
Tf(t) = T0 + Rb·q(t) + Σᵢ (qᵢ − qᵢ₋₁)·g(t − tᵢ₋₁)
```

The optimization uses `Optimization.jl` with the `Optim.jl` `Fminbox(LBFGS())` backend and finite
differences by default; bounds, initial guesses, the optimizer and the AD backend are all keyword
arguments.

| Function | Model | Unknowns | Extra inputs |
|---|---|---|---|
| `fit_ils` | infinite line source | `k`, `Rb` | `rb, T0, Cs` |
| `fit_ics` | infinite cylindrical source | `k`, `Rb` | `rb, T0, Cs` |
| `fit_fls` | finite line source | `k`, `Rb` | `rb, H, D, T0, Cs` |
| `fit_mils` | moving infinite line source | `k`, `Rb`, `vD` | `rb, T0, Cs, Cf` |
| `fit_mfls` | moving finite line source | `k`, `Rb`, `vD` | `rb, H, D, T0, Cs, Cf` |

```julia
res = fit_ils(t, T, q, rb, T0, Cs)                 # res.k, res.Rb, res.sol
res = fit_mfls(t, T, q, rb, H, D, T0, Cs, Cf)       # res.k, res.Rb, res.vD, res.sol

# DataFrame overloads take the (heating) TRT data and a depth H for q = power/H
res = fit_fls(dataset.heating, H, rb, 4.0, T0, Cs)
```

> **Note** The superposition forward model assumes **uniformly spaced** time steps (a standard TRT
> logging rate). The moving models require `vD > 0`.

## Scripts

Runnable, plotted examples live in `script/`. Run them from the package root:

```
julia --project=script -e 'using Pkg; Pkg.instantiate()'
julia --project=script script/script_first_order_approximation.jl
```

Each script saves its figure to `script/figures/`.

| Script | What it shows |
|---|---|
| `script_first_order_approximation.jl` | The three FOA methods (heating, recovery, derivative) on a real TRT |
| `script_inversion.jl` | Model inversion with all five ground models; measured vs fitted |
| `script_load_trt.jl` | Loading, decomposition, and plotting of a TRT data file |
| `script_mean_fluid_temperature.jl` | Comparison of the mean-fluid-temperature averaging methods |
| `script_temperature_derivative.jl` | Centered finite-difference derivative vs the analytical ILS derivative |
| `script_step_signal.jl` | k-means step segmentation of a variable heating power |

## Installation

The GeothermalJL packages are not yet registered. Clone the ecosystem packages side by side and
develop them locally:

```julia
using Pkg
Pkg.develop(path = "../BoreholeResistance.jl")
Pkg.develop(path = "../GroundResponse.jl")
Pkg.develop(path = "../GroundHeatExchanger.jl")
Pkg.develop(path = ".")
Pkg.instantiate()
```

## Dependencies

### Library

| Package | Used for |
|---|---|
| [GroundHeatExchanger.jl](https://github.com/GeothermalJL/GroundHeatExchanger.jl) | Ground models (re-exported from GroundResponse.jl) and temporal-superposition `convolution` |
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
