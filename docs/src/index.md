# ThermalResponseTest.jl

*Interpretation of thermal response tests, in Julia.*

`ThermalResponseTest.jl` interpretes a measured **thermal response test (TRT)** (power, inlet and
outlet fluid temperature logged over time on a ground heat exchanger) to provide the ground and
borehole properties that govern its long-term thermal performance:

- the **ground thermal conductivity** ``k_s``,
- the **effective borehole thermal resistance** ``R_b^*``,
- and, with the moving ground models, the **groundwater Darcy velocity** ``v_D``.

Two complementary families of interpretation methods are provided:

1. **First-order approximation (FOA)** — fast, closed-form regressions built on the infinite line
   source, following Pasquier (2018). Four variants cover the heating and recovery phases, using
   either the temperature itself or its time derivative.
2. **Model inversion** — bounded least-squares fitting of a full ground-response model (infinite
   line/cylindrical source, finite line source, or their moving groundwater-advection variants) to
   the entire measured signal, using [Optimization.jl](https://github.com/SciML/Optimization.jl).

`ThermalResponseTest.jl` is the interpretation layer of the GHE-jl ecosystem. It depends on
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl), which re-exports the
[GroundResponse.jl](https://github.com/GHE-jl/GroundResponse.jl) ground models and the
temporal-superposition `convolution` that both the model inversions and the synthetic tests in this
package's test suite are built on.

## Installation

The GHE-jl packages are not yet registered. Clone the ecosystem packages side by side and develop
them locally:

```julia
using Pkg
Pkg.develop(path = "../GroundHeatExchanger.jl")
Pkg.develop(path = ".")
Pkg.instantiate()
```

## Quick start

```julia
using ThermalResponseTest

# Load a TRT file [Time, Power (W), T_in (°C), T_out (°C)] and split heating / recovery
trt     = load_trt_data("data/TRT_CL_2Phases.csv")
dataset = decompose_trt(trt)

H, rb, Cs = 138.0, 0.075, 2.0e6          # depth [m], radius [m], ground heat capacity [J/m³K]

# First-order approximation (ILS, Pasquier 2018)
k_H, Rbₑ, _, _ = fit_ils_foa_T(dataset, H, rb, Cs)            # heating  (UFOA-T-H)
k_R,  _,  _    = fit_ils_foa_T_recovery(dataset, H, rb)       # recovery (UFOA-T-R)
k_dH, _,  _    = fit_ils_foa_dT(dataset, 30 / 60000, H, 0.02) # derivative (CFOA-Ṫ-H)

# Model inversion (Optimization.jl + Optim.jl)
res = fit_fls(dataset, H, rb, 4.0, trt.T_mean[1], Cs)         # FLS: res.k, res.Rbₑ
```

## Manual outline

- **[Tutorial](@ref)** — loading a TRT log, running every interpretation method, step by step.
- **Interpretation theory**:
  - [Overview](@ref) — the superposition model both method families are built on, and how to
    choose between them.
  - [First-order approximation](@ref) — the four FOA regressions and their validity windows.
  - [Model inversion](@ref) — the generic bounded least-squares core and its five model wrappers.
- **[API reference](@ref)** — the complete docstring reference for every exported function.
- **[References](@ref)** — the bibliography underpinning the implementation.

## Conventions used throughout

| Symbol | Meaning | Unit |
|---|---|---|
| ``k`` / ``k_s`` | Ground thermal conductivity | W/m·K |
| ``C_s`` | Volumetric heat capacity of the ground | J/m³·K |
| ``C_f`` | Volumetric heat capacity of the groundwater (moving models) | J/m³·K |
| ``R_b^*`` (`Rbₑ`) | Effective borehole thermal resistance | m·K/W |
| ``v_D`` | Groundwater Darcy velocity (moving models) | m/s |
| ``T_0`` | Undisturbed ground temperature | °C |
| ``q`` | Heat injection rate per unit borehole length (``Q/H``) | W/m |
| ``H``, ``D``, ``r_b`` | Borehole depth, buried depth, radius | m |
| ``t_c`` | ILS critical time, ``5 r_b^2 / (k/C_s)`` | s |
| ``t_r`` | Fluid residence time in the borehole | s |
| ``\bar t`` (`t̄`) | Heating-phase duration | s |

>**Note** Time convention: 
>Every fitting function expects `t` measured **from the start of heating** (`t = 0` the instant
>the heater turns on, strictly positive, uniformly spaced). [`decompose_trt`](@ref) produces this
>automatically as `:t_rel`, rebasing away any recirculation phase logged before heating. Hence, no
>manual pre-treatment is needed when working through a `TRTDataset`.

## Ecosystem

`ThermalResponseTest.jl` sits downstream of the simulation stack, alongside
`GroundHeatExchangerSizing.jl` and `GroundSourceHeatPumpDesign.jl`:

```
GroundResponse.jl          BoreholeResistance.jl
        │                           │
        └────────────┬──────────────┘
                     ↓
          GroundHeatExchanger.jl
                     |
        ─────────────────────────────────────────────
        |                        |                   |
ThermalResponseTest.jl   GroundHeatExchangerSizing.jl   GroundSourceHeatPumpDesign.jl
```

Calling `using GroundHeatExchanger` (a dependency of this package) already re-exports the
`GroundResponse.jl` ground models and the `convolution` superposition, so both are available
without an extra `using` statement.
