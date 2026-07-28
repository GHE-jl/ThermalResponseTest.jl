# Overview

A thermal response test injects a (roughly) constant heat rate `Q` into a ground heat exchanger and
logs the fluid temperature response over one to several days, then — often — lets the borehole
recover with circulation but no heating. Interpreting it means recovering the ground and borehole
properties that produced that response.

## The forward model

Every method in this package, FOA or inversion, is built on the same superposition of the measured
load with a ground g-function `g(t)` [°C·m/W] (Eq. 1 of Pasquier 2018):

```math
T_f(t) = T_0 + R_b^\ast \, q(t) + \sum_i \big(q_i - q_{i-1}\big)\, g(t - t_{i-1})
```

``T_0`` is the undisturbed ground temperature, ``q = Q/H`` the heat rate per unit borehole length,
and ``R_b^\ast`` the effective borehole thermal resistance (see `BoreholeResistance.jl`'s
[effective resistance](https://GHE-jl.github.io/BoreholeResistance.jl/theory/effective/) page for
its own derivation — here it is simply fit as an unknown alongside the ground properties). The sum
is the temporal superposition of every load *change* with the ground's step response, exactly what
`GroundHeatExchanger.jl`'s `convolution` computes for a full simulation.

The two method families differ in how they extract parameters from this equation:

1. **First-order approximation** linearizes `g(t)` in its long-time (large Fourier number) limit,
   where the infinite line source reduces to a `log(t)` (heating) or `log(t/(t-t̄))` (recovery)
   form. That linearization turns the fit into an ordinary least-squares regression — fast, no
   optimizer needed — but only valid once enough time has passed for the approximation to hold, and
   only for a single borehole with the ILS ground model.
2. **Model inversion** does not linearize anything: it evaluates the *exact* superposition above
   for whatever `AbstractGroundModel` is supplied (ILS, ICS, FLS, or a moving variant with
   groundwater advection) and fits it by bounded nonlinear least squares over the *entire* signal,
   heating and recovery together. It is more general — any model, any borefield — at the cost of
   needing an optimizer.

See [First-order approximation](@ref) and [Model inversion](@ref) for each family's equations.

## Choosing a method

- Reach for **FOA** for a quick estimate or a sanity check with no dependency on optimizer
  convergence — it is the standard, ubiquitous approach for conventional single-borehole tests.
- Reach for **model inversion** when the FOA's validity assumptions don't hold well (e.g. a short
  test where little data lies past the critical time), when a groundwater Darcy velocity estimate is
  needed (the moving models), or when the finite borehole depth or a multi-borehole field matters
  enough that the infinite line source is not an adequate ground model.
- Running both and comparing `k`/`Rbₑ` is a reasonable sanity check either way, since the two
  families make different simplifying assumptions and generally won't agree by construction if one
  of them is invalid for the test at hand.

## Conventions

| Symbol | Meaning | Unit |
|---|---|---|
| ``k`` | Ground thermal conductivity (the primary unknown) | W/m·K |
| ``C_s``, ``C_f`` | Volumetric heat capacity of the ground / groundwater | J/m³·K |
| ``R_b^\ast`` | Effective borehole thermal resistance | m·K/W |
| ``v_D`` | Groundwater Darcy velocity (moving models only) | m/s |
| ``g(t)`` | Ground thermal response (g-function) | °C·m/W |
| ``q`` | Heat rate per unit borehole length, ``Q/H`` | W/m |
| ``t_c`` | ILS critical time | s |
| ``t_r`` | Fluid residence time | s |

All fitting functions expect `t` measured from the start of **heating** — `t = 0` the instant the
heater turns on — strictly positive, and (for the model-inversion methods, whose forward model is an
FFT convolution) uniformly spaced. [`decompose_trt`](@ref) enforces this automatically via `:t_rel`.
