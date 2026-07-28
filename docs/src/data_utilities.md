# Data & utilities

Before any interpretation method runs, a raw TRT log needs to become a clean, uniformly-sampled,
correctly-phased signal. This page covers that pipeline and the small standalone utilities the
interpretation methods build on.

## Loading a TRT log

[`load_trt_data`](@ref) reads a four-column CSV — time, power, inlet and outlet temperature — into
a `DataFrame`. It handles two real-world data issues automatically:

- **Non-uniform time steps.** If any step drifts from the median by more than `dt_tol` (default
  0.1 s), the whole series is PCHIP-interpolated onto a uniform grid at the median step. This is
  common with data-acquisition systems that use adaptive time stepping.
- **A recirculation phase before heating.** Many rigs log a period of pump circulation, with the
  heater still off, before the test proper starts. `trim_recirculation = true` (the default) keeps
  only the single sample immediately before heating turns on — the undisturbed ground temperature
  `T0` — and rebases `:elapsed_time` so that sample is `t = 0`.

The mean fluid temperature `:T_mean` is added using the arithmetic mean; use
[`mean_fluid_temperature`](@ref) directly afterwards for a different averaging convention.

## Splitting heating and recovery

[`decompose_trt`](@ref) splits the loaded `DataFrame` into a [`TRTDataset`](@ref) on a heating-power
threshold (default 100 W). `heating` and `recovery` each carry a `:t_rel` column — time rebased so
`t_rel = 0` at the instant heating starts — which is what every `TRTDataset`-based fitting method
in [First-order approximation](@ref) and [Model inversion](@ref) actually consumes. `data` is the
full test (both phases), used by the model-inversion methods, which fit the whole signal at once
rather than one phase at a time.

`TRTDataset` assumes exactly one heating phase followed by one recovery phase — a **conventional**
TRT. An unconventional test (several on/off cycles) does not fit that shape; interpret it with the
vector-based `(t, T, q, ...)` method signatures directly on `load_trt_data`'s output.

## Temperature derivatives

The derivative-based FOA methods ([`fit_ils_foa_dT`](@ref), [`fit_ils_foa_dT_recovery`](@ref)) need
`dT/dt`. Two estimators are provided, both returning the same physical quantity and interchangeable
wherever a derivative is required:

- [`centered_finite_difference`](@ref) — the plain three-point finite difference.
- [`bourdet_derivative`](@ref) — brackets each point at a fixed logarithmic spacing `δ` instead of
  its immediate neighbors (linearly interpolating the temperature there). This trades local
  resolution for robustness to measurement noise, which matters on real TRT signals far more than on
  synthetic ones, and is why it is the default derivative estimator for the `TRTDataset`-based FOA
  overloads.

## Validity windows

Both FOA phases have a validity window past which the underlying infinite-line-source approximation
holds:

- [`critical_time`](@ref) — the Fourier-number threshold `tc = 5 rb²/(k/Cs)` for the
  temperature-based methods.
- [`residence_time`](@ref) — the time for fluid to travel down and back up the borehole, the unit
  the derivative methods express their window in (`[64 tr, 512 tr]` by default).

## Variable-power signals

[`step_signal`](@ref) segments a noisy, variable-power signal into a chosen number of constant
steps by k-means clustering — useful when a test doesn't hold power perfectly flat and the fitting
methods expect (or are more accurate with) a piecewise-constant load.

## Functions on this page

```@docs
TRTDataset
load_trt_data
decompose_trt
mean_fluid_temperature
centered_finite_difference
bourdet_derivative
critical_time
residence_time
step_signal
```
