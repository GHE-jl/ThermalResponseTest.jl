# API reference

This page is an alphabetical index of every exported symbol. Each entry links to its full
docstring, which lives next to the relevant theory or usage page.

```@index
Modules = [ThermalResponseTest]
```

## By topic

### Data loading and pre-processing

- [`TRTDataset`](@ref) — heating/recovery/full-test container
- [`load_trt_data`](@ref) — read a TRT CSV into a `DataFrame`
- [`decompose_trt`](@ref) — split into heating and recovery phases

### Utilities

- [`mean_fluid_temperature`](@ref) — p-linear mean of inlet/outlet temperature
- [`centered_finite_difference`](@ref), [`bourdet_derivative`](@ref) — `dT/dt` estimators
- [`critical_time`](@ref), [`residence_time`](@ref) — FOA validity-window helpers
- [`step_signal`](@ref) — k-means segmentation of a variable-power signal

### First-order approximation (Pasquier 2018)

- [`fit_ils_foa_T`](@ref) — UFOA-T-H, temperature, heating phase
- [`fit_ils_foa_T_recovery`](@ref) — UFOA-T-R, temperature, recovery phase
- [`fit_ils_foa_dT`](@ref) — CFOA-Ṫ-H, temperature derivative, heating phase
- [`fit_ils_foa_dT_recovery`](@ref) — CFOA-Ṫ-R, temperature derivative, recovery phase

### Model inversion wrappers

- [`fit_ground_response`](@ref) — generic core, any `AbstractGroundModel`
- [`fit_ils`](@ref) — infinite line source
- [`fit_ics`](@ref) — infinite cylindrical source
- [`fit_fls`](@ref) — finite line source
- [`fit_mils`](@ref) — moving infinite line source (groundwater advection)
- [`fit_mfls`](@ref) — moving finite line source (groundwater advection)
