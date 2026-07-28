# First-order approximation

All four methods on this page are derived in Pasquier (2018) from the infinite line source (ILS)
g-function, linearized in its long-time limit. They differ in which phase of the test (heating or
recovery) and which quantity (temperature or its time derivative) they regress.

## Temperature — heating phase (UFOA-T-H)

Past the ILS critical time ``t_c``, the mean fluid temperature is linear in ``\log t``:

```math
T_f(t) = \frac{q}{4\pi k}\Big[\ln t + \ln\!\big(\tfrac{4k}{C_s r_b^2}\big) - \gamma\Big] + T_0 + R_b^\ast q
```

with ``\gamma`` the Euler–Mascheroni constant. [`fit_ils_foa_T`](@ref) regresses ``T`` against
``\ln t`` on the samples past ``t_c`` to get a slope and intercept, recovers `k` from the slope
(``k = q/(4\pi\,\text{slope})``), and iterates — since ``t_c`` itself depends on `k` — until both
converge. `Rbₑ` then follows directly from the intercept, the fitted slope, and the known `T0`
(Eq. 3 of Pasquier 2018). `T0` must be supplied; it is not observable from the heating phase alone.

If fewer than 3 points ever fall past the (moving) ``t_c``, the iteration falls back to the last 30%
of the data and emits a warning rather than failing outright.

## Temperature — recovery phase (UFOA-T-R)

Once the heater switches off (at duration ``\bar t`` after the start of heating) while circulation
continues, superposition with a second, negative pulse gives a temperature linear in
``\ln\!\big(t/(t-\bar t)\big)`` instead of ``\ln t`` (Eq. 13 of Pasquier 2018):

```math
T_f(t) = \frac{q\bar t}{4\pi k}\,\ln\!\left(\frac{t}{t - \bar t}\right) + T_0, \qquad t > \bar t
```

[`fit_ils_foa_T_recovery`](@ref) mirrors the heating-phase iteration, offsetting the critical-time
check by ``\bar t``. Its intercept estimates ``T_0`` directly — a convenient cross-check against the
value used (or assumed) for the heating-phase fit — so, unlike [`fit_ils_foa_T`](@ref), it takes no
`T0` argument. `q` here is the **heating-phase mean** power per length, a scalar, not a vector: the
recovery-phase temperature only depends on the total heat injected, not its detailed time history
during heating.

## Temperature derivative — heating phase (CFOA-Ṫ-H)

Differentiating the heating-phase expression with respect to time removes the additive `T0` and
`Rbₑ q` terms entirely, leaving a pure power law (Eq. 6 of Pasquier 2018):

```math
\dot T(t) = \frac{q}{4\pi k\, t}
```

[`fit_ils_foa_dT`](@ref) works in log space, ``\ln t + \ln\dot T = \ln\!\big(q/(4\pi k)\big)``,
averaging the left-hand side over the chosen window to solve for `k` directly — no iteration needed,
since the critical-time dependence on `k` has been eliminated along with the terms it multiplied.

The window is expressed in fluid residence times `tr` rather than the physical critical time,
capturing the different physics that bounds a derivative fit: too early and the signal is
contaminated by borehole thermal-mass effects the ILS doesn't model; too late and measurement noise
dominates a vanishingly small derivative. Pasquier's original window is `[4 tr, 16 tr]`; this
package defaults to the wider `[64 tr, 512 tr]` because on real (as opposed to synthetic) TRT
signals, the derivative's log-log trend has often not yet settled onto its asymptotic unit slope by
`16 tr` — the original window is, in Pasquier's own words, "purely arbitrary," chosen for being
computable from readily available parameters rather than derived from first principles. Pass an
explicit `tr` or `indices` to override either window.

## Temperature derivative — recovery phase (CFOA-Ṫ-R)

The recovery-phase derivative follows the same idea, differentiating the recovery temperature
expression (Eqs. 14–18 of Pasquier 2018):

```math
\dot T(t) = -\frac{q\bar t}{4\pi k\, t(t-\bar t)}, \qquad t > \bar t
```

[`fit_ils_foa_dT_recovery`](@ref) again solves for `k` in one shot from the log-space average, with
the same `[64 tr, 512 tr]` default window (now measured from the *start of recovery*, i.e.
``t > \bar t + 64\,t_r``).

## Choosing a derivative estimator

Both derivative methods need `dT/dt` as an input, computed ahead of time by either
[`centered_finite_difference`](@ref) or [`bourdet_derivative`](@ref) — see
[Data & utilities](@ref) for the two estimators. The `TRTDataset`-based overloads default to the
Bourdet estimator, whose noise robustness is what makes the wider `[64 tr, 512 tr]` window usable at
all: the raw finite difference's sign-flipping noise on real signals would otherwise break the
log-space regression well before `512 tr`.

## Functions on this page

```@docs
fit_ils_foa_T
fit_ils_foa_T_recovery
fit_ils_foa_dT
fit_ils_foa_dT_recovery
```
