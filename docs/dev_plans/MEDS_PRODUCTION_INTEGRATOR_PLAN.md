# MEDS production numerical scheme — the plan

**Status:** design, rewritten 2026-07-31 (second revision, same day). The first revision replaced the
retired default-choice plan with a production-efficiency plan; this one adds the **stability analysis
of §1**, which changes what the workstreams should be, and reorganises into two parts:

- **Part I — testing and comparing integrators** (§3–§4): the measurement program. Nothing in Part II
  is built before the relevant Part I number exists.
- **Part II — numerical improvement and multi-core** (§5–§8): the build program.

**Question this answers:** production runs were bounded by a `dt_fast` of 150 s set by *stability*,
not accuracy — 6× more fast steps than the model used to take, each paying a full frozen pre-pass.
What raises that bound, what makes each step cheaper, and what parallelises — in what order, and gated
on which measurement?

**Answered, 2026-07-31 → 2026-08-02.** The stability bound is **gone**: one frozen coefficient (the
CAS↔atmosphere conductance) carried 99% of it, and re-solving it per stage — 2% of the pre-pass —
removes it while *reducing* integrator work. That landed as **PR #90** (merged), together with
`dt_fast = 900 s` and the non-stomatal water-stress limb defaulted off.

**The stability question is now closed, and closed on the right evidence.** The claim used to rest on
`Phi'` for a *single* state (canopy-air enthalpy), with the standing caveat that this is one diagonal
entry and not a spectral radius. §1i measures the **full six-state outer-map Jacobian**: `rho < 1` at
every cadence tested and *decreasing* with `dt_fast` (0.990 at 150 s → 0.943 at 900 s). Nothing
further is needed on the `gah` seam — see the N2b/N2c retirement in §5.

**`dt_fast` is now bounded by accuracy in exactly one variable.** At 900 s every state and every flux
is converged — GPP 1.0005, ET 0.9997, NEE 1.0002; `T_cas` 0.053 K, `T_gnd` 0.018 K, `CO2_cas`
0.29 ppm — **except `psi_leaf`, at 0.76 MPa**, an outlier by ~1.5 orders of magnitude. And §1i finds
that error is **inherited from the canopy air and amplified ~4×**, not generated autonomously by the
hydraulics. That reframes N2b (§5) and is the top of the queue.

**Two proposals died on measurement this round:** the ground-flux analogue of N2a (§1i.1 — the 0.70 K
soil-T figure does not reproduce, and the mechanism contradicts the scaling), and the two cheaper
`gah` variants (§5 N2b/N2c-on-`gah`). Multi-core (§7) is untouched and is now the largest open item.

**Companions.** `docs/science/numerical_scheme.md` is the user-facing description of the two schemes
and the freeze. `MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` §9–14 is the measurement record for the
oscillation. `MEDS_NUMERICS_SCOPING.md` §11 is the bare-array kernel convention. The archived parity
document is history and is not evidence — see §2b.

---

## 0. What is void from the previous version

The original subject of this file — which integrator family becomes the default — is **closed**.
`split` is retired, `ark` (ESDIRK2) is the default, `rk45` is the accuracy baseline, and its stiff
rescue redoes the step on `ark` (a within-family *stepper* swap, not a model swap). PR #88 closed the
last physics rows: `numerical_scheme.md` §4 now reads "there are no *physics* rows left."

| previous item | disposition |
|---|---|
| F3 — the `split` ↔ ARK/RK45 ~0.45 K gap; P1-a "attribute it" | **VOID.** The family was deleted, so there is nothing to attribute it to. |
| F4 — `ark @ 900 s` beats `split @ 225 s` on the efficiency frontier | **VOID** (one arm retired) and **stale** (both cells inside the oscillating regime). |
| P0-a…e — the five correctness blockers | **All landed.** `[fast].snow_on` and the `with_mass`/`with_theta` norm switches deleted. |
| P1-a′ — the three Class-1 physics rows | **CLOSED** (PR #88). |
| P2 — retire the cross-family fallback | **DONE** (`meds_fast_step.f90:97-113`). |
| P3 — a validator whose main rule was "reject `split` below 900 s" | **Rule void**; the cost-ceiling fragment survives into §6. |
| The layer-3 `dt_fast` → expected-error table | **Deleted** — every entry was measured inside the oscillating regime (§2b). |
| Workstream A variants A1/A2 (aerodynamics / surface coupling explicit) | **Refuted** before being built; retained only as §5 N6. |
| Workstream B3 as first drafted ("quasi-steady CAS, biggest stability win") | **Inverted by §1.** Naive `wcap → 0` is the *maximally unstable* configuration. Rewritten as N4. |

The one structural rule that survives intact — *adapt the step size within a run, choose the scheme at
configuration time* — is now honoured by construction, since the RK45 rescue is a stepper swap on one
semi-discretisation. §5 N5 proposes the first principled exception, at layer 3 rather than layer 2.

---

## 1. The governing diagnosis

### 1a. The cost model

```
  cost  =  n_patch  x  (dt_slow / dt_fast)  x  [ C_freeze  +  n_substep x C_stage x 2 ]
                        ^^^^^^^^^^^^^^^^^        ^^^^^^^^     ^^^^^^^^^^^^^^^^^^^^^^^
                        freeze cadence           pre-pass      the ESDIRK march
                        (stability-bounded)      (once per     (adaptive, ~2-3 substeps
                                                  dt_fast)      x 2 implicit stages)
```

Raising `dt_fast` is the only lever that cuts *both* bracketed terms, which is why §5 leads with it.
`C_freeze` is very probably dominant and **has never been measured** — T1.

### 1b. The canopy-air storage term is not what is failing

The canopy air is already implicit. `column_be_stage` commits it in flux form,

```
H_{n+1} = ( wcap*H_n + dt*(src_enth + gah*H_atm) ) / ( wcap + dt*gah )
```

inside an L-stable ESDIRK2, with the leaf↔CAS pair closed by a 2×2 Newton every stage
(`newton_surface_solve`). The fast mode is handled by an unconditionally stable method; small `wcap`
is not breaking the time discretisation.

Write the memory weight the old state carries:

```
a  =  wcap / (wcap + dt*gah)  =  1 / (1 + dt/tau),      tau = wcap/gah
```

**Measured** on the reference fixture at midday (§1g): `wcap = 24 kg/m2`, `gah = 1.68e-2 kg/m2/s`,
so `tau ≈ 1430 s` — an order of magnitude longer than a first estimate of 80–120 s assumed, because
`u_ref = 2 m/s` over an 18 m canopy gives `ustar = 0.167`, not the 0.3–0.6 of a windy day.

| `dt_fast` | `dt/tau` | `a` — fraction of `H_n` surviving |
|---|---|---|
| 900 s | 0.63 | **0.61** |
| 225 s | 0.16 | 0.86 |
| 150 s | 0.10 | 0.91 |
| 100 s | 0.07 | 0.93 |

So at production settings the canopy air is **not** nearly-diagnostic — it retains 60–90% of its state
across a step, and the storage term is doing real damping. The instability therefore cannot be blamed
on the storage treatment at all; it is entirely in the lag (§1g).

### 1c. The instability is an outer map on the frozen coefficients

`gah`, `g_tr_f`, `h_coeff_f` and `t_emit` are frozen at state `n`, so across one `dt_fast` the
canopy air obeys a map `H_{n+1} = Phi(H_n)`. Linearised, it splits into storage memory plus lag
feedback:

```
Phi'  =  a  +  (1 - a) * L        L = the quasi-steady map's coefficient-lag gain
```

`a` is positive and bounded by 1 — it is the **contraction**. `L` is the destabiliser.

**An earlier draft of this section fitted `L ≈ −2.3` to the amplitude scan. That estimate is
REFUTED — `Phi'` has now been measured directly (§1g) and the true lag gain is an order of magnitude
larger.** The structural decomposition survives; the number did not. Amplitude is a property of the
*nonlinear saturation*, not of the multiplier, and fitting one to the other was the error.

### 1d. Consequences that invert two earlier proposals

**Naive `wcap → 0` (a diagnostic canopy air with frozen coefficients) is the worst case, not the fix.**
It sets `a = 0`, deleting the damping the measured `a = 0.61–0.91` is currently supplying, and pins
the multiplier at `L ≈ −40 to −60`. This inverts the first draft's B3.

**Inner sub-stepping is very slightly harmful, which explains a previously unexplained observation.**
`n` sub-steps of `dt/n` give `a_eff = (1 + dt/(n·tau))^-n`, which *decreases* with `n`, driving the
canopy air closer to the frozen-coefficient equilibrium. `MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` §9
recorded that "no amount of inner adaptivity fixes it" as an empirical surprise; it is a prediction of
this model.

**"Diagnostic" has two meanings and they land on opposite sides.** The second one works:

- `wcap → 0` with coefficients still frozen at state `n` — an explicit map with no damping. Above.
- An **index-1 DAE**: solve `S(H, c(H)) + gah(H)·(H_atm − H) = 0` with the coefficients evaluated at
  the *same* `H` as the balance. There is no lag, so `L` does not exist. This is what CLM5 does — its
  canopy air *is* diagnostic, and it works because the Monin–Obukhov and stomatal updates sit inside
  the surface-energy iteration. ED2 gets the same effect by recomputing `canopy_turbulence8` at every
  RK stage.

**Therefore "solve the CAS algebraically" and "put `∂c/∂y` into the Newton" are the same fix.** The
2×2 Newton already treats `H` and `q` implicitly; what it holds fixed are the Jacobian entries coming
through the coefficients. Adding them converts the lagged map into a genuine implicit solve. That is
the unifying statement of §5.

### 1e. Why the adaptive controller cannot find this by itself

A fair objection: MEDS already runs an embedded-error adaptive march — should it not discover that it
needs smaller steps? **No, and the reason is structural.**

The controller integrates `y' = f(y, c_frozen)` over `[0, dt_fast]` and controls the **local
truncation error of that ODE**. The frozen problem is smooth and the method L-stable, so the embedded
pair agrees closely, `err ≈ 0`, and the controller grows the step. The error actually being made is
`c_frozen ≠ c(y(t))` — a **semi-discretisation error, not a truncation error**. Both members of the
embedded pair solve the *same wrong right-hand side*, so their difference is blind to it by
construction. And the oscillation is a property of the map *between* `dt_fast` steps, which no
within-step estimator ever samples.

This is the same blindness as the conservation ledgers, for the same reason: every detector currently
in the code measures a within-step property of the frozen problem. Three independent instruments —
the ledgers, the embedded estimator, and the Newton's own convergence test — all report health while
the canopy air swings 8 K step to step.

**The corollary is the design:** give the freeze its own error estimator. The next step's pre-pass
computes `c(y_{n+1})` anyway, so the coefficient drift `||c(y_{n+1}) − c(y_n)||` is available at zero
extra cost, and it is the natural controlled quantity for an adaptive *freeze cadence* (N5).

### 1f. What this predicted, and how it scored

| # | prediction | outcome |
|---|---|---|
| 1 | `Phi'` crosses −1 between `dt_fast` 150 s and 225 s | **WRONG** — it crosses between **75 s and 100 s** (§1g T2a). The default 150 s is unstable, `Phi' = −2.74`. |
| 2 | `Phi' = a + (1−a)L` with `L` roughly constant | **PARTLY** — the form holds and `Phi'` is monotone in `dt_fast`, but the implied `L` drifts −37 → −62 across the sweep, so the one-mode model is a caricature. |
| 3 | under-relaxation with `θ < 0.5` stabilises 900 s | **WRONG** — with the measured `L ≈ −62`, `θ < 0.07` would be needed at 900 s. N1 is demoted. |
| 4 | a young stand is *stricter* than 150 s | **REFUTED** — the relationship is non-monotone. The worst stand is **mid-height (10 m, `Phi' = −8.2` at 150 s)**; the 2 m regenerating stand is the *most* stable of the five (§1h). |
| 5 | the tissue store makes 150 s pessimistic | **UNTESTED** — moot once N2a lands, since the CAS is then stable at every cadence measured. |
| 6 | the frozen ground flux carries the 0.70 K soil-T error at 900 s; the ground is "the natural next target" | **REFUTED** (§1i.1) — 0.70 K does not reproduce across a 6× sweep of ground forcing (0.054–0.117 K), and the error fails to scale with the driving. `ggnet` swings 81% within a step and it still does not matter: the store it feeds is heavy and `ggnet` is ~7× smaller than `ggbare`. No ground refresh built. |
| 7 | `Phi'` on canopy-air enthalpy is a good proxy for the map's stability | **HELD, but incompletely** (§1i.2) — the full six-state spectral radius is `rho < 1` everywhere and tracks the diagonal closely, so the proxy was not misleading *for stability*. It was blind to something else: an **amplifying-but-stable** mode (`∂(leaf_water)/∂(cas_shv)` = −3.75, loop gain 0.09) that no spectral radius and no single-state `Phi'` can see. Stability and accuracy needed separate instruments. |
| 8 | `psi_leaf` is non-converged because its own update is an explicit step with frozen sapflow | **INVERTED** (§1i.4) — its own diagonal is ≈0.09, i.e. it barely remembers its own state. The error is **inherited from the canopy air and amplified ≈4×**. The fix is coefficient-side, not integrator-side. |

### 1h. MEASURED, T3 and T5

**T3 — the `ark` ↔ `rk45` gap is an artifact of the instability, not a scheme difference.**
5 cohorts, 24 h, scored against `ark` + refresh at `dt_fast = 12.5 s`, `rtol_all = 1e-9`,
`atol_scale = 1e-3`:

| | `dt_fast` | 25 s | 50 s | 75 s | 100 s | 150 s |
|---|---|---|---|---|---|---|
| **frozen** | `\|ark − rk45\|` CAS-T, default tol | 0.0001 | 0.0002 | 0.0002 | 0.014 | **0.491** |
| **frozen** | ...tightened tol | 0.0001 | 0.0002 | 0.0002 | 0.010 | **0.483** |
| **refreshed** | `\|ark − rk45\|` CAS-T, default tol | 0.0009 | 0.0018 | 0.0027 | 0.0036 | 0.0074 |
| **refreshed** | ...tightened tol | 0.0005 | 0.0009 | 0.0012 | 0.0013 | **0.0013** |

With the refresh and a tight tolerance the two schemes — one implicit, one explicit — agree to
**1.3 mK at the production cadence**. Frozen, the same comparison reads 0.49 K, **370× larger**, and
tightening the tolerance does not shrink it. That 0.49 K is not a scheme difference: it is two
steppers sampling different phases of the same period-2 oscillation. This is the same trap that
produced the retired `split` gap, and it is worth stating as a rule: **a cross-scheme difference
measured in an unstable regime measures the instability, not the schemes.**

**The freeze dominates and the integrator tolerance is irrelevant — confirmed at a *stable* cadence.**
Refreshed `ark` at 150 s: 0.0174 K at default tolerance, 0.0180 K tightened. The time-stepping error
is ~0.6 mK against a 17 mK total, so ~97% of the error is the Category-0 freeze. (Frozen, the same
test is meaningless — tightening moves `ark` from 0.391 to 0.421 K, i.e. *worse*, because there is no
convergence to measure.)

**Convergence.** Refreshed `ark` errors 0.0006 / 0.0015 / 0.0024 / 0.0052 / 0.0174 K over
25 → 150 s: roughly **second order**. Frozen: 0.0011 / 0.0021 / 0.0032 / **0.153** / **0.391** — a
discontinuity between 75 and 100 s, which is the stability boundary, not a convergence curve. The
frozen scheme has no order across it.

**Carbon and water at 150 s** (against the same reference): GPP −3.8% refreshed vs −4.1% frozen; ET
**−2.6% refreshed vs −9.0% frozen**. ET improves 3.5×. GPP barely moves, because it is limited by the
*leaf gas-exchange* freeze, which N2a does not touch — that, and the ground fluxes, are what bound
`dt_fast` once the canopy air is stable.

`rk45_rescue = 0` in every cell, so every RK45 number above is genuinely RK45.

**T5 — the default is unsafe on most stands, and prediction #4 was backwards.**
`Phi'` at the probe state, `can_depth` following the stand (`max(min_canopy_depth, h + freeboard)`):

| stand | `wcap` | `Phi'` @ 50 s | @ 100 s | @ 150 s | @ 900 s | p2p @ 150 s |
|---|---|---|---|---|---|---|
| 2 m regen, LAI 1 | 8.4 | +0.16 | −0.44 | **−0.90** | −2.48 | 0.008 K |
| 5 m, LAI 2 | 12.0 | −0.93 | −2.51 | **−3.84** | −10.3 | 3.36 K |
| **10 m, LAI 3** | 18.0 | **−2.50** | −5.53 | **−8.22** | −25.9 | 2.07 K |
| 18 m, LAI 3 | 27.6 | −1.15 | −3.08 | −4.84 | −19.6 | 1.70 K |
| 30 m, LAI 4 | 42.0 | −1.10 | −3.05 | −4.88 | −23.3 | 0.79 K |

**Four of five stands are unstable at the default `dt_fast = 150 s`**, and the worst is *worse* than
the 18 m fixture the default was calibrated on (−8.2 vs −4.8). The frozen scheme would need
`dt_fast ≈ 50 s` — 3× more steps than today — to be stable on a 10 m stand. So N2a is not only an
optimisation: **without it, the current production default is numerically unsafe on ordinary stands.**

The prediction that a *shorter* stand is stricter is refuted. Capacity does fall with height, but so
does the flux driving it: over this height–LAI sequence the ratio LAI/`wcap` peaks in the middle
(0.17 at 5–10 m against 0.12 at 2 m and 0.10 at 30 m), and `Phi'` tracks that, not `wcap` alone. The
mechanism in §1c survives; the monotone corollary does not. **Caveat: LAI was varied with height, so
height and leaf area are confounded by construction — a 2 m stand carrying LAI 3 was not tested.**

**With the refresh, every stand is stable or marginal at every cadence measured**: `Phi'` at 150 s
runs −0.07 / −0.47 / −0.75 / −0.89 / −1.16 across the five stands, and peak-to-peak stays ≤ 0.04 K at
150 s and ≤ 0.44 K at 900 s. The 30 m stand is marginal (`Phi' = −1.16` at 150 s, −1.02 at 100 s) but
its amplitude is 0.017 K, so it is a weak mode rather than a live oscillation — worth re-checking on a
real tall-canopy configuration.

Also visible, and a caution for anyone reading amplitudes: the 2 m stand has the *mildest* multiplier
at 900 s (−2.48) and the *largest* excursion (20.4 K p2p). Multiplier and amplitude are decoupled by
the saturating nonlinearity. Rank stability by `Phi'`, not by peak-to-peak.

**Caveat on the probe.** The T5 base state was spun up at 50 s with the refresh off, which is itself
inside the unstable regime for the 5–30 m stands, so each probe is taken somewhere on a limit cycle
rather than at a fixed point. `Phi'` is a local derivative and the comparison across stands is
like-for-like, but the absolute values carry that.

### 1g. MEASURED, 2026-07-31 — the multiplier, its cause, and the cost split

Harness: the `test_column_dynamics` high-LAI sunlit column (LAI 3, 18 m canopy, single cohort),
spun up once to 12:00 solar and probed from that **one common base state**; central differences on
`bio%cas%can_enthalpy` applied *before* `build_column_frozen`, so what is measured is the lagged
outer map. `build_column_frozen` and `adaptive_ark_march` are public, so the coefficient
decomposition needed no source change.

**(i) The map is smooth, and there is no stable `dt_fast` near the default.**

| `dt_fast` | 25 s | 50 s | 75 s | 100 s | 150 s | 225 s | 450 s | 900 s |
|---|---|---|---|---|---|---|---|---|
| `Phi'` | +0.36 | −0.26 | −0.88 | **−1.49** | **−2.74** | **−4.64** | **−10.6** | **−23.2** |

The stability boundary is at **~75–80 s**, not 150 s. `Phi'` converges cleanly as the perturbation
shrinks (−23.40, −23.37, −23.21 at δ = 10, 50, 100 J/kg, forward and backward differences agreeing to
1%), so **the map is differentiable** — a Jacobian-based fix is viable, and the switching-nonlinearity
worry is dismissed. Nonlinearity appears only above δ ≈ 0.5 K, where `|Phi'|` *falls* — a saturating
nonlinearity, which is what bounds the limit-cycle amplitude and is why amplitude is a poor proxy for
the multiplier.

**(ii) One coefficient carries essentially all of it: the CAS↔atmosphere conductance.** Holding each
frozen coefficient at its *unperturbed* value (i.e. "this coefficient did not see the perturbation"):

| pinned | `Phi'` @ 900 s | `Phi'` @ 150 s |
|---|---|---|
| nothing | −23.24 | −2.74 |
| **`gah`/`gaw`/`gac`** | **+0.80** | **+0.88** |
| `g_tr_f` (stomatal) | −23.33 | −2.77 |
| `h_coeff_f` (leaf sensible) | −23.29 | −2.76 |
| `abs_lw` (`t_emit`) | −23.24 | −2.74 |
| `f_wet_c` | −23.24 | −2.74 |
| all of the above | +0.66 | +0.83 |

Removing the aerodynamic conductance's lag alone takes the map from `Phi' = −23.2` to **+0.80 —
stable at 900 s**. Everything else contributes ≤ 0.7%. The mechanism is the Monin–Obukhov stability
feedback: a warmer canopy air is a more unstable surface layer, which raises `ustar` and the transfer
factor, which vents the canopy air harder — a negative feedback of very high gain (measured
`d ln gah/dT ≈ 2.2 K⁻¹`: +0.5 K gives `gah` ×3.1, +2 K gives ×12.2), lagged by exactly one `dt_fast`.

**This reconciles with, rather than contradicts, the earlier "aerodynamics only amplifies" result.**
That experiment pinned `ustar`. But `gah = rho·ustar·temp1`, and `temp1` — the MO scalar-transfer
factor — is *also* stability-dependent and was left free. At +0.5 K the measurement gives `ustar` ×1.75
against `gah` ×3.06, so pinning `ustar` alone removes only about half the log-sensitivity. The prior
test did not pin what it intended to pin.

**(iii) It is a light-to-moderate-wind phenomenon, which is the common case.** `Phi'` (full / with
`gah` pinned):

| `u_ref` | `ustar` | `Phi'` @ 150 s | `Phi'` @ 900 s |
|---|---|---|---|
| 1 m/s | 0.100 (at the floor) | +0.31 / +0.89 | −3.63 / +0.94 |
| **2 m/s** | 0.167 | **−2.74** / +0.88 | **−23.2** / +0.80 |
| 5 m/s | 0.630 | −0.48 / +0.48 | −2.42 / −0.02 |
| 10 m/s | 1.640 | −0.066 / +0.04 | −0.30 / −0.16 |

Strong winds are stable (high ventilation, weak stability sensitivity). Very light winds are stable
too, but only because `ustar` is pinned at its 0.1 floor, which incidentally zeroes the feedback — an
accidental stabiliser, not a designed one. The danger band is **light-to-moderate wind**, which is
the ordinary condition over a forest canopy.

**(iv) T1 — the pre-pass dominates, and the aerodynamics is a rounding error inside it.**
Copy-overhead subtracted, `dt_fast = 150 s`, 10 soil layers, 20 000 repeats:

| cohorts | `C_freeze` | whole sub-step | freeze share |
|---|---|---|---|
| 1 | 14.7 µs | 37.4 µs | 39% |
| 3 | 30.0 µs | 53.4 µs | **56%** |
| 10 | 86.9 µs | 137.2 µs | **63%** |
| 30 | 325.7 µs | 366.7 µs | **89%** |

The freeze scales with cohort count (leaf gas exchange, ~11 µs/cohort); the march does not. At any
realistic patch occupancy **`C_freeze` is the majority of the cost**, so cheapening the march is
low-value and §6 is confirmed as minor. And within the pre-pass, `canopy_aerodynamics` is **2.0%**
(0.28 µs of 13.9 µs) — so refreshing it at four stage points costs **~8% of one pre-pass**, roughly
3% of a sub-step.

**(v) The synthesis.** The single coefficient that causes the instability is also one of the cheapest
things in the pre-pass. Refreshing `gah`/`gaw`/`gac` per stage — and nothing else — is predicted to
move the stability limit from ~75 s to beyond 900 s for ~3% more work per sub-step, i.e. **up to a 6×
reduction in fast-loop cost.** That is now N2 and it is the plan's headline.

**Caveats.** `Phi'` is `∂H_{n+1}/∂H_n`, one diagonal entry of the map's Jacobian, not its spectral
radius; it is a strong indicator, not a proof of stability. One fixture, one time of day, one cohort,
one soil column, constant-forcing MVP (no live met, no RT join). The pinned-`gah` values (+0.80 at
900 s) are close enough to 1 that a different stand could still sit outside. Reproduce with
`scratchpad/meas_fastloop{2,3,4}.f90`. **The spectral-radius caveat is discharged in §1i.2.**

### 1i. MEASURED 2026-08-02 — the ground analogue, and the full-Jacobian checklist

Harnesses `scratchpad/meas_ground.f90` and `scratchpad/meas_checklist.f90` (+ `build.sh`), linked
against the built libraries with zero source instrumentation. Fixture: 3 cohorts, LAI 5.1, 18 m stand,
`u_ref` 2 m/s, 10 soil layers, spun up **once** to a common base state.

#### 1i.1 The ground analogue of N2a is REFUTED. Do not build it.

The N2a block below asserted soil-T 0.70 K at 900 s "because the ground fluxes remain frozen", and
recommended `dt_fast` 225–450 s on that basis. Measured across three regimes spanning a 6× range of
ground forcing:

| regime | `T_gnd` diurnal amplitude | soil-T RMSE @ 900 s | err/amp |
|---|---|---|---|
| damp soil, shaded ground | 4.5 K | **0.054 K** | 1.2% |
| dry soil, shaded ground | 3.6 K | 0.117 K | 3.2% |
| dry soil, **open** ground | **29.0 K** | 0.115 K | **0.4%** |

**0.70 K does not reproduce in any regime** — it is 6–13× smaller — and the attributed mechanism
contradicts the scaling: the ground forcing varies 6× across regimes while the error barely moves. The
original fixture could not be recovered, so the defensible statement is *"does not reproduce across a
6× sweep"*, not *"the old number was wrong"*.

Two supporting measurements, which cut in opposite directions and must be read together:

- **`ggnet` genuinely swings up to 81% within one step** (vs `gah` 47–131%). "The coefficient moves"
  is TRUE. But its swing is nearly **`dt_fast`-INDEPENDENT** (0.769 at 150 s vs 0.812 at 900 s),
  whereas `gah`'s grows properly with `dt_fast` (0.47 → 1.31). A `dt`-independent swing is the diurnal
  cycle plus the `stab` switch inside `ggveg` — **not a lag**.
- **The ground does care about `ggnet`.** Permanently scaling `cs_dense` by 5× moves `T_gnd` by 0.83 K
  (damp) to 4.10 K (dry/open). So this is not insensitivity: the *lag* error is transient and averages
  out, while a permanent bias does not.

Mechanism: `ggnet ≈ 3.4e-4 m/s`, about **7× smaller than `ggbare`** — `ggveg = cs_dense·ustar`
dominates the CLM blend, so the ground is nearly decoupled from the canopy air aerodynamically. That
is why an 80% conductance error buys so little temperature error.

> **Generalisable lesson.** *A coefficient moving a lot is not evidence that freezing it costs
> anything.* `gah` mattered because the canopy air has tiny capacity **and** the loop closes on
> itself. `ggnet` moves just as much and does not matter, because the store it feeds is heavy and the
> coupling is weak. Measure the **error**, never the coefficient swing.

**Consequence:** the contradiction between this document's two `dt_fast` recommendations resolves in
favour of **900 s**. The sentence "recommend `dt_fast` 225–450 s, not 900 s, until something
equivalent is done for the ground surface" is **void**.

#### 1i.2 Full outer-map Jacobian — `rho < 1` everywhere

Six states (`cas_enthalpy`, `cas_shv`, `cas_co2`, `soil_energy(1)`, `theta(1)`, `leaf_water_mass`),
central differences, non-dimensionalised by a per-state characteristic scale (≈0.2 K, 2e-4 kg/kg,
2 ppm, ≈0.2 K, 2e-3 m³/m³, 1% of the leaf pool) so the entries are comparable. Spectral radius by
power iteration, plus `sqrt(rho(J²))` to catch a complex pair — which is precisely the period-2
signature the CAS oscillation had.

| `dt_fast` | 150 s | 300 s | 450 s | 900 s |
|---|---|---|---|---|
| max abs diagonal | 0.990 | 0.980 | 0.971 | 0.943 |
| **spectral radius** | **0.990** | **0.980** | **0.971** | **0.943** |

**`rho < 1` at every cadence, and it decreases with `dt_fast`.** The map contracts; the largest
eigenvalue is the slow soil/`theta` mode, not the canopy air. N2a survives the full-matrix test, so
the §1g caveat ("one diagonal entry, not the spectral radius") is discharged.

The CAS-enthalpy diagonal is **non-monotone** — −0.30, −0.13, −0.10, **−0.51** at 150/300/450/900 s.
Worth watching; far from 1. (§1g reported −0.14 at 900 s on a 5-cohort fixture; −0.51 here on 3.)

#### 1i.3 The checklist: everything converges at 900 s except `psi_leaf`

24 h, scored against a 25 s reference, `wstress_nonstomatal` OFF (the shipped default):

| `dt_fast` | `T_cas` [K] | `q_cas` [g/kg] | `CO2_cas` [ppm] | `T_gnd` [K] | `theta1` | **`psi_leaf` [MPa]** |
|---|---|---|---|---|---|---|
| 150 s | 0.005 | 0.008 | 0.034 | 0.007 | 1e-5 | **0.100** |
| 900 s | 0.053 | 0.065 | 0.288 | 0.018 | 4e-5 | **0.764** |

Day-total fluxes at 900 s: **GPP 1.0005, ET 0.9997, NEE 1.0002** — converged to 0.05%. (ET is taken
from `budg%cas_water%outflux`, the same number the conservation ledger charges — "one flux, both
sides". Only the *ratio* is validated; the absolute was not checked against a physical ET.)

`psi_leaf` is the sole outlier, by roughly 1.5 orders of magnitude against every other state.

#### 1i.4 `psi_leaf`'s error is INHERITED and AMPLIFIED, not autonomous

This changes what N2b should be. The two largest entries in the entire Jacobian are both **into**
`leaf_water` **from** the canopy air, and both grow superlinearly with `dt_fast`:

| entry | 150 s | 300 s | 450 s | 900 s |
|---|---|---|---|---|
| `∂(leaf_water)/∂(cas_shv)` | −0.39 | −1.01 | −1.67 | **−3.75** |
| `∂(leaf_water)/∂(cas_enthalpy)` | +0.12 | +0.29 | +0.47 | **+1.03** |
| `∂(leaf_water)/∂(leaf_water)` *(own)* | 0.085 | 0.092 | 0.094 | 0.099 |

Leaf water is a **one-way amplifier of canopy-air error** — gain ≈4× at 900 s and rising — with almost
no memory of its own state (diagonal ≈0.09: it nearly fully relaxes within one step). The return path
is weak (`∂(cas_shv)/∂(leaf_water)` = +0.023), so the loop gain is ≈0.09 and **the system stays stable
while amplifying heavily**. That is exactly why the spectral radius never saw this, and why a
stability-only assessment would have declared the integrator healthy and stopped.

So `psi_leaf` is a **fast, near-algebraic variable being integrated as a prognostic one under frozen
coefficients** — the same category as the canopy air before N2a, one seam over. See N2b in §5.

**Caveats.** One fixture, one base state (midday), one wind speed; the Jacobian is a one-step
linearisation about that state. The leaf-water diagonal being nearly flat in `dt_fast` (0.085 → 0.099)
is **not** a simple exponential relaxation and is currently unexplained — it is the loose thread.

---

## 2. Evidence

### 2a. What is measured and current

> **SUPERSEDED BY N2a — this subsection is the PRE-REFRESH record.** Every number below was taken
> with `gah`/`gaw`/`gac` frozen for the whole `dt_fast`, which PR #90 removed. It is kept because it
> is the evidence that *identified* the defect and because §1g's reconciliation of the "aerodynamics
> only amplifies" result depends on it — not because it describes the shipped scheme. Current status:
> `dt_fast` is **no longer a stability boundary** (§1i.2, `rho < 1` everywhere), the default is
> **900 s**, and `meds_config`'s warning is about carbon, not stability. Read §1g/§1h/§1i for what is
> true now.

`dt_fast` is a stability boundary. Canopy-air temperature, peak-to-peak over *consecutive* steps,
high-LAI sunlit fixture (`MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` §9–10):

| `dt_fast` | 900 s | 450 s | 300 s | 225 s | 150 s | 100 s |
|---|---|---|---|---|---|---|
| CAS-T p2p | ~8 K | ~8 K | ~5 K | ~4 K | ~1 K | smooth |

Every conservation ledger closed to ~1e-6 J at every row. Default is 150 s; `meds_config` warns above
300 s.

**Pinning the entire surface coupling** — aerodynamic conductances, `g_tr_f` and `h_coeff_f` together
— takes 9.5 K to 6.7 K at 900 s. Pinning is strictly more aggressive than any refresh (it removes the
lag *and* the response), so per-sub-step refresh of the surface coupling cannot recover `dt_fast`.
Also eliminated, each still oscillating at 900 s: plant hydraulics (10.2 K), soil thermal column
(9.1 K), soil water column (10.0 K).

**Untested suspects for the residual `L`:** the leaf gas-exchange pre-pass (`gsw` solved once per
`dt_fast` at state `n`; pinning `g_tr_f` removed its *variation* while the pre-pass still ran at a
lagged state), the longwave emission base `t_emit`, and `f_wet`.

**A methodological gap that T4 closes.** Every experiment so far *pinned* a coefficient, which asks
"is X **necessary**?" — answered *no* for everything tried. It never asked "would refreshing X **cure**
it?" Under §1c those are different questions: pinning sets that coefficient's contribution to `L` to
zero *and* removes its physical response, whereas the perturbation decomposition in T2/T4 measures its
contribution to `L` directly.

Both schemes close whole-column water and energy to machine precision, with snow, condensation and
the tissue store active.

### 2b. What is void as evidence

Anything measured at `dt_fast` ∈ {900, 1800} s in a sunlit high-LAI regime is a **phase-sample of an
oscillation**, not a convergence measurement. That voids: the cross-scheme accuracy tables in the
retired parity document; "ARK time-stepping error 0.039 K vs total 1.10 K at 1800 s" and the 8×
refinement ratios; the deleted `dt_fast` → error table; ARK's rejection-storm counts; and the 30-year
`runs/ithaca_ark30/` beds.

### 2c. Code facts this plan assumes

Verified on `main` at `0f95e67`:

- `ark` is a 2-solve **ESDIRK2**, `f_E ≡ 0`. There is no explicit tableau branch to populate — a true
  IMEX is new machinery, not a flag (`ark2_column_step`).
- Inside the tableau: CAS enthalpy/humidity/CO₂ + the soil thermal column. **Operator-split out:** soil
  water (committed once from the scratch `column_hydrology_flux`), plant internal water mass, canopy
  surface water, plant hydraulics (an exact matrix exponential — deliberately not a stage).
- The leaf↔CAS pair is a direct 2×2 Newton with a **numerical** Jacobian (`jac_surface`):
  `NEWT_MAX = 4`, `LS_MAX = 6`, `FEVAL_CAP = 24`. A numerical Jacobian is why enlarging the implicit
  block (N2) is mechanical rather than a derivation exercise.
- `TISSUE_STORE_SCALE = 1.0` — the tissue store is on, carried as a b-weighted time integral, adding
  ~1.5e4 J/m2/K against the canopy air's ~3.0e4.
- Every state combinator (`state_init`, `state_axpy`, `state_accum`, `state_extrap`, `state_sub`,
  `state_err_diff`) allocates four per-cohort arrays on entry; one sub-step calls several. §6.
- The patch loop (`fast_dynamics`, `meds_fast_dynamics.f90:322`) is serial, and its per-patch scratch
  was deliberately hoisted *out* of the loop (BB1 phase 1) — exactly the structure that blocks
  threading. §7.
- `t_sub` depends only on `isub`, yet `met_advance` / `met_instant` / `apply_met_to_ctx` run **inside**
  the patch loop, repeating site-uniform work `n_patch` times per sub-step. §7 C1.
- CMake gives OpenMP flags to NVHPC only, and only on `meds_core`. The fast loop is in `meds_aux`.
- **Dead split residue:** `picard_max_iter` / `picard_tol_temp` / `picard_tol_shv` / `picard_relax` /
  `picard_fixed_iter` are plumbed to no reader; `budg%picard_iters` / `picard_nonconv` are read and
  never written; `[fast].ark_relax` is vestigial; `[fast].ark_niter` is only tested as `np <= 1`, so
  it is a boolean in disguise (the real cap is `NEWT_MAX = 4`).

---

# PART I — TESTING AND COMPARING INTEGRATORS

The model in §1 is a hypothesis fitted to five points. Part I is what turns it into a measurement, and
what any candidate scheme is then judged against. **Nothing in Part II is built before its gating
number exists.**

## 3. The measurement program

**Progress, 2026-07-31: T1, T2 and T4 are DONE — results in §1g.** They were run with a standalone
harness (`scratchpad/meas_fastloop{2,3,4}.f90`) against the built libraries, exploiting the fact that
`build_column_frozen` and `adaptive_ark_march` are public, so no source instrumentation was needed.
T3 and T5 remain outstanding. Two protocol errors were made and corrected on the way, both worth
recording because they are generic: a first pass spun up **separately per `dt_fast`**, so each row
probed a different state (`tau` varied 40× across rows and `Phi'` was non-monotone); and it timed
`adaptive_ark_march` in isolation with a **cold start**, so the isolated march did more work than the
warm-started in-context call and freeze + march exceeded the total.

### T1 — per-phase cost profile of one `dt_fast`  *(gates all of Part II)* — **DONE, §1g(iv)**

**Answer: the pre-pass is 56–89% of a sub-step at realistic cohort counts, and `canopy_aerodynamics`
is 2% of the pre-pass.** So cheapening the march is low-value, and the one coefficient that must stop
being frozen is also one of the cheapest to refresh.

Still worth building as permanent instrumentation, behind `[output].numerics`: wall-clock and call
counts for `column_prepass` (split out:
leaf gas exchange, aerodynamics, respiration), `canopy_radiation`, `solve_plant_water_batch`,
`column_hydrology_flux`, the ESDIRK march as a whole, and within it `surface_derivs` evaluations and
`soil_energy_step_implicit` solves. Report per patch and per site, summer and winter month, at
`dt_fast = 150 s`.

*Why it gates everything:* if `C_freeze` is 80% of the step, then a cheaper march is irrelevant, §6 is
noise, and threading's speedup is set by how well the *pre-pass* threads. If it is 30%, priorities
invert. This is the single most valuable missing number in the document.

### T2 — the stability multiplier `Phi'`  *(gates §5)* — **DONE, §1g(i)–(iii)**

**Answer: the boundary is ~75–80 s, so the default `dt_fast = 150 s` is unstable (`Phi' = −2.74`);
the map is differentiable; and the CAS↔atmosphere conductance carries 99% of the lag gain.**

Protocol, for reproduction: on the high-LAI sunlit fixture, at a
fixed step, perturb `bio%cas%can_enthalpy` by δ **before** `column_fast_step` (so `build_column_frozen`
sees the perturbed state — that is what makes it the *lagged* map), advance one `dt_fast`, and form

```
Phi'  =  ( H_{n+1}(H_n + delta) - H_{n+1}(H_n) ) / delta
```

Repeat at `dt_fast` ∈ {100, 150, 225, 450, 900} s and test §1f predictions 1–2. Then **decompose `L`**
by repeating with each frozen coefficient (`gah`, `g_tr_f`, `h_coeff_f`, `t_emit`, `f_wet`) held at its
*unperturbed* value — that isolates each one's contribution to the lag gain, which is the sufficiency
question §2a says has never been asked.

Deliverables: `Phi'(dt_fast)`; the fitted `L`; a per-coefficient decomposition of `L`; and the
predicted stabilising under-relaxation factor `θ_crit`.

**Secondary deliverable worth keeping permanently:** a `test_cas_stability` fixture asserting
`|Phi'| < 1` at the default `dt_fast`. The project currently has **no detector** for this failure
mode (§1e) — the ledgers, the embedded estimator and the Newton convergence test are all blind. A
regression test on the multiplier is the missing instrument.

### T3 — the accuracy baseline at a stable cadence — **DONE, §1h**

**Answer: with the refresh on, `ark` and `rk45` agree to 1.3 mK at 150 s; frozen, the "gap" is 0.49 K
and is an artefact of the oscillation. The freeze is ~97% of the remaining error, so integrator
tolerance is irrelevant. Refreshed `ark` converges at roughly second order; the frozen scheme has no
order across its stability boundary.** Protocol, for reproduction:

`ark` vs `rk45` at `dt_fast` ∈ {150, 100, 75, 50} s with `rtol_all`/`atol_scale` tightened on the
reference, scoring CAS temperature, soil-surface temperature, GPP and ET. There is currently **no**
valid quantitative statement of how far apart the two schemes are, or of what any `dt_fast` costs in
accuracy — §2b deleted them all. Guard the reference cell as `numerics_sweep.py` already does (the
reference `dt_fast` must divide both `dt_slow` and 3600 s), and add the `parity_fidelity.py`
record-count / duplicate-key guard.

### T4 — the sufficiency question — **DONE, §1g(ii)**, and it answers differently than expected

Asked as a *derivative* decomposition rather than an amplitude sweep: build the frozen bundle at the
perturbed state, then overwrite one coefficient with its unperturbed value. **`gah`/`gaw`/`gac`:
`Phi'` −23.2 → +0.80. Everything else: ≤0.7%.** The earlier "pinning" experiments pinned `ustar`,
which leaves the stability-dependent transfer factor `temp1` free and therefore removes only about
half of `gah`'s sensitivity — which is why they concluded the aerodynamics merely amplifies.

**Still worth running as an amplitude confirmation** once N2 exists: the ED2 arm (refresh everything
per stage) as an upper bound, and a `dt_fast` amplitude scan with only `gah` refreshed, to confirm
that a multiplier below 1 actually produces a smooth trace rather than a differently-shaped one.

### T5 — the stability threshold where it is worst — **DONE (a), §1h**

**Answer: four of five stands are unstable at the default `dt_fast = 150 s`, and the worst is a
MID-HEIGHT canopy (10 m, LAI 3, `Phi' = −8.2`), not a regenerating one — prediction #4 was backwards.
The frozen scheme would need `dt_fast ≈ 50 s` to be safe there.** This is the correctness issue the
item was written to catch, and it outranks every optimisation in Part II: it means the *current*
default is unsafe, not merely slow. With the refresh, every stand is stable or marginal at every
cadence measured.

(b) — the tissue-store contribution — is **not run and is now moot**: with N2a the canopy air is
stable at every cadence measured, so separating the store's ~50% capacity contribution no longer
gates anything. Re-open only if N2a is rejected.

### T6 — evidence-base requirements for any cell reported

One site, one month per cell, one forcing year, no confidence intervals is the current standard and it
is too weak to rank anything within ~2×. Any claim promoted out of Part I carries: a dry-down window,
a cold season with snow, and a wet saturated window, in addition to the summer/winter pair; both
compilers; and `[output].numerics = true`.

### T7 — the ground-flux decomposition — **DONE, §1i.1, and it REFUTED the hypothesis**

Was: "apply §1g's decomposition to the ground surface; the soil-T error at 900 s should follow the
frozen `ggnet`/`soil_evap`." Ran across three regimes (6× range of ground forcing). The 0.70 K figure
does not reproduce (0.054–0.117 K), and the ground forcing varies 6× while the error does not. **No
ground refresh is built.** The measurement also produced the reusable rule in §1i.1: *a coefficient
moving a lot is not evidence that freezing it costs anything* — measure the error, not the swing.

### T8 — the full-state stability + convergence checklist — **DONE, §1i.2–§1i.4**

Was implicit in every previous `Phi'` caveat and never scheduled. Measures the **six-state outer-map
Jacobian** (canopy-air enthalpy / humidity / CO₂, top soil energy, top soil water, leaf water) and the
`dt_fast` convergence of both the states and the **fluxes** (GPP, ET, NEE). Outcome: `rho < 1` at every
cadence (discharging T2's standing caveat), all fluxes converged to 0.05% at 900 s, and `psi_leaf`
isolated as the single non-converged variable — with its error attributed to inheritance from the
canopy air rather than to the hydraulics.

**Why this belonged in the program from the start.** A stability-only assessment on one diagonal would
have declared the integrator healthy and stopped. The amplifying-but-stable leaf-water mode
(`∂(leaf_water)/∂(cas_shv)` = −3.75, loop gain 0.09) is invisible to a spectral radius *and* invisible
to a single-state `Phi'`. **Rank by both: `rho` for stability, per-state convergence for accuracy.**

## 4. How a candidate integrator or stabilisation is judged

A Part II option may become the default only when **all** hold:

1. It raises the measured stable `dt_fast` by ≥ 2× on **both** the high-LAI sunlit fixture and the
   young-stand fixture (T5), at CAS-T p2p ≤ 1 K — or, if it is a cost option rather than a stability
   option, it leaves `Phi'` unchanged within measurement error.
2. Net fast-loop cost per simulated day falls against today's 150 s configuration on the T1 profile.
   A stability win that costs more per step than it saves in steps is not a win.
3. **It does not change the limit.** The T3 refinement curve still converges, and to the same answer
   within the stated tolerance. A damping scheme that converges somewhere else is a *different model*
   — the lesson `split` cost this project, and it must not be repeated silently. This is the criterion
   most likely to fail for N1.
4. Every whole-column ledger still closes to machine precision with snow, condensation and the tissue
   store active, `debug_error = true`, on both schemes — **and** a consecutive-step CAS-T trace is
   inspected, because conservation is not stability.
5. Full `ctest` green on ifx Release, ifx Debug and nvfortran multicore (the issue-#7 rule: a green ifx
   run is not sufficient).
6. T6's scenario coverage.

A 30-year Ithaca run at the shipped configuration completes, conserves, and its AGB trajectory is
explicable. It is a **new baseline, not a comparison**: `runs/ithaca_ark30/` predates the recycle-phase
fix (PR #69), the snowfall routing fix, PR #88 and the `dt_fast` change.

---

# PART II — NUMERICAL IMPROVEMENT AND MULTI-CORE

## 5. Raising the stable `dt_fast` — the lever that cuts both cost terms

T2/T4 have run, so this section is no longer a menu of guesses — it has one indicated fix and a set
of alternatives that the measurement demoted.

> **SECTION STATUS, 2026-08-02.** Its purpose is **served**. `dt_fast` was raised 150 s → 900 s, the
> stability bound is gone (§1i.2), and the enabling fix (N2a) merged as PR #90. What is left here is:
> **N2b** (leaf-water potential — the sole non-converged variable, and correctness debt rather than a
> production blocker since PR #90 unwired it from carbon); and **N5** (adaptive freeze cadence), the
> only remaining item in this section with a plausible *efficiency* case. Everything else below is
> either landed, retired, or demoted. **The next real cost lever is §7, not §5.**

**N2 — Unfreeze `gah`/`gaw`/`gac`. THE headline item.** *(indicated by §1g(ii)+(iv); do this first)*

The CAS↔atmosphere conductance carries 99% of the lag gain and costs 2% of the pre-pass. Refreshing
it — and nothing else — is predicted to take `Phi'` from −23.2 to +0.80 at `dt_fast = 900 s` for ~3%
more work per sub-step. If that holds end to end it is up to a **6× cut in fast-loop cost**, which is
larger than everything else in this document combined.

> ### N2a IS IMPLEMENTED AND MEASURED — 2026-07-31. It works, and by more than predicted.
>
> *Read this block as a chronological record. The `aero_refresh` flag it describes was a measurement
> scaffold and **no longer exists in the code** — it was deleted before merge (see "LANDED" below), so
> do not go looking for it. What ships is the unconditional refresh, gated only by
> `surface_frozen_t%mo_live`, which is a populated-inputs marker and not a user switch.*
>
> `[ccfg].aero_refresh` (default **0 = frozen**, bit-identical; `1` = refresh per ESDIRK stage).
> `refresh_cas_conductances` in `meds_fast_ark.f90` re-solves **only** `mo_surface_layer` — the bulk
> MO solve — at the stage's own `(T_cas, q_cas)`, leaving every per-cohort quantity frozen. The
> refreshed `gah/gaw/gac` are the *local* variables `column_be_stage` already used for the BE commit,
> the Newton and the boundary-flux ledger, so "one flux, both sides" holds with no extra work.
>
> **Stability — the oscillation is gone at every `dt_fast` tested.**
>
> | `dt_fast` | 25 s | 50 s | 75 s | 100 s | 150 s | 225 s | 450 s | 900 s |
> |---|---|---|---|---|---|---|---|---|
> | `Phi'` frozen | +0.36 | −0.26 | −0.88 | −1.49 | −2.74 | −4.64 | −10.6 | −23.2 |
> | `Phi'` refreshed | +0.35 | −0.21 | −0.48 | −0.74 | **−0.79** | **−0.67** | **−0.42** | **−0.14** |
> | CAS-T p2p frozen [K] | 0.12 | 0.03 | 0.25 | 0.73 | 2.01 | 4.04 | 6.60 | 7.69 |
> | CAS-T p2p refreshed [K] | 0.12 | 0.08 | 0.04 | 0.02 | **0.012** | **0.015** | **0.032** | **0.100** |
>
> `|Phi'| < 1` everywhere, and it *decreases* with `dt_fast` — the residual lag stops mattering as the
> storage term shrinks. Peak-to-peak at 900 s falls **77×**, from 7.7 K to 0.10 K.
>
> **Production economics — a full 24 h day, 5 cohorts, scored against `aero_refresh=1 @ 25 s`:**
>
> | cell | wall | RMSE CAS-T | RMSE soil-T | sub-steps | rejections |
> |---|---|---|---|---|---|
> | **frozen @ 150 s (today's default)** | 24.8 ms | 0.390 K | 0.079 K | 931 | 4 |
> | refreshed @ 150 s | 22.5 ms | **0.017 K** | **0.063 K** | 576 | **0** |
> | refreshed @ 225 s | 18.2 ms | 0.036 K | 0.145 K | 664 | 1 |
> | refreshed @ 450 s | 10.1 ms | 0.085 K | 0.372 K | 348 | 4 |
> | refreshed @ 900 s | **6.4 ms** | 0.160 K | 0.701 K | 239 | 39 |
>
> At **matched** `dt_fast = 150 s` the refresh is a strict improvement on every axis at once: 9%
> cheaper, 23× better on canopy-air temperature, better on soil temperature, and it removes the
> rejections. Against today's default, refreshed @ 900 s is **3.9× cheaper and 2.4× more accurate on
> CAS-T**; refreshed @ 450 s is 2.5× cheaper and 4.6× more accurate.
>
> **It does not move the limit (§4 criterion 3).** At `dt_fast = 50 s` the *frozen* scheme scores
> 0.0016 K against the refreshed reference — the two semi-discretisations converge to the same answer,
> so this is a consistency improvement, not a different model.
>
> **Conservation holds.** `n_fail = 0` across all six ledgers in every cell. Full `ctest` is green on
> ifx Release **and** nvfortran multicore both with the flag off *and* with it forced on (36/36 each).
>
> **The cost model needs correcting, in the plan's favour.** §1g(iv) predicted ~3% overhead from the
> refresh itself. Measured at a single midday probe it looked like +50%, but that was the *march*
> taking an extra sub-step at that one state. Over a whole day the refresh **reduces** integrator work
> at every `dt_fast ≥ 150 s` (576 vs 931 sub-steps at 150 s), because a non-oscillating column is an
> easier problem for the error controller. The refresh pays for itself before any `dt_fast` change.
>
> ~~**Soil temperature is now the binding constraint.** CAS-T error is nearly `dt_fast`-independent
> once the feedback is live, but soil-T error still grows (0.06 → 0.70 K over 150 → 900 s) because the
> *ground* fluxes remain frozen. So the recommendation is **`dt_fast` 225–450 s**, not 900 s, until
> something equivalent is done for the ground surface.~~
>
> **VOID — 2026-08-02, §1i.1.** The 0.70 K does not reproduce across a 6× sweep of ground forcing
> (measured 0.054–0.117 K, worst case 0.4% of a 29 K ground diurnal amplitude), and the ground-lag
> mechanism is inconsistent with the scaling. `dt_fast = 900 s` stands; the ground refresh is not
> built. Kept struck through rather than deleted because the recommendation it carried was live for
> two days and shaped the §1i measurement plan.
>
> ### RK45 WIRED — 2026-07-31. The two schemes now agree to 1%.
>
> The family asymmetry is closed. RK45's CAS tendency is built in `surface_derivs` (not
> `column_derivs`), so that is where its refresh lives; `stage_bnd` re-solves the same `pure` kernel
> at the same stage state so its boundary flux is charged at the conductance the tendency used.
>
> **ARK pays for one surface-layer solve per stage, not one per residual evaluation.** `surface_derivs`
> is called up to 24 times per stage by the Newton, to fill a CAS tendency ARK never reads — it commits
> the CAS through its own backward-Euler denominator. So `column_be_stage` refreshes once, writes the
> result into its stage-local `fs`, and **clears `fs%aero_refresh`**, after which `surface_derivs` reads
> the already-refreshed values. RK45, whose every RHS evaluation *is* one stage, keeps the flag set.
>
> **5 cohorts, 24 h day, scored against `ark`/refresh/25 s:**
>
> | scheme | refresh | `dt_fast` | wall | RMSE CAS-T | RMSE soil-T | sub-steps | rescue | `n_fail` |
> |---|---|---|---|---|---|---|---|---|
> | ark | off | 150 s | 39.8 ms | 0.390 K | 0.079 K | 931 | 0 | 0 |
> | ark | **on** | 150 s | 25.3 ms | **0.0171 K** | 0.063 K | 576 | 0 | 0 |
> | rk45 | off | 150 s | 32.6 ms | 0.375 K | 0.077 K | 576 | 0 | 0 |
> | rk45 | **on** | 150 s | 41.5 ms | **0.0169 K** | 0.061 K | 838 | 0 | 0 |
>
> **ARK and RK45 land within 1% of each other** (0.0171 vs 0.0169 K) — an implicit and an explicit
> method agreeing, which is the strongest available evidence that this is a consistency fix rather
> than a scheme-specific tweak. `n_fail = 0` on all six ledgers in every cell, which is also what
> validates the `stage_bnd` plumbing: had the ledger used a different `gah` than the tendency, whole-
> column energy would have failed immediately.
>
> **The refresh helps ARK's cost and hurts RK45's** — sub-steps 931 → 576 on ARK, 576 → 838 on RK45.
> Expected, and it sharpens the tier split: making `gah` state-dependent puts a strong negative
> feedback into the RHS, which an implicit method *absorbs* (the column stops oscillating, so the
> controller relaxes) and an explicit method must *resolve* (its stability limit tightens). At 900 s
> refreshed, `ark` is 2.3× cheaper than `rk45` and slightly more accurate. `rk45` is the accuracy
> baseline, not the production tier, so paying there to buy consistency is the right trade.
>
> **Bit-identity of the default path: verified, after initially losing it.** The first version reported
> the used conductances on `surface_tend_t` so `stage_bnd` could read them back — the more obviously
> safe design. Measured against a `git stash` of the pre-change tree it was **not** bit-identical:
> ~1e-12 relative, because three extra fields on that type perturb ifx's inlining/FMA decisions inside
> `surface_derivs`. Harmless numerically, but it silently shifted derived RMSEs by up to 8% in the
> cells where the frozen map is unstable and round-off gets amplified, and it would break the
> stash-and-`cmp` protocol the project verifies with. Removing the fields and re-solving in
> `stage_bnd` restores **bit-for-bit** identity against the pre-change tree, confirmed by `cmp` on a
> 288-sample trajectory across both schemes at three `dt_fast`.
>
> ### LANDED 2026-07-31 — the flag is deleted and the refresh is unconditional.
>
> User decision, on the T3/T5 evidence: **delete the switch, always refresh.** The frozen alternative
> is unstable at the production cadence on four of five stand heights (`Phi'` to −8.2), so it is
> known-wrong physics rather than a supported configuration — the same reasoning that deleted
> `snow_on` and the `with_theta` norm switch.
>
> `surface_frozen_t%mo_live` is what remains, and it is **not** that switch: it records whether the
> `mo_*` inputs are populated, and therefore whether a re-solve is possible at all. Three uses, none
> of them a user choice — default `.false.` ("use `gah/gaw/gac` as given", which is what a hand-built
> unit-test bundle means); `build_column_frozen` sets `.true.`; ARK's stage copy clears it after
> refreshing once, so the Newton's ≤24 residual evaluations per stage do not each re-solve.
>
> **That default direction was learned the hard way.** With the marker defaulting to "refresh",
> `test_column_derivs` — which builds a `surface_frozen_t` by hand and never populates the `mo_*`
> fields — called the surface-layer solve with zero roughness and zero wind and **hung**. The safe
> default is the one a bundle that knows nothing about this mechanism will get.
>
> `dt_fast` default moved 150 s → **900 s** with it, per the same decision. See the warning below.
>
> Validation: `ctest` 36/36 on ifx Release and nvfortran multicore.

> **MERGED as PR #90, 2026-08-01** (`feature/cas-conductance-refresh` → `main`, commit `8b7ee1b`), in
> three commits: the per-stage conductance refresh; the `wstress_nonstomatal` gate + `dt_fast = 900 s`
> (breaking — the TOML key is required); and the documentation rewrite.
>
> **A 50-year Ithaca pair was run with the limb on and off at 900 s and it does NOT support 900 s** —
> a distinction easy to lose. The two arms differ by 0.39% in mean GPP and 0.56% in final AGB, which
> says only that this mesic site rarely drives `psi_leaf` far below `psi_open`. The `dt_fast` argument
> rests entirely on the single-column probes; **no site-level `dt_fast` scan exists.**

### The carbon bias at long `dt_fast` — measured while landing the above, and it constrains the default

Moving `dt_fast` to 900 s costs far more in carbon than in temperature. Against a 12.5 s reference on
the high-LAI sunlit stand:

| `dt_fast` | 150 s | 225 s | 300 s | 450 s | 900 s |
|---|---|---|---|---|---|
| GPP | −3.8% | −8.4% | −12.6% | −19.8% | **−33.1%** |
| ET | −2.6% | −5.8% | −8.7% | −14.5% | **−23.8%** |
| CAS-T RMSE | 0.017 K | 0.037 K | 0.055 K | 0.085 K | 0.16 K |
| soil-T RMSE † | 0.06 K | 0.14 K | 0.22 K | 0.37 K | 0.70 K |

† **The soil-T row did not reproduce (§1i.1).** A three-regime re-measurement spanning a 6× range of
ground forcing put the 900 s figure at 0.054–0.117 K, not 0.70 K. The original fixture could not be
recovered, so this row is flagged rather than deleted — treat it as unconfirmed. The GPP/ET/CAS-T rows
were independently re-derived and stand.

**Attributed in two steps, and the second changes the conclusion.**

*Step 1 — it is the plant hydraulics.* Freezing the plant-hydraulics store (`mask%hydraulics`) removes
the GPP bias entirely — **−33.1% → +0.93% at 900 s, and flat in `dt_fast`**.

*Step 2 — it is not a quadrature bug, and the stress limb is an amplifier rather than the source.* GPP
is accumulated as `Σ gpp(state n)·dt`, a left rectangle in the state, so an averaging artefact was the
obvious suspect. Disabling the **non-stomatal (capacity) water-stress limb** — a linear ramp
`beta = (ψ − ψ_close)/(ψ_open − ψ_close)` on Vcmax/Jmax/TPU — settles both questions at once:

| `dt_fast` | 12.5 s | 150 s | 300 s | 450 s | 900 s |
|---|---|---|---|---|---|
| GPP, limb on | 1.000 | 0.962 | 0.874 | 0.802 | **0.667** |
| GPP, limb **off** | 1.000 | 1.000 | 1.000 | 1.000 | **0.9995** |
| ET, limb off | 1.000 | 1.000 | 1.000 | 0.993 | 0.988 |
| daytime-mean ψ_leaf [MPa] | −0.23 | −0.50 | −0.72 | −0.88 | **−1.19** |
| max per-step \|Δψ\| [MPa] | 0.0001 | 0.010 | 0.034 | 0.071 | **0.84** |
| fraction of daylight steps with beta < 1 | 0.00 | 0.61 | 0.77 | 0.83 | 0.91 |
| fraction with beta clamped at 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |

So: **the quadrature is sound** (0.05% at 900 s with the limb off), **ψ_leaf itself is badly
non-converged** (−0.23 → −1.19 MPa), and **the limb converts that linearly into carbon** — slope
~0.5 MPa⁻¹ on the shipped PFT file, so ~1 MPa of ψ error ≈ 50% of Vcmax. The `beta = 0` clamp is never
reached, so this is straight linear amplification, not a saturation/Jensen artefact.

**Decision (issue #47, user, 2026-07-31): the capacity limb is OFF by default.**
`[leaf_physiology].wstress_nonstomatal`, default `.false.`. It is rarely measured directly and its two
parameters are weakly constrained, so wiring an unconverged ψ into carbon through it is not a trade
worth making. The stomatal limb (ψ_soil → `g1`/λ) is better constrained and stays on. The kernel keeps
the term and `test_leaf_physiology` still exercises it explicitly, plus a new assertion that ψ_leaf is
**inert** in the default configuration.

**Consequence for the `dt_fast = 900 s` default: it is now defensible.** With the limb off, 900 s costs
0.05% of GPP and 1.2% of ET. The `meds_config` warning is conditional on the limb accordingly.

**But the ψ error is unfixed — it is only unwired.** Anything keyed to leaf water potential still
inherits it, and it will come back the moment ψ feeds anything else (hydraulic mortality, ψ-driven
allocation, or the limb being re-enabled for a study). That is N2b, and it is still the top item.

### N2b — the leaf-water-potential lag. THE remaining accuracy item, and §1i.4 changed what it is.

`psi_leaf` is diagnosed once per `dt_fast` from `leaf_water_mass` at state `n` (`column_prepass`) and
feeds `g_sw` through `leaf_gas_exchange_batch`. It is the **only** state that fails to converge at
900 s (§1i.3: 0.76 MPa, against `T_cas` 0.053 K and every flux inside 0.05%).

**Its verification step has now run, and it inverted the diagnosis.** This section used to read "ψ is
badly non-converged *because* the water-mass update is an explicit step with frozen sapflow/uptake" —
i.e. an autonomous integration error of the leaf water ODE. The Jacobian (§1i.4) says otherwise:

- the two largest entries in the whole six-state matrix are both **into** `leaf_water` **from** the
  canopy air — `∂(leaf_water)/∂(cas_shv)` = **−3.75** and `∂(leaf_water)/∂(cas_enthalpy)` = **+1.03**
  at 900 s, both growing superlinearly with `dt_fast`;
- `leaf_water`'s **own** diagonal is ≈0.09 — it retains almost nothing of its own state across a step;
- the return path is weak (+0.023), so loop gain ≈0.09: it **amplifies without oscillating**, which is
  why §1i.2's spectral radius is blind to it.

So the error is **inherited from the canopy air and amplified ≈4×**, not generated in the hydraulics.
`psi_leaf` is a *fast, near-algebraic variable being integrated as a prognostic one under frozen
coefficients* — the same defect category as the canopy air before N2a, one seam over.

> **SUPERSEDED 2026-08-01 — the mechanism is now known exactly, and it is neither "inherited and
> amplified" nor a coefficient problem. It is a FLUX INCONSISTENCY, provable by algebra.** See
> §N2b-RESOLVED immediately below. The Jacobian reading above is not *wrong* — the dependence on
> `cas_shv`/`cas_enthalpy` is real, and it is real *because transpiration depends on them* — but the
> prescriptions it motivated (tangent-linear `g_sw`, algebraic leaf store, "advance ψ within the step")
> all target the wrong object. None of them were needed. The three bullets below are **retired**; they
> are kept only because §1f scores predictions against them.

**What that implies for the fix.** The target is the **coefficient**, not the integrator: a better ODE
solver for `leaf_water_mass` addresses a term whose diagonal is already ≈0.09. Candidates, in order:

- **tangent-linear `g_sw`** — evaluate `∂g_sw/∂ψ` (and `∂ψ/∂W`) once per `dt_fast` alongside the kernel
  call already happening, then correct at stage points. The mass update is a closed-form Euler step, so
  ψ at a stage point is free. Re-solving the *whole* leaf kernel per stage is not affordable (~89% of
  the pre-pass at 30 cohorts) — this is the affordable analogue of what N2a did for `gah`.
- **treat the leaf store algebraically** — a diagonal of ≈0.09 says it is already nearly slaved within
  one step, so solving it as an algebraic closure (as `veg_energy_diagnostic` does for tissue heat) may
  be both cheaper and more accurate than integrating it. This is N4's structure applied to a store that
  measurement says is genuinely fast.
- **advance ψ within the step** so the value the next step inherits is a step-average, not an endpoint.

**Open first question, before any of the above.** The leaf-water diagonal is nearly flat in `dt_fast`
(0.085 → 0.099 over 150 → 900 s). A store that relaxes with a fixed `tau` would decay as `exp(−dt/tau)`
— 0.085 at 150 s implies `tau ≈ 60 s`, which would put the 900 s diagonal at ~1e-7, not 0.099. The
measured flatness is **unexplained**, and it is the thing to understand before choosing between
"refresh the coefficient" and "make it algebraic": the two prescriptions differ precisely in what they
assume about that relaxation.

**Priority note.** With `wstress_nonstomatal` defaulted off (PR #90), nothing currently consumes
`psi_leaf` in a way that reaches carbon — GPP is converged to 0.05% at 900 s. So this is **correctness
debt, not a production blocker**: it returns the moment ψ feeds hydraulic mortality, ψ-driven
allocation, or the limb is re-enabled for a study. Rank it against §7 accordingly.

### N2b-RESOLVED — MEASURED 2026-08-01. The mechanism, the fix that landed, and what it exposed.

**The mechanism, by algebra rather than by inference.** `solve_plant_water` defines its own reported
fluxes *from* its own storage change (`meds_plant_hydraulics.f90`):

```fortran
flux%sapflow     = dw_l/dt + e_transp          ! dw_l = W(psi_end) - W(psi_start)
flux%root_uptake = (dw_l + dw_w)/dt + e_transp
```

Therefore `W + dt*(sapflow − transp)` reproduces the kernel's matrix-exponential endpoint **exactly, at
any dt** — *provided the transpiration pricing the sapflow is the transpiration actually debited*. It
was not. The pre-pass priced `sapflow_frozen` at `transp_pp` (state-*n* full demand); `advance_water_
mass_full` debited `transp_c_bw` (the b-weighted realised demand, chosen deliberately so the mass debit
matches the CAS vapour credit). The leaf error is therefore **exactly `dt*(transp_pp − transp_bw)`** —
a pure flux inconsistency, growing *linearly* in dt, never converging away within a step.

Independent check: 8% transpiration mismatch × 8.14e-4 kg/plant/s × 900 s ÷ ~1 kg pool = **5.9%**
pool drawdown predicted, **5.7%** measured. This is the defect class in
`project_meds_frozen_flux_defect_class`, not a stiffness or coefficient problem. The claim in the
superseded text that "the leaf is a pure integrator with no restoring force" is **wrong** — the
restoring force is present, inside `sapflow_frozen`.

**The fix (LANDED, in tree, ifx 36/36 + nvfortran multicore 36/36).** Re-solve the kernel in
`advance_water_mass_full` on the *realised* transpiration, on the SAME Category-0 frozen conductances
and soil potentials as the pre-pass — one semi-discretisation, not a second linearisation point. Cost
is one extra hydraulics solve per step, measured at **+0.6–7% of a `column_fast_step`**.

3 h midday probe, common **dt = 5 s** reference (`<scratchpad>/meas_psi2.f90`):

| dt_fast | ψ_leaf err before | ψ_leaf err after | ψ_wood err before | ψ_wood err after |
|---|---|---|---|---|
| 25 s  | 1.44e-2 | 3.07e-4 | 4.21e-4 | 6.57e-4 |
| 450 s | 4.98e-1 | 5.27e-4 | 2.27e-3 | 4.66e-2 |
| 900 s | **1.454** | **4.63e-3** | 4.47e-3 | **1.743e-1** |

Both versions agree at 5 s to 8e-4 MPa, so this changes the **discretisation, not the dt→0 limit**.
`T_cas` differs by 6e-5 K at 25 s and 0.007 K at 900 s — also vanishing with dt, as a consistent scheme
must. (This also **voids** the earlier "T_cas regression that survives to dt = 12.5 s" blocker: that
0.032 K belonged to a hand-rolled leaf relaxation, not to correcting the transpiration.)
Regression guard: `test_psi_dt_convergence` in `test/test_column_ark.f90` — its window **must** start
at solar noon or there is no transpiration signal and it passes on broken code.

**What it exposed: the error RELOCATES to the wood, and no plant-side BC can stop it.** `uptake_frozen`
is the number the soil column already committed as its root sink, so the corrector cannot re-price it
without breaking whole-column water. A specified-flux (**Neumann**) root BC was fully built and
measured to test whether the BC choice could fix this — it cannot:

| root BC | ψ_leaf err @900 s | ψ_wood err @900 s | plant total @900 s |
|---|---|---|---|
| potential (Dirichlet, **landed**) | **4.6e-3** | 1.74e-1 | 1.9139 |
| specified-flux (Neumann, **not landed**) | 1.61e-1 | 1.60e-1 | 1.9139 |
| 5 s reference | — | — | 1.9392 |

The plant **total** is identical under both BCs and ~1.3% low under both: the deficit is set by
`uptake_frozen` being priced at the pre-pass transpiration, and the root BC only chooses **where the
deficit is parked**. Park it in the wood — ψ_leaf drives the stomatal water-stress feedback, while
ψ_wood at these potentials is far from `wood_psi50` and drives almost no PLC. Patch kept at
`<scratchpad>/neumann_root_bc.patch`. Implementation note worth keeping: the Neumann branch needs its
own exact solution because with the flux fixed at *both* ends the 2×2 `M` is **singular**
(`det M = keff²/(cl·cw) − keff²/(cl·cw) = 0`); the gradient `g = ψ_w − ψ_l − grav` obeys a scalar
linear ODE and the total obeys `d(W_l+W_w)/dt = U − T`, which integrate exactly.

Second implementation note: the kernel's reported `flux%root_uptake` equals a specified `u_bc` only to
**linearisation** error — sub-steps conserve in `dW = C·dψ`, but `dw` is re-derived from the *nonlinear*
PV curve at the endpoints. That residual breaks the whole-column ledger (`test_fast_loop` catches it).
Always take the committed number verbatim, never the kernel's echo of it.

### N2d — the uptake seam. THE remaining hydraulic item. Design fixed, staged, IN PROGRESS.

Closing the 1.3% deficit requires the **soil** side re-priced, i.e. the soil must see the realised
transpiration. The right structure is a **true Godunov split** — vegetation/CAS first at θⁿ, then soil
water advanced *once* with the realised uptake as its sink:

```
Act 1 (θⁿ):  snow → interception → surface_derivs → soil_evap(θⁿ) → supply throttle f(θⁿ)
                                 → hydraulics pre-pass (seeds the stages' mass ODE + qloss)
Act 2:       ESDIRK stages — CAS + soil heat, θ frozen at θⁿ      [unchanged]
Act 3:       advance_water_mass_full (transp_bw) → REALISED uptake
             column_hydrology_flux with THAT uptake as its sink → θ¹ committed
```

The seam then closes by construction — no booking, no post-hoc edit of a committed θ, no over-drain
risk. Expected: ψ_wood 1.74e-1 → ~1.6e-3 (measured by forcing the corrector's uptake as a diagnostic),
ψ_leaf unchanged at 4.6e-3; max-norm ~35×.

**The ordering constraint that motivated today's design has been measured, and it dissolves.** The only
reason the soil must currently run *first* is the supply throttle
`scale = min(1, hflux%uptake_total / total_uptake_b)`. Measured with public `column_hydrology_flux`
(`<scratchpad>/meas_scale.f90`, dt_fast 900 s, θ_res 0.078):

| θ | ψ_soil | scale @0.25× demand | @site | @2× | @4× |
|---|---|---|---|---|---|
| 0.35–0.20 | −0.29 … −1.78 MPa | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| 0.15 | −4.69 | 0.99115 | 0.99105 | 0.99091 | 0.99063 |
| 0.12 | −12.35 | 0.93982 | 0.93940 | 0.93878 | 0.93759 |
| 0.10 | −39.25 | 0.75943 | 0.75737 | 0.75438 | 0.74859 |
| 0.085 | −303 | 0.0 | 0.0 | 0.0 | 0.0 |

The throttle is **inactive above θ ≈ 0.20**, and it is **near-independent of demand** (17× demand range
moves it 1.3%) — so the limitation is a *proportional* wilting factor `f(ψ_k)` on the sink, not a
ceiling, and `scale = Σ_k root_frac_k · f(ψ_k)` is a function of the soil profile **alone**. It can
therefore be evaluated from θⁿ *before* the plant runs, no less accurately than today (today's solve
also starts from θⁿ).

**The real obstacle, found while scoping: a circular dependency in the ENERGY terms.** The soil-*heat*
stage consumes hydrology outputs — `fro%clip_enth`, `fro%floor_enth`, `fro%w_flux_frozen` all pair
enthalpy with water movement inside `column_be_stage`. So soil heat (in the stages) needs soil water's
fluxes, while soil water (after the stages) needs the plant's uptake. Committing a θ from a *different*
solve than the one that supplied the enthalpy pairing is precisely
`project_meds_frozen_flux_defect_class`, so this cannot be waved through.

**Staged plan.**

- **P0 — the transpiration corrector. ✅ DONE**, §N2b-RESOLVED, in tree, 36/36 both back ends.
- **P1 — RK45 state-dependent root sink. ❌ BUILT, MEASURED, REVERTED 2026-08-01. Do not retry on an
  explicit path.** The idea: `column_derivs` reads `fro%uptake` as a frozen constant, but RK45
  integrates its **own** θ in-stage, so it could take an *algebraic* sink
  `rhizo_k·(ψ_soil_k − ψ_wood(stage))` with the identical expression feeding `d_wood_water_mass` —
  "one flux, both sides" per stage, algebraic rather than a sub-solve. It was implemented and it
  **broke `test_column_dynamics`** (RK45-vs-split snow melt disagreed by 0.47 against a 10% band).
  That was not a tolerance problem. Measuring the mode it introduces
  (`<scratchpad>/meas_tauw.f90`, `C_wood = 0.488 kg/plant/MPa`):

  | θ | rhizo_total [kg/plant/s/MPa] | τ_w = C/rhizo | dt_fast/τ_w @900 s |
  |---|---|---|---|
  | 0.35 | 1.17e+0 | **0.42 s** | **2162** |
  | 0.25 | 5.35e-2 | 9.12 s | 98.7 |
  | 0.15 | 2.48e-4 | 1963 s | 0.5 |
  | 0.10 | 1.17e-6 | 4.15e5 s | 0.002 |

  A Darcy sink gives the wood node its **own** relaxation time, and in wet soil it is **sub-second**.
  An explicit stage cannot carry that at 900 s: RK45 would shrink to ~1 s steps (a ~900× blow-up) and
  the RK4/IMEX oracle would simply be unstable. This **vindicates the original design comment** in
  `column_derivs` — *"mass adds no stiff mode because its inflow is frozen and its outflow moves at
  the CAS timescale"* — which is a load-bearing statement, not a description.

  **What it implies for P2, and it sharpens it:** the stiff wood↔soil relaxation must stay **inside**
  the matrix-exponential kernel, which handles it with its own adaptive sub-stepping. Only the
  kernel's *output* may cross the seam. P2 is unaffected — its uptake comes from the corrector's
  kernel call, not from an explicit stage — but any future "refresh the uptake in-stage" proposal on
  an explicit tableau is dead on arrival, and this table is why.

  Secondary consequence: RK45 as the **accuracy oracle** needs no fix here. The oracle is exercised at
  small `dt_fast`, where `transp_pp → transp_bw` and the frozen-uptake deficit vanishes by
  construction. P1 was justified as a prerequisite for validating P2; that justification was wrong.
- **P2 — ARK: the Act-1/Act-3 reordering above. NEXT.** The blocker is the **enthalpy-pairing
  circularity**: water carries heat, so the soil-*heat* stage needs the water movements
  (`fro%clip_enth`, `fro%floor_enth`, `fro%w_flux_frozen`) that the soil-*water* solve produces — but
  P2 moves that solve to *after* the stages. Soil heat needs soil water; soil water needs the plant's
  uptake; the uptake needs the stages. Two resolutions:

  - **Option A — two solves.** Keep a predictor solve up front purely to supply the stages'
    heat-transport fluxes; run a corrector at the end with the realised uptake and commit *its* θ.
    Cost +5–15% (§ cost table below). **Requires the paired enthalpy to be re-sourced from the
    corrector as well** — otherwise heat moves per one solve and water per another, which is
    `project_meds_frozen_flux_defect_class` verbatim.
  - **Option B — defer the advection. ← RECOMMENDED.** Take water-borne enthalpy out of the stages
    entirely: stages do soil-heat *conduction* only, and all advective heat is applied once at the end
    alongside the committed θ. One solve, one set of fluxes, mismatch structurally impossible. It
    mirrors how soil *water* is already treated. Rationale for preferring it over A: the frozen-flux
    defect has bitten this codebase three separate times in this work alone, and A's correctness rests
    on wiring the pairing correctly by hand, while B removes the failure mode by construction.

  **Do T9 (below) before choosing.** B changes a soil-heat stage that currently has passing
  conservation tests, and if the advective term is large at 900 s then deferring it is itself an
  approximation needing justification. If it is negligible, B is straightforward.

  Also needed for P2 regardless of A/B: `ground_evaporation` exported from `meds_soil_water` (only
  `column_hydrology_flux`, `soil_water_time_deriv`, `soil_water_step_implicit` are public today), a
  decision on `q_top` (conductivity-limited infiltration is computed *inside* the Richards solve, so
  it is not available pre-stages), and the soil-water `bf` ledger terms rewired.

- **P3 — verify the θⁿ throttle** against the solve-derived one at θ ≈ 0.10, where `scale` = 0.76 and
  the discrepancy is largest.

### T9 — the advective soil-heat term — **DONE 2026-08-01. Option B is CONFIRMED.**

Two results, and the second is structural.

**(i) The advective Courant number is ≪ 1 at 900 s** (`<scratchpad>/meas_advect.f90`, reading
`fro%w_flux_frozen` straight out of the public `build_column_frozen` — no instrumentation):

| regime | max Courant `|w|·dt/dz` | max `clip−floor` source |
|---|---|---|
| dry, θ 0.25 | 8e-5 | 0 |
| 10 mm/h rain, θ 0.25 | **0.087** | 2.32e3 W/m² |
| 10 mm/h onto θ 0.38 (near sat) | 0.013 | 2.02e3 W/m² |

Explicit application of the interior advection is stable everywhere at production cadence.
*(A `dT_src` column in that probe reports 11–13 K for the clip term; it is MISLEADING and must not be
quoted as a temperature — it treats the saturation clip as a pure heat source at fixed heat capacity,
whereas the clip removes mass and energy TOGETHER and the temperature barely moves. What the number
does show is that the clip term is LARGE in magnitude.)*

**(ii) The advection is ALREADY explicit — it was never inside the implicit operator.** Reading
`soil_energy_step_implicit`: `soil_heat_be_solve` takes `kappa, c_eff, q_src, g_top, geothermal` and
solves **conduction only**. The water-enthalpy fluxes `qwf(0..n)` are then formed from the
*post-solve* `t_new` (upwind) and applied in the conservative update
`col%soil_energy(k) += dt/dz*((hf(k)−hf(k−1)) + (qwf(k)−qwf(k−1)))`. So `w_flux` never enters the
tridiagonal system.

**This removes the main objection to Option B.** Deferring the advection does not convert an implicit
operator into an explicit one — it moves an *already-explicit* post-solve correction from per-stage to
once-per-step, at Courant ≤ 0.087. Combined with (i), and with the fact that Option A's
predictor/corrector enthalpy mismatch would be **2.3e3 W/m² in the rain regimes** — not a small
residual — **Option B is the choice.**

**Scope consequence, and it shrinks P2.** If *all* advective channels are deferred together
(`w_flux` interior, `w_flux_top`/`qwf(0)`, `w_flux_bot`, `clip_enth`, `floor_enth`), then the ESDIRK
stages need **no hydrology input for energy at all** — they reduce to conduction + surface heat flux +
the transpiration root-heat sink (`coh_qsoil`, `qloss`). The only remaining pre-stage dependency on
the soil solve is `fro%surf%soil_evap` (the CAS vapour source), which is an algebraic
`ground_evaporation` evaluation at θⁿ and does not need the Richards solve. `q_top` stops being a
pre-stage requirement because the top-face water enthalpy is deferred with everything else.

So P2 becomes: **export `ground_evaporation` → defer all advective enthalpy out of the stages → move
the single Richards solve after `advance_water_mass_full` and feed it the realised uptake.** One soil
solve, not two.

#### P2 implementation steps (mapped 2026-08-01; step 1 DONE)

1. **✅ DONE — the soil-evaporation seam.** `ground_evap_from_state(theta1, params, forcing, opts, dt)`
   is now public in `meds_soil_water`, and `column_hydrology_flux` **delegates to it**, so there is one
   implementation and the driver cannot drift from the solve. Verified behaviour-neutral (36/36).
   This is **exact, not an approximation**: the alpha_soil/DSL formula and the storage cap both read
   `theta0(1)` — the entry moisture — so the value never depended on the Richards transport at all.

2. **NEXT — factor the advective application out of `soil_energy_step_implicit`.** The `qwf(0..n)`
   block plus its share of the conservative update should become a public
   `apply_water_enthalpy_advection(col, forcing, soil, dt, ...)`, called by the kernel itself (so the
   upwind rule stays single-sourced) and callable by the driver afterwards. **This is surgery on a
   conservation-critical kernel — do it as its own commit, verified bit-identical before anything else
   moves.** The natural test is that the existing suite stays green with the kernel merely delegating.

3. **Then — strip the advective inputs from `column_be_stage`:** `fro%clip_enth`/`fro%floor_enth` out
   of `eforc%root_heat_sink`, `eforc%w_flux` → 0, `eforc%w_flux_top`/`w_flux_bot` → 0, and `e_drain`
   off `root_heat_sink(nsl)`. The stage's `bf` ledger loses `e_infil`, `e_floor`, `e_drain`, `e_clip`
   — **these move to the post-step application, they are not deleted**; the whole-column ledger must
   still close, which is what will catch a mistake here.

4. **Then — reorder** in `column_fast_step_ark`: `build_column_frozen` stops calling
   `column_hydrology_flux` (it computes `soil_evap` via step 1 and the supply throttle from θⁿ);
   after `advance_water_mass_full` the single Richards solve runs with the realised uptake, commits
   θ¹, and its advective enthalpy is applied via step 2.

5. **Verify:** `psi_wood` 1.74e-1 → ~1.6e-3 expected on `<scratchpad>/meas_psi2.f90`; whole-column
   water and energy ledgers still machine-precision; ifx + nvfortran multicore both 36/36.

### RK45 after P1's reversal — MEASURED 2026-08-01, and it needs a WARNING not a fix

`meds_fast_rk45` never calls `advance_water_mass_full`; it integrates plant water *inside* its own
stages via `column_derivs`. **So RK45 did not receive the P0 corrector, and still carries the full
pre-P0 defect.** Same 3 h midday probe, same 5 s reference:

| dt_fast | RK45 ψ_leaf err | ARK (post-P0) ψ_leaf err | RK45 ψ_wood err | ARK ψ_wood err |
|---|---|---|---|---|
| 450 s | 4.95e-1 | 5.27e-4 | 2.26e-3 | 4.66e-2 |
| 900 s | **9.551e-1** | **4.63e-3** | 2.39e-3 | 1.74e-1 |

RK45 now behaves like pre-fix ARK, and the two schemes have **mirror-image** error structure: RK45
parks everything in the leaf (its wood ODE is two frozen constants, so the wood converges); ARK parks
it in the wood. **ARK is therefore more accurate than its own accuracy oracle for `psi_leaf`.**

**RECOMMENDATION: leave RK45 alone; document the constraint.**

1. *It is still valid where it is used.* At 5 s the two agree to 8e-4 MPa. The oracle is exercised at
   small `dt_fast`, where `transp_pp → transp_bw` and the defect vanishes by construction.
2. *But it is a trap.* Anyone comparing ARK against RK45 at production cadence will see ~1 MPa of
   `psi_leaf` disagreement and conclude ARK is broken — when ARK is the correct one. **Do not use RK45
   as a `psi_leaf` reference above `dt_fast` ≈ 150 s.**
3. *Do not try to fix it in-stage.* The stiffness wall that killed P1 applies to the leaf too:
   `tau_leaf ≈ 15.5 s` against `dt = 900 s`, so an explicit stage needs `dt <~ 43 s`. The only way to
   give RK45 the benefit is to pull plant water out of its tableau and share ARK's operator-split
   corrector — which would cost the very thing that makes an oracle useful, being an **independent
   construction**. If both schemes share the same operator-split plant water, RK45 can no longer
   validate that choice. Keep them different; be explicit about the valid range.

**Cost context for P2** (`<scratchpad>/meas_cost.f90`, dt_fast 900 s): `column_hydrology_flux` is
8.6/7.6/7.8 µs at 1/3/20 cohorts against a 56/87/165 µs `column_fast_step` — i.e. **4.7–15.3%**. A
second (corrector) soil solve is affordable. Note this is *not* an argument for moving soil water into
the ARK stages: that was measured separately and the split error it would buy out is **0.011% of column
water dry, 0.018% under 10 mm/h rain** (θ(1) error 8.5e-5 → 2.0e-3 m³/m³), ~10× below sensor accuracy.
**Do not move soil water into the tableau; there is no prize.** Keep these two distinct — *soil-water
split error* is negligible, *uptake-freeze error* is the one that matters.

Three implementations of the `gah` refresh were scoped, cheapest first. **One landed; the other two
are retired — the `gah` seam is closed.** (Note: `N2b` unqualified means the *leaf-water* item above.
The two `gah` variants are suffixed `-gah` to end a naming collision that was live in this file.)

- **N2a — per-stage refresh. ✅ IMPLEMENTED AND MERGED (PR #90).** Re-solve `mo_surface_layer` inside
  `column_be_stage` at the stage's own CAS state, leaving every other frozen quantity alone. This is
  what ED2 does (`update_diagnostic_vars` → `canopy_turbulence8` at every RK stage). Deliberately
  *not* a full `canopy_aerodynamics` call — the per-cohort boundary layers stay frozen, since they
  carry ≤0.7% of the lag gain and are what makes a full refresh expensive.
- **N2b-gah — tangent-linear `gah`. ❌ RETIRED 2026-08-02, do not build.** Evaluate `∂gah/∂H` once per
  `dt_fast` and apply `gah(H) ≈ gah_n + (∂gah/∂H)(H − H_n)` at stage points. **Its entire attraction
  was being cheaper than N2a — and N2a measured cost-NEGATIVE** (931 → 576 sub-steps at 150 s, because
  a non-oscillating column is an easier problem for the error controller). There is no saving to
  capture by approximating something that already pays for itself. Worse, this file already recorded
  that `d ln gah/dT ≈ 2.2 K⁻¹` is a strong nonlinearity, so a linear model over a multi-K excursion
  would trade accuracy for a saving that does not exist.
- **N2c-gah — `gah` inside the Newton. ❌ RETIRED 2026-08-02, do not build.** Add `∂gah/∂H` to
  `jac_surface`, making ventilation implicit rather than merely refreshed. Its attraction was
  robustness — it cannot overshoot. **§1i.2 shows there is no overshoot to guard against:** `rho < 1`
  at every cadence and *decreasing* with `dt_fast`, and rejections already fell to zero at 150 s under
  N2a. It would add one `surface_derivs` evaluation per Newton iteration to fix a problem the model
  does not have.

> **The `gah`/CAS seam is closed.** Stability is solved (§1i.2, full Jacobian), the fix is cheaper than
> what it replaced (§1g), and both cheaper-or-safer alternatives are motivated by properties the
> measurements refute. Remaining efficiency levers are **not** on this seam: they are §7 (multi-core,
> untouched, the only lever costing no accuracy) and possibly N5 (adaptive freeze cadence). Reopen only
> if a stand geometry or wind regime outside the measured sweep puts `rho` back near 1.

**Watch items for all three.** `gah` appears in the CAS commit denominator `(wcap + dt·gah)` and in
the boundary-flux ledger terms (`bf%cas_enth_in/out`, `whole_enth_out`), so a per-stage `gah` must be
the *same* value in the state update and in the ledger or conservation breaks — the "one flux, both
sides" discipline the ARK already relies on. And the `ustar` floor (0.1) is currently an accidental
stabiliser at very light wind (§1g(iii)); once the feedback is live, check that removing the lag does
not simply move the chattering onto the floor.

**N1 — Under-relax the frozen coefficients. DEMOTED.** `c_{n+1} = (1−θ)·c_n + θ·F(y_n)` gives
`Phi' = a + (1−a)θL`. With the measured `L ≈ −62` at 900 s this needs **`θ < 0.07`** — a coefficient
time constant of ~15 steps, which would lag sunrise and rain onset by hours and is nearly certain to
fail §4 criterion 3. At 150 s (`L ≈ −38`, `a = 0.91`) `θ ≈ 0.5` would suffice, so N1 remains viable as
a *cheap safety net at the current default* if N2 is delayed — but it cannot buy a longer `dt_fast`,
which was its entire attraction. Prediction 3 of §1f is withdrawn.

**N3 — Rosenbrock / W-method (ROS3P, RODAS3).** Linearly implicit: one Jacobian per step, linear
solves only, no nonlinear iteration. The relevant property is that **W-methods tolerate an inexact
Jacobian by design** — which legitimises carrying an approximate `∂c/∂y` (the exact thing N2 adds by
finite differences) without harming order. L-stable and stiffly accurate variants exist at order 3.
This is the best *formal* fit to the problem as §1 characterises it, and the natural landing place if
N2 works but its Newton cost is unattractive. Larger change than N2; do not start here.

**N4 — A self-consistent algebraic (index-1 DAE) canopy air.** The correct form of "make it
diagnostic": fast variables (CAS twins, leaf and wood temperature, ground skin) solved as an algebraic
block with coefficients evaluated at the same state, slow variables (soil heat, soil water, plant
water) integrated. This is CLM5's structure. It is N2 taken to completion, and it should be reached
*through* N2 rather than as a separate build. Two standing requirements: the CAS storage term appears
in the ledgers (`wcap*enth0`/`wcap*enth1` and the `cas_*` budgets in `column_fast_step_ark`), so those
must be re-derived rather than patched; and `wcap` stops being negligible at night and in calm
conditions, so validation must happen there, not at midday.

**N5 — An adaptive freeze cadence, with a real error estimator.** `tau` spans ~80 s (windy midday) to
~2 h (calm night), so the stability limit varies by more than an order of magnitude within a day, and
`dt_fast = 150 s` is set by the worst hour and paid for all 24. Per §1e the controller cannot find this
because it has no estimator for the freeze — so build one: the next step's pre-pass computes
`c(y_{n+1})` anyway, making the drift `||c(y_{n+1}) − c(y_n)||` (or the induced `a`, from `wcap` and
`gah`, both already on `fro%surf`) free. Choose `n_fast_per_slow` per slow step from the previous
step's indicator, keeping `dt_slow/dt_fast` an exact integer. Long steps at night and in calm periods,
short in windy sun.

This is the first proposed exception to the layer-3 rule that `dt_fast` is user-set, and it is a
narrow one: the *cadence* adapts, the *scheme* does not, and the quantity controlled is a measured
error indicator rather than a heuristic. It may well be the largest production win available, since it
recovers the cost of the worst hour across the other 23. Gated on T2 supplying an indicator that
actually tracks `Phi'`.

**N6 — Refuted or inverted; do not rebuild.**

- *A true IMEX tableau built for its own sake.* The general point stands: moving a term from *frozen*
  to *explicit* costs **more** evaluations per step, so a true IMEX can never cheapen the current step
  — its only route to a win is through `dt_fast`. **But note that N2a is the useful part of that idea,
  now targeted.** The first draft's A1/A2 proposed moving the *whole* surface coupling explicit, at a
  cost that scales with cohort count (leaf gas exchange, 89% of the pre-pass at 30 cohorts); §1g says
  only `gah`/`gaw`/`gac` need to move, at 2%. The refuted item is the *scope*, not the mechanism.
- *Partial refresh of the surface coupling per sub-step, as previously scoped* — the earlier
  refutation was measured with `ustar` pinned, which does not pin `gah` (§1g(ii)); it did not test
  what it claimed to test. What is refuted is refreshing the **expensive** parts.
- *Naive `wcap → 0`* — §1d: maximally unstable.
- *Coefficient extrapolation* (`2c_n − c_{n−1}`) — cancels the lag to first order but amplifies an
  already-oscillating state. Under-relaxation is the safe direction; the asymmetry is real. Revisit
  only after N1/N2 have restored contraction.
- *Inner sub-stepping as a stability measure* — §1d: slightly harmful.

## 6. Cheaper steps

Ordered by expected value once T1 says where the time goes:

1. **Allocation churn in the state combinators.** Six routines allocate four per-cohort arrays each on
   entry; one adaptive sub-step calls several. A pre-allocated workspace carried through the march (or
   fixed capacity sized to the site-wide max cohort count, as the driver's scratch already is) removes
   it. Also a §7 prerequisite — concurrent threads hammering the allocator is a classic scaling wall.
2. **`bflux_zero`'s per-call allocation** of the tissue integrals — same treatment.
3. **The RK45 rescue snapshot** copies whole `patch_biophys_t` / `column_budget_t` values per sub-step
   (`meds_fast_step.f90:101`) whether or not a rescue fires. Baseline tier only; measure first.
4. **A cost ceiling wired to a diagnostic, not a control-flow switch.** If a `dt_fast` needs more than
   *N* sub-steps, report it through the work counters and halt under L2. The machinery exists
   (`RK45_WORK_CAP`, the `dt_fast/64` floor); the reporting path does not.
5. **Delete the split residue** (§2c): the `picard_*` surface and its always-zero counters; give
   `ark_niter` real meaning (wire it to `NEWT_MAX`) or make it an honest boolean; resolve `ark_relax`
   (N1 may claim it). These make `[fast]` read as if it offers cost knobs it does not.

## 7. Multi-core

The only lever that costs no accuracy — **provided output is bit-identical regardless of thread
count**. That is a hard requirement, not a nicety: verification here is byte-for-byte netCDF
comparison (`MEDS_NUMERICS_SCOPING.md` §11.2), and a non-deterministic reduction destroys it.

**C0 — the axis.** Patches within a site: columns are independent within a `dt_fast`, coupled only
through the slow loop. Cohorts within a patch are too fine (and the batch kernels are the natural
vectorisation target); sub-steps are sequential; polygons/sites are the real long-run scaling axis but
the runtime is single-polygon today and that is out of scope. Speedup is bounded by `n_patch` and by
**load imbalance** — each patch runs its own adaptive march — so `schedule(dynamic, 1)`.

**C1 — hoist the met sampling out of the patch loop.** `t_sub` depends only on `isub`, yet
`met_advance` (file reader, `intent(inout)`, may reload a netCDF bracket), `met_instant` (solar
geometry, shortwave disaggregation, Weiss–Norman) and `apply_met_to_ctx` all run inside the patch
loop. Precompute the `n_fast_per_slow` samples once per slow step and index them. Correctness-neutral,
independently a win, **and it removes the file reader from the region about to become parallel.**
First, on its own commit, verified byte-identical.

**C2 — per-thread scratch.** BB1 phase 1 hoisted `coh`/`bio`/`aero`/`forc`/`budg`/`gpp_coh…` out of the
loop and sized them to the site-wide max cohort count, to cut `O(n_patch)` allocations per slow step.
That optimisation and threading are in direct conflict — those buffers are the loop-carried shared
state. Keep the hoist but make it per thread (a pool of `n_thread` copies, indexed by thread id), and
verify the allocation count does not regress; that was BB1's whole point.

**C3 — deterministic reductions.** Classify every write:

- *Disjoint, safe as-is:* `site%patch%{cas,soil_e,soil_w,snow}(ip)`, `site%patch%xi_accum(ip)`, and the
  per-cohort slices `site%cohort%*_accum(i0:i0+ncoh-1)`.
- *Reductions:* `site%et_accum`, `site%pheno_tair_sum/n`, the nine `site%work_*`, and local
  `we`/`ww`/`nfail`. OpenMP `reduction(+:…)` sums in **thread-arrival order**, so the answer moves with
  thread count. Accumulate into per-patch arrays and sum them **in patch order** after the loop.
- *Output staging:* `mgr%fast(isub)%…` and `mgr%fast_soil_*` are area-weighted over patches — same
  treatment. `mgr%fast_coh_*` is written by global cohort slot and is disjoint.
- *Not thread-safe:* `write_fast_probe` holds `save`d `unit`/`opened` and writes a shared file. Either
  serialize it or make `fast_probe` and threading mutually exclusive at config time, and say why.

**C4 — build wiring.** `-qopenmp` (IntelLLVM), `-fopenmp` (GNU), `-mp` (NVHPC multicore), applied to
the target owning `meds_fast_dynamics`, threading **off by default** so no existing result moves
without opt-in. Expose it once — a `[run].n_threads` key *or* the environment variable, not both
silently.

**C5 — validation.** Byte-identical diagnostic netCDF at 1/2/4/8 threads; allocation count not
regressed against BB1; `ctest` green on all three back ends; speedup reported against the measured
parallel fraction from T1, not against wall-clock alone.

## 8. Deliberately not in this plan

- **A higher-order tableau.** The freeze dominates the error at production `dt_fast`; an ARK4(3) over a
  first-order semi-discretisation buys nothing. Revisit only if §5 makes the tableau binding.
- **Folding plant hydraulics into the tableau.** Measured 4×–170× *worse*, and a category error: ψ is
  advanced by an exact matrix exponential and the coupled subsystem is ψ-independent within a step. ED2
  makes the same choice for the same reason.
- **Folding soil water into the ARK tableau.** Real work — it would retire the structurally-zero
  dilution in the ARK's error norm — but it buys correctness, not speed. Track separately.
- **GPU offload of the fast loop.** Gated on the bare-array kernel conversion
  (`MEDS_NUMERICS_SCOPING.md` §11). Note the ordering: §7's per-thread scratch and deterministic
  reductions are prerequisites for offload too, so C is not wasted if the GPU path resumes.
- **MPI / multi-polygon parallelism.** The right long-run scaling axis, explicitly not this document's.
- **Re-litigating `split`.** Archived with a tombstone; not evidence for anything here.

## 9. Where the code is

| piece | file |
|---|---|
| dispatch + the RK45→ARK stiff rescue | `src/driver/meds_fast_step.f90` (`column_fast_step`) |
| ESDIRK2 scheme (config name `ark`) | `src/driver/meds_fast_ark.f90` (`ark2_column_step`, `adaptive_ark_march`) |
| the frozen pre-pass — §5 acts here | `src/driver/meds_fast_ark.f90` (`column_prepass`, `build_column_frozen`) |
| the per-stage implicit surface solve — N2 acts here | `src/driver/meds_fast_ark.f90` (`column_be_stage`, `newton_surface_solve`, `jac_surface`) |
| the shared right-hand side | `src/driver/meds_fast_time_derivs.f90` (`surface_derivs`, `column_derivs`) |
| adaptive Cash–Karp scheme (baseline tier) | `src/driver/meds_fast_rk45.f90` |
| tolerances, error norm, step controller | `src/driver/meds_fast_control.f90` |
| the patch loop — §7 acts here | `src/driver/meds_fast_dynamics.f90` (`fast_dynamics`) |
| cadence owner + the `dt_fast` warning | `src/driver/meds_stepper.f90`, `src/shared/config/meds_config.f90` |
| the per-stage conductance refresh (N2a, shipped) | `src/driver/meds_fast_time_derivs.f90` (`refresh_cas_conductances`, `cas_conductances`); armed in `build_column_frozen`, disarmed for the Newton in `column_be_stage` |
| the ground↔CAS conductance — **not** refreshed, and §1i.1 says it need not be | `src/biophysics/meds_canopy_aerodynamics.f90` (`ggbare`/`ggveg`/`ggnet`), consumed via `surface_frozen_t%ggnet` |
| the leaf-water store — N2b acts here | `src/driver/meds_fast_ark.f90` (`advance_water_mass_full`), `column_state_t%leaf_water_mass` |
| benchmark + scoring harness | `scripts/numerics_sweep.py`, `scripts/parity_fidelity.py` |

**Measurement harnesses.** The `Phi'`/Jacobian/decomposition probes are *not* in the repository — they
are standalone programs linked against the built libraries (`build-ifx/*.a` + `build-ifx/modules`)
with **zero source instrumentation**, which is what makes them cheap to write and safe to throw away.
That works because `build_column_frozen`, `adaptive_ark_march`, `column_prepass`, `column_be_stage`
and `aero_bottom_to_top` are all `public`. Rebuild them per session; the recipe is one `ifx` line
against the library list. Probes used: `meas_fastloop{2,3,4}.f90` (§1g), `meas_t3t5.f90` (§1h),
`meas_gpp2.f90` (the carbon attribution), `meas_ground.f90` + `meas_checklist.f90` (§1i).

> **Harness trap, cost me a fake result.** The driver reads **`ccfg%aero`**, not `cfg%aero`. Perturbing
> `cfg%aero%cs_dense` by 25× produced *exactly* zero response, which reads as "`ggnet` is not wired"
> and was really "I perturbed the wrong config object". **An exact zero from a sensitivity probe is a
> harness bug until proven otherwise** — a genuine insensitivity is small, not identically zero.
